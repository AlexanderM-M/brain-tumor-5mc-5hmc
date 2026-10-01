"""Validate CNA completion without promoting rejected fits to accepted values."""
import os
from pathlib import Path
import csv, hashlib, json, math, re, shutil
from pathlib import Path

def nonempty(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f'Missing or empty CNA output: {path}')
    return path

def table(path, required):
    nonempty(path)
    with path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not set(required).issubset(reader.fieldnames or []):
            raise RuntimeError(f'Invalid CNA table header: {path}')
        n = 0
        for row in reader:
            n += 1
            if None in row or any(row[k] is None for k in required):
                raise RuntimeError(f'Malformed CNA row {n}: {path}')
            yield row
        if n == 0:
            raise RuntimeError(f'No data rows: {path}')

def finite(value):
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f'Non-finite CNA value: {value}')
    return result

def validate_cna(label, folder, log, exit_code):
    if exit_code != 0:
        raise RuntimeError(f'{label}: CNA command exited {exit_code}; see {log}')
    text = nonempty(log).read_text()
    if 'Total time to perform copy number calling' not in text:
        raise RuntimeError(f'{label}: missing CNA completion footer')
    if 'Traceback (most recent call last)' in text:
        raise RuntimeError(f'{label}: traceback in CNA log')
    rel = folder / f'{label}_read_counts_self_log2r_segmented.tsv'
    raw = nonempty(folder / f'{label}_raw_read_counts.tsv')
    alleles = nonempty(folder / f'{label}_allele_counts_hetSNPs.bed')
    segments = {}; chromosomes = set(); last = {}; bins = 0
    for r in table(rel, ['chromosome', 'start', 'end', 'log2r_copynumber', 'seg_id', 'seg_log2r_copynumber']):
        ch = r['chromosome']; start = int(r['start']); end = int(r['end'])
        if start < 0 or end < start or start <= last.get(ch, -1):
            raise RuntimeError(f'{label}: invalid/overlapping relative CNV coordinates')
        last[ch] = end; chromosomes.add(ch); bins += 1
        finite(r['log2r_copynumber']); value = finite(r['seg_log2r_copynumber'])
        key = (ch, r['seg_id'])
        if key not in segments:
            segments[key] = [ch, start, end, r['seg_id'], 0, value]
        s = segments[key]
        if s[5] != value:
            raise RuntimeError(f'{label}: inconsistent value within segment')
        s[2] = end; s[4] += 1
    if chromosomes != {f'chr{i}' for i in range(1,23)}:
        raise RuntimeError(f'{label}: incomplete autosomal segmentation')
    fit = folder / f'{label}_fitted_purity_ploidy.tsv'
    absolute = folder / f'{label}_segmented_absolute_copy_number.tsv'
    params = folder / 'No_fit_found_PARAMS.tsv'
    no_fit = 'WARNING: NO FITS FOUND.' in text
    outputs = [raw, alleles, rel]
    result = dict(sample_id=label, command_exit_code=exit_code, cnv_output_valid=True,
                  cnv_validation='structural integrity; scientific interpretation awaits QC',
                  relative_cnv_bins=bins, relative_cnv_segments=len(segments),
                  autosomes=len(chromosomes), purity=None, ploidy=None,
                  preliminary_cellularity_is_final_purity=False)
    if no_fit:
        nonempty(params)
        if fit.exists() or absolute.exists():
            raise RuntimeError(f'{label}: contradictory no-fit and accepted-fit outputs')
        if 'No viable solution' not in params.read_text():
            raise RuntimeError(f'{label}: invalid no-fit diagnostic')
        result.update(outcome='COMPLETED_NO_ACCEPTED_FIT', purity_ploidy_status='unresolved_no_accepted_fit', absolute_copy_number_available=False)
        outputs.append(params)
    else:
        rows = list(table(fit, ['purity', 'ploidy']))
        if len(rows) != 1:
            raise RuntimeError(f'{label}: expected one accepted fit')
        purity, ploidy = finite(rows[0]['purity']), finite(rows[0]['ploidy'])
        if not 0 < purity <= 1 or ploidy <= 0:
            raise RuntimeError(f'{label}: invalid accepted fit')
        for r in table(absolute, ['chromosome', 'start', 'end', 'copyNumber']):
            if int(r['end']) < int(r['start']) or finite(r['copyNumber']) < 0:
                raise RuntimeError(f'{label}: invalid absolute CNV')
        result.update(outcome='COMPLETED_ACCEPTED_FIT', purity_ploidy_status='accepted_fit_awaiting_qc', purity=purity, ploidy=ploidy, absolute_copy_number_available=True)
        outputs.extend([fit, absolute])
    result['outputs'] = [str(p) for p in outputs]
    result['output_evidence'] = [dict(path=str(p), bytes=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns) for p in outputs]
    result['relative_segmentation_sha256'] = hashlib.sha256(rel.read_bytes()).hexdigest()
    result['preliminary_cellularity_estimates'] = [float(v) for v in re.findall(r'estimated cellularity using hetSNPs = ([0-9]+(?:\.[0-9]+)?)', text)]
    return result, list(segments.values())

def save_cna_outcome(result, segments, metadata):
    summary = metadata / 'unreviewed_summaries'; summary.mkdir(exist_ok=True)
    label = result['sample_id']
    with (summary / f'{label}_relative_copy_number_segments.tsv').open('w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['chromosome','start','end','segment_id','bin_count','seg_log2r_copynumber'])
        writer.writerows(segments)
    preliminary = metadata / 'preliminary_estimates'; preliminary.mkdir(exist_ok=True)
    (preliminary / 'cellularity.json').write_text(json.dumps(dict(values=result['preliminary_cellularity_estimates'], method='hetSNPs', accepted_final_purity=False, source_log=str(metadata/'copy_number_purity_ploidy.workflow.log')), indent=2)+'\n')
    if result['outcome'] == 'COMPLETED_NO_ACCEPTED_FIT':
        rejected = metadata / 'rejected_fit_diagnostics'; rejected.mkdir(exist_ok=True)
        params = next(Path(p) for p in result['outputs'] if Path(p).name == 'No_fit_found_PARAMS.tsv')
        shutil.copy2(params, rejected / params.name)
    (metadata / 'cna_outcome.json').write_text(json.dumps(result, indent=2)+'\n')
