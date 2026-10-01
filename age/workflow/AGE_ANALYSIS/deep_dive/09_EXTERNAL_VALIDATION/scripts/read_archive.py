#!/usr/bin/env python3
"""Read verified ZIP/GZIP members, retaining only prespecified clock rows.
Run through the AGE workspace launcher; --cache accepts a processed public-cohort cache.
No whole-methylome statistical processing; original archive remains read-only.
"""
import os
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import json
import zipfile
from pathlib import Path
import pandas as pd

PROJECT = Path(os.environ['AGE_WORKSPACE'])
ROOT = PROJECT / 'AGE_ANALYSIS/deep_dive/09_EXTERNAL_VALIDATION'
LOCI = ['cg07158339', 'cg22947000', 'cg01511567', 'cg19761273', 'cg13460409']


class HashReader:
    def __init__(self, raw):
        self.raw = raw
        self.digest = hashlib.sha256()

    def read(self, n=-1):
        data = self.raw.read(n)
        self.digest.update(data)
        return data


def metadata_frame(rows):
    geo = next(r[1:] for r in rows if r[0] == '!Sample_geo_accession')
    records = [{'sample_id': s} for s in geo]
    for r in rows:
        key = r[0].replace('!Sample_', '')
        if not r[0].startswith('!Sample_'):
            continue
        assert len(r[1:]) == len(geo), (key, len(r), len(geo))
        for out, value in zip(records, r[1:]):
            if key.startswith('characteristics') and ': ' in value:
                k, v = value.split(': ', 1)
                if k in out and out[k] != v:
                    raise ValueError(f'Conflicting metadata: {k}')
                out[k] = v
            elif not key.startswith('characteristics'):
                if key in out and out[key] != value:
                    out[key] += ' | ' + value
                else:
                    out[key] = value
    return pd.DataFrame(records).set_index('sample_id')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--archive', type=Path, default=PROJECT / 'Age_Analysis.zip')
    ap.add_argument('--cache', type=Path, default=Path(os.environ['AGE_EXTERNAL_CACHE']))
    args = ap.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    weights = pd.read_csv(PROJECT / 'data/reference/horvath2013/clock353_hg38.tsv', sep='\t')
    clock = set(weights.CpGmarker)
    assert len(clock) == 353 and set(LOCI) <= clock
    # Deployment-specific manuscript archive checks are outside data ingestion.
    hashes = {}
    before = args.archive.stat()
    audit = {'archive_path': str(args.archive), 'members': [], 'canonical_files_verified': len(hashes)}
    with zipfile.ZipFile(args.archive) as z:
        expected = {}
        for line in z.read('Age_Analysis/SHA256SUMS').decode().splitlines():
            checksum, name = line.split(None, 1)
            expected[name.strip()] = checksum
        audit['transfer_manifest'] = json.loads(z.read('Age_Analysis/download_manifest.json'))
        members = [i for i in z.infolist() if not i.is_dir()]
        for info in members:
            base = Path(info.filename).name
            print('READ', base, info.file_size, flush=True)
            if base in ('README_TRANSFER.txt', 'SHA256SUMS', 'download_manifest.json', 'download_urls.txt'):
                content = z.read(info)
                if base in expected:
                    assert hashlib.sha256(content).hexdigest() == expected[base]
                continue
            if base == 'GBM_clinicalMatrix':
                content = z.read(info)
                assert hashlib.sha256(content).hexdigest() == expected[base]
                df = pd.read_csv(io.BytesIO(content), sep='\t', dtype=str)
                df.to_csv(args.cache / 'TCGA_clinical.tsv', sep='\t', index=False)
                audit['members'].append({'file': base, 'zip_crc': info.CRC, 'n_rows': len(df)})
                continue
            if not base.endswith('.gz'):
                raise ValueError(f'Unexpected archive member: {base}')
            targets = set(LOCI) if base.startswith('GSE74193') else clock
            records, header, metadata, n_rows, delim = {}, None, [], 0, None
            digest = hashlib.sha256()
            with z.open(info) as compressed:
                hashed = HashReader(compressed)
                g = gzip.GzipFile(fileobj=hashed)
                with io.TextIOWrapper(g, encoding='utf-8-sig') as text:
                    for line in text:
                        if line.startswith('!'):
                            if line.startswith(('!Sample_', '!Series_')):
                                metadata.append(next(csv.reader([line], delimiter='\t')))
                            continue
                        if not line.strip() or line.startswith('#'):
                            continue
                        if header is None:
                            delim = '\t' if '\t' in line else ','
                            header = next(csv.reader([line], delimiter=delim))
                            continue
                        n_rows += 1
                        probe = line.split(delim, 1)[0].strip('"\r\n')
                        if probe not in targets:
                            continue
                        if probe in records:
                            raise ValueError(f'Duplicate probe {base}: {probe}')
                        row = next(csv.reader([line], delimiter=delim))
                        assert len(row) == len(header), (base, probe)
                        records[probe] = row[1:]
                        digest.update(line.encode())
                assert hashed.digest.hexdigest() == expected[base], ('Transfer SHA256 mismatch', base)
            stem = base.replace('.txt.gz', '').replace('.csv.gz', '').replace('.gz', '')
            if records:
                df = pd.DataFrame.from_dict(records, orient='index', columns=header[1:])
                df.index.name = 'CpGmarker'
                df.to_csv(args.cache / (stem + '.targets.tsv'), sep='\t')
                assert set(LOCI) <= set(records), (base, 'missing required CpGs')
            if metadata:
                metadata_frame(metadata).to_csv(args.cache / (stem + '.metadata.tsv'), sep='\t')
                (args.cache / (stem + '.metadata.json')).write_text(json.dumps(metadata, indent=2))
            audit['members'].append({'file': base, 'zip_crc': info.CRC, 'gzip_and_zip_crc_checked': True, 'transfer_sha256_verified': expected[base],
                                     'matrix_rows_streamed': n_rows, 'target_rows_retained': len(records),
                                     'target_text_sha256': digest.hexdigest(), 'missing_targets': sorted(targets - set(records))})
            print('DONE', base, 'targets', len(records), 'metadata', len(metadata), flush=True)
    after = args.archive.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'Archive changed during reading'
    h = hashlib.sha256()
    with args.archive.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    audit['archive_sha256'] = h.hexdigest()
    audit['archive_bytes'] = after.st_size
    (args.cache / 'ingest_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print('COMPLETE', args.cache, flush=True)


if __name__ == '__main__':
    main()
