#!/usr/bin/env python3
"""Stream public processed matrices; retain five CpGs (optionally fixed 353).

This is acquisition/preparation only, not a completed validation analysis.
Output belongs in a temporary cache, not the retained results directories.
No normalization, clock fitting, sample selection or methylome-wide analysis.
Dependencies: Python >=3.9, requests (only for URLs).
"""
import os
from pathlib import Path
import argparse
import contextlib
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlparse

LOCI = ('cg07158339', 'cg22947000', 'cg01511567', 'cg19761273', 'cg13460409')
PROJECT = Path(os.environ['AGE_WORKSPACE'])
SOURCES = {
    'GSE74193_beta': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74193/suppl/GSE74193_GEO_procData.csv.gz',
    'GSE74193_metadata': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74193/matrix/GSE74193_series_matrix.txt.gz',
    'GSE60274': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60274/matrix/GSE60274_series_matrix.txt.gz',
    'GSE195684_beta': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/GSE195684/suppl/GSE195684_Matrix_Processed_BADY_GBMnordic.csv.gz',
    'GSE195684_metadata': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/GSE195684/matrix/GSE195684_series_matrix.txt.gz',
    'TCGA_GBM_beta': 'https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/HumanMethylation450.gz',
}


@contextlib.contextmanager
def source_text(source):
    with contextlib.ExitStack() as stack:
        if source.startswith(('https://', 'http://')):
            import requests
            response = stack.enter_context(requests.get(source, stream=True, timeout=(20, 90)))
            response.raise_for_status()
            raw = response.raw
            raw.decode_content = True
        else:
            raw = stack.enter_context(open(source, 'rb'))
        if urlparse(source).path.endswith('.gz'):
            raw = stack.enter_context(gzip.GzipFile(fileobj=raw))
        text = stack.enter_context(io.TextIOWrapper(raw, encoding='utf-8-sig'))
        yield text


def extract(source, targets, delimiter=None, metadata_only=False):
    selected, metadata, header = {}, [], None
    with source_text(source) as stream:
        for line in stream:
            if line.startswith('!'):
                if line.startswith(('!Sample_', '!Series_')):
                    metadata.append(next(csv.reader([line], delimiter='\t')))
                if metadata_only and line.startswith('!series_matrix_table_begin'):
                    break
                continue
            if line.startswith('#') or not line.strip():
                continue
            if metadata_only:
                continue
            if delimiter is None:
                delimiter = '\t' if '\t' in line else ','
            row = next(csv.reader([line], delimiter=delimiter))
            if header is None:
                header = row
                continue
            probe = row[0].strip()
            if probe not in targets:
                continue
            if probe in selected:
                raise ValueError(f'Duplicate target row: {probe}; inspect source before proceeding')
            if len(row) != len(header):
                raise ValueError(f'Column count differs for {probe}')
            selected[probe] = row
            if len(selected) == len(targets):
                break
    if metadata_only:
        if not metadata:
            raise ValueError('No GEO series/sample metadata found')
    else:
        missing = set(LOCI).difference(selected)
        if missing:
            raise ValueError(f'Five-locus coverage not satisfied: missing {sorted(missing)}')
    return header, selected, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', help='Listed source name, local processed matrix, or public URL')
    parser.add_argument('--cache', type=Path, required=True, help='Temporary output directory')
    parser.add_argument('--name', required=True, help='Cache filename stem')
    parser.add_argument('--metadata-only', action='store_true')
    parser.add_argument('--clock353', action='store_true', help='Retain existing clock coefficients\' 353 probes; no training')
    parser.add_argument('--delimiter', choices=('tab', 'comma'))
    args = parser.parse_args()
    if Path(args.name).name != args.name or args.name in ('.', '..'):
        parser.error('--name must be a filename stem')
    retained = Path(__file__).resolve().parents[1]
    cache = args.cache.resolve()
    if cache == retained or retained in cache.parents:
        parser.error('Use a temporary cache outside 09_EXTERNAL_VALIDATION')
    if PROJECT.resolve() == cache or PROJECT.resolve() in cache.parents:
        parser.error('Use /tmp or another cache outside the validated project')
    targets = set(LOCI)
    if args.clock353:
        with open(PROJECT / 'data/reference/horvath2013/clock353_hg38.tsv') as f:
            targets = {r['CpGmarker'] for r in csv.DictReader(f, delimiter='\t')}
        assert len(targets) == 353 and set(LOCI) <= targets
    source = SOURCES.get(args.source, args.source)
    delimiter = {'tab': '\t', 'comma': ','}.get(args.delimiter)
    header, rows, metadata = extract(source, targets, delimiter, args.metadata_only)
    cache.mkdir(parents=True, exist_ok=True)
    output = cache / (args.name + '.targeted.tsv')
    if output.exists() or (cache / (args.name + '.metadata.json')).exists():
        raise FileExistsError('Refusing to overwrite a previous extraction')
    if not args.metadata_only:
        with output.open('w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(header)
            writer.writerows(rows[p] for p in sorted(rows))
    audit = {
        'source': source, 'requested_probes': sorted(targets),
        'retained_probes': sorted(rows), 'missing_probes': sorted(targets.difference(rows)),
        'sample_selection': 'NOT_YET_CURATED', 'analysis_status': 'NOT_ANALYSED',
        'metadata': metadata,
        'targeted_output_sha256': hashlib.sha256(output.read_bytes()).hexdigest() if not args.metadata_only else None,
        'full_source_checksum': None,
        'note': 'Stream can stop after all requested rows; hash refers only to retained rows. '
                'Probe coverage, beta scale, sample identities, diagnoses, ages, repeated donors, '
                'processing, IDH/G-CIMP and overlap across cohorts require review before analysis.'
    }
    (cache / (args.name + '.metadata.json')).write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps({k: v for k, v in audit.items() if k != 'metadata'}, indent=2))


if __name__ == '__main__':
    main()
