#!/usr/bin/env python3
"""
Workflow script for automated registry updates.

Called by GitHub Actions. Searches GEO and ArrayExpress for new datasets,
checks each against the existing registry, and APPENDS only new entries.
Never rewrites or reformats existing entries.

Usage:
    python3 update_registry_workflow.py
"""

import sys
import os

# Add script directory to path so we can import dataset_search
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dataset_search import search_geo, search_arrayexpress
from disease_lookup import check_dataset_relevance

REGISTRY_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dataset_registry.csv")

# Terms too vague to produce useful GEO results
SKIP_TERMS = {'healthy', 'developmental', 'cancer'}

# Minimum sample count for registry inclusion
MIN_SAMPLES = 10


def load_existing_ids(path):
    """Read the registry and return sets of known dataset_ids and accessions.
    Uses plain file reading to avoid pandas reformatting anything."""
    dataset_ids = set()
    accessions = set()

    if not os.path.exists(path):
        return dataset_ids, accessions

    with open(path, 'r') as f:
        header = f.readline().strip().split(',')
        # Find column positions
        try:
            id_col = header.index('dataset_id')
            acc_col = header.index('accession')
        except ValueError:
            print("ERROR: registry CSV missing dataset_id or accession column")
            return dataset_ids, accessions

        for line in f:
            # Simple CSV parsing: split by comma, but respect quotes
            fields = []
            current = ''
            in_quotes = False
            for ch in line.strip():
                if ch == '"':
                    in_quotes = not in_quotes
                elif ch == ',' and not in_quotes:
                    fields.append(current.strip().strip('"'))
                    current = ''
                else:
                    current += ch
            fields.append(current.strip().strip('"'))

            if len(fields) > max(id_col, acc_col):
                did = fields[id_col].strip()
                acc = fields[acc_col].strip()
                if did:
                    dataset_ids.add(did)
                if acc:
                    accessions.add(acc)

    return dataset_ids, accessions


def format_csv_field(value):
    """Quote a field if it contains commas, quotes, or newlines."""
    s = str(value)
    if ',' in s or '"' in s or '\n' in s:
        return '"' + s.replace('"', '""') + '"'
    return s


def format_row(entry):
    """Format a dict as a CSV row string matching the registry columns."""
    cols = [
        'dataset_id', 'disease', 'data_type', 'source', 'accession',
        'description', 'n_samples', 'species', 'tissue', 'cell_types',
        'access_type', 'download_url', 'download_instructions',
        'paper_doi', 'paper_citation', 'last_verified'
    ]
    return ','.join(format_csv_field(entry.get(c, '')) for c in cols)


def main():
    # --- Step 1: Load existing IDs (plain file read, no pandas) ---
    print("=== Reading existing registry ===")
    known_ids, known_accs = load_existing_ids(REGISTRY_PATH)
    print(f"  Known dataset_ids: {len(known_ids)}")
    print(f"  Known accessions: {len(known_accs)}")

    # --- Step 2: Extract disease terms from registry ---
    diseases = set()
    with open(REGISTRY_PATH, 'r') as f:
        header = f.readline().strip().split(',')
        try:
            disease_col = header.index('disease')
        except ValueError:
            disease_col = 1  # fallback
        for line in f:
            fields = []
            current = ''
            in_quotes = False
            for ch in line.strip():
                if ch == '"':
                    in_quotes = not in_quotes
                elif ch == ',' and not in_quotes:
                    fields.append(current.strip().strip('"'))
                    current = ''
                else:
                    current += ch
            fields.append(current.strip().strip('"'))

            if len(fields) > disease_col:
                for term in fields[disease_col].split(';'):
                    term = term.strip()
                    if term and term not in SKIP_TERMS:
                        diseases.add(term)

    search_queries = [f"{d} differential expression" for d in sorted(diseases)]
    search_queries.append("CRISPR perturbation screen gene expression")

    print(f"\n=== Searching {len(search_queries)} queries ===")
    for q in search_queries:
        print(f"  - {q}")

    # --- Step 3: Search and collect truly new entries ---
    new_entries = []
    seen_this_run = set()

    from datetime import datetime

    for query in search_queries:
        print(f"\n--- {query} ---")

        # Extract the disease term from the query for relevance checking
        # Queries look like "bipolar differential expression" or "CRISPR perturbation screen..."
        query_disease = query.replace('differential expression', '').replace('gene expression', '').strip()

        try:
            geo_results = search_geo(query, max_results=10, human_only=True)
            print(f"  GEO: {len(geo_results)} results")
        except Exception as e:
            print(f"  GEO error: {e}")
            geo_results = []

        try:
            ae_results = search_arrayexpress(query, max_results=5)
            print(f"  ArrayExpress: {len(ae_results)} results")
        except Exception as e:
            print(f"  ArrayExpress error: {e}")
            ae_results = []

        for r in geo_results + ae_results:
            did = str(r.get('dataset_id', '')).strip()
            acc = str(r.get('accession', '')).strip()

            # Skip if already in registry
            if did and did in known_ids:
                continue
            if acc and acc in known_accs:
                continue
            # Skip if already found in this run
            if did and did in seen_this_run:
                continue

            # Skip small datasets
            try:
                n = int(str(r.get('n_samples', '0')).replace('~', '').replace(',', '').replace('>', ''))
            except (ValueError, TypeError):
                n = 0
            if 0 < n < MIN_SAMPLES:
                continue

            # Relevance gate: check if the description matches the search disease
            description = r.get('description', '')
            relevance = check_dataset_relevance(description, query_disease)
            disease_tag = r.get('disease', '')

            if relevance['confidence'] == 'mismatch':
                # Tag the disease field to flag the mismatch
                actual_subject = ', '.join(relevance['detected_diseases'])
                disease_tag = f"FLAGGED:{disease_tag} (description suggests: {actual_subject})"
                print(f"  FLAGGED: {acc} — tagged as {query_disease} but description suggests {actual_subject}")

            entry = {
                'dataset_id': did,
                'disease': disease_tag,
                'data_type': r.get('data_type', ''),
                'source': r.get('source', ''),
                'accession': acc,
                'description': description,
                'n_samples': r.get('n_samples', ''),
                'species': r.get('species', ''),
                'tissue': r.get('tissue', ''),
                'cell_types': '',
                'access_type': r.get('access_type', 'open'),
                'download_url': r.get('download_url', ''),
                'download_instructions': r.get('download_instructions', ''),
                'paper_doi': r.get('paper_doi', ''),
                'paper_citation': r.get('paper_citation', ''),
                'last_verified': datetime.now().strftime('%Y-%m-%d'),
            }
            new_entries.append(entry)

            if did:
                seen_this_run.add(did)
                known_ids.add(did)
            if acc:
                known_accs.add(acc)

    # --- Step 4: Append only new entries (never touch existing lines) ---
    print(f"\n=== Results ===")
    print(f"  New datasets found: {len(new_entries)}")

    if new_entries:
        with open(REGISTRY_PATH, 'a') as f:
            for entry in new_entries:
                f.write(format_row(entry) + '\n')
        print(f"  Appended {len(new_entries)} new rows to registry")
    else:
        print("  No new datasets to add")

    return 0


if __name__ == '__main__':
    sys.exit(main())
