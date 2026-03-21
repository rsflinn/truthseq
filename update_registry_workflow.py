#!/usr/bin/env python3
"""
Workflow script for automated registry updates.

Called by GitHub Actions. Reads the existing registry, searches GEO and
ArrayExpress for each disease term, deduplicates, and writes the result.
Everything happens in one process with one read and one write to avoid
intermediate file corruption or quoting inconsistencies.

Usage:
    python3 update_registry_workflow.py
"""

import sys
import os
import csv
import pandas as pd

# Add script directory to path so we can import dataset_search
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dataset_search import search_geo, search_arrayexpress, REGISTRY_PATH

# Terms too vague to produce useful GEO results
SKIP_TERMS = {'healthy', 'developmental', 'cancer'}

# Minimum sample count for registry inclusion
MIN_SAMPLES = 10


def main():
    # --- Step 1: Read the existing registry exactly once ---
    print("=== Reading existing registry ===")
    existing = pd.read_csv(REGISTRY_PATH)
    print(f"  Existing entries: {len(existing)}")

    # Normalize identifiers
    existing['dataset_id'] = existing['dataset_id'].astype(str).str.strip()
    existing['accession'] = existing['accession'].astype(str).str.strip()

    # Remove any pre-existing duplicates (by dataset_id, which is unique per entry)
    before = len(existing)
    existing = existing.drop_duplicates(subset='dataset_id', keep='first')
    if len(existing) < before:
        print(f"  Cleaned {before - len(existing)} pre-existing duplicates")

    # Build lookup sets
    known_dataset_ids = set(existing['dataset_id'].tolist())
    known_accessions = set(existing['accession'].tolist())

    # --- Step 2: Extract disease terms from registry ---
    diseases = set()
    for d in existing['disease'].dropna().unique():
        for term in d.split(';'):
            term = term.strip()
            if term and term not in SKIP_TERMS:
                diseases.add(term)

    # Always include a general perturbation search
    search_queries = [f"{d} differential expression" for d in sorted(diseases)]
    search_queries.append("CRISPR perturbation screen gene expression")

    print(f"\n=== Searching {len(search_queries)} queries ===")
    for q in search_queries:
        print(f"  - {q}")

    # --- Step 3: Run all searches, collect results ---
    all_new = []
    seen_in_this_run = set()  # track dataset_ids added during this run

    for query in search_queries:
        print(f"\n--- {query} ---")

        # Search GEO
        try:
            geo_results = search_geo(query, max_results=10, human_only=True)
            print(f"  GEO: {len(geo_results)} results")
        except Exception as e:
            print(f"  GEO error: {e}")
            geo_results = []

        # Search ArrayExpress
        try:
            ae_results = search_arrayexpress(query, max_results=5)
            print(f"  ArrayExpress: {len(ae_results)} results")
        except Exception as e:
            print(f"  ArrayExpress error: {e}")
            ae_results = []

        for r in geo_results + ae_results:
            did = str(r.get('dataset_id', '')).strip()
            acc = str(r.get('accession', '')).strip()

            # Skip if already in existing registry
            if did and did in known_dataset_ids:
                continue
            if acc and acc in known_accessions:
                continue

            # Skip if already added in this run
            if did and did in seen_in_this_run:
                continue

            # Skip small datasets
            try:
                n = int(str(r.get('n_samples', '0')).replace('~', '').replace(',', '').replace('>', ''))
            except (ValueError, TypeError):
                n = 0
            if 0 < n < MIN_SAMPLES:
                continue

            all_new.append({
                'dataset_id': did,
                'disease': r.get('disease', ''),
                'data_type': r.get('data_type', ''),
                'source': r.get('source', ''),
                'accession': acc,
                'description': r.get('description', ''),
                'n_samples': r.get('n_samples', ''),
                'species': r.get('species', ''),
                'tissue': r.get('tissue', ''),
                'cell_types': '',
                'access_type': r.get('access_type', 'open'),
                'download_url': r.get('download_url', ''),
                'download_instructions': r.get('download_instructions', ''),
                'paper_doi': r.get('paper_doi', ''),
                'paper_citation': r.get('paper_citation', ''),
                'last_verified': pd.Timestamp.now().strftime('%Y-%m-%d'),
            })

            if did:
                seen_in_this_run.add(did)
            if acc:
                known_accessions.add(acc)

    # --- Step 4: Merge and write ---
    print(f"\n=== Results ===")
    print(f"  New datasets found: {len(all_new)}")

    if all_new:
        new_df = pd.DataFrame(all_new)
        combined = pd.concat([existing, new_df], ignore_index=True)
    else:
        combined = existing

    # Final safety dedup on dataset_id
    before_final = len(combined)
    combined['dataset_id'] = combined['dataset_id'].astype(str).str.strip()
    combined = combined.drop_duplicates(subset='dataset_id', keep='first')
    if len(combined) < before_final:
        print(f"  Final dedup removed {before_final - len(combined)} entries")

    print(f"  Total entries: {len(combined)}")

    # Write with consistent quoting (quote fields containing commas)
    combined.to_csv(REGISTRY_PATH, index=False, quoting=csv.QUOTE_NONNUMERIC)
    print(f"  Written to {REGISTRY_PATH}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
