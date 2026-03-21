#!/usr/bin/env python3
"""
TruthSeq Dataset Search: Find publicly available expression datasets.
=====================================================================

Searches GEO (NCBI), ArrayExpress (EBI), and the local registry for
datasets matching a disease, tissue, or perturbation query. Returns
ranked results with download instructions.

This is a discovery tool, not a downloader. It tells you what exists
and how to get it. For datasets that require login (Synapse, dbGaP),
it provides the access instructions.

Usage:
    # Search for autism datasets
    python3 dataset_search.py --query "autism spectrum disorder"

    # Search for perturbation data in neurons
    python3 dataset_search.py --query "CRISPR perturbation neurons"

    # Search local registry only (no API calls)
    python3 dataset_search.py --query "autism" --registry-only

    # Update registry with new GEO results
    python3 dataset_search.py --query "autism single-cell RNA-seq brain" --update-registry
"""

import os
import sys
import argparse
import logging
import json
import time
import csv
from datetime import datetime

import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

REGISTRY_PATH = os.path.join(os.path.dirname(__file__), "dataset_registry.csv")


# ============================================================
# Local Registry Search
# ============================================================

def search_registry(query, registry_path=REGISTRY_PATH):
    """Search the local dataset registry for matching datasets."""
    if not os.path.exists(registry_path):
        log.warning(f"Registry not found: {registry_path}")
        return []

    df = pd.read_csv(registry_path)
    query_terms = query.lower().split()

    scores = []
    for idx, row in df.iterrows():
        searchable = ' '.join(str(v).lower() for v in row.values)
        score = sum(1 for term in query_terms if term in searchable)
        if score > 0:
            scores.append((score, idx))

    scores.sort(reverse=True)
    results = []
    for score, idx in scores:
        row = df.iloc[idx]
        results.append({
            'source': 'registry',
            'dataset_id': row.get('dataset_id', ''),
            'accession': row.get('accession', ''),
            'description': row.get('description', ''),
            'disease': row.get('disease', ''),
            'data_type': row.get('data_type', ''),
            'n_samples': row.get('n_samples', ''),
            'species': row.get('species', ''),
            'tissue': row.get('tissue', ''),
            'access_type': row.get('access_type', ''),
            'download_url': row.get('download_url', ''),
            'download_instructions': row.get('download_instructions', ''),
            'paper_doi': row.get('paper_doi', ''),
            'paper_citation': row.get('paper_citation', ''),
            'match_score': score,
        })

    return results


# ============================================================
# GEO Search (NCBI E-utilities)
# ============================================================

def search_geo(query, max_results=20):
    """
    Search NCBI GEO for datasets matching the query.
    Uses E-utilities API (free, no authentication required).
    """
    import requests

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Build a GEO-specific query
    # Filter for expression profiling datasets
    geo_query = f'({query}) AND ("expression profiling by high throughput sequencing"[DataSet Type] OR "expression profiling by array"[DataSet Type])'

    # Step 1: Search for matching GEO DataSets
    log.info(f"Searching GEO for: {query}")
    try:
        search_resp = requests.get(f"{base_url}/esearch.fcgi", params={
            'db': 'gds',
            'term': geo_query,
            'retmax': max_results,
            'retmode': 'json',
            'sort': 'relevance',
        }, timeout=30)
        search_resp.raise_for_status()
        search_data = search_resp.json()

        id_list = search_data.get('esearchresult', {}).get('idlist', [])
        total_count = int(search_data.get('esearchresult', {}).get('count', 0))
        log.info(f"  Found {total_count} GEO datasets ({len(id_list)} retrieved)")

        if not id_list:
            return []

        # Step 2: Fetch summaries for each dataset
        time.sleep(0.5)  # Rate limiting
        summary_resp = requests.get(f"{base_url}/esummary.fcgi", params={
            'db': 'gds',
            'id': ','.join(id_list),
            'retmode': 'json',
        }, timeout=30)
        summary_resp.raise_for_status()
        summary_data = summary_resp.json().get('result', {})

        results = []
        for gds_id in id_list:
            info = summary_data.get(gds_id, {})
            if not info:
                continue

            accession = info.get('accession', '')
            title = info.get('title', '')
            summary = info.get('summary', '')
            n_samples = info.get('n_samples', '')
            gpl = info.get('gpl', '')
            gds_type = info.get('gdstype', '')
            taxon = info.get('taxon', '')

            # Build GEO URL
            if accession.startswith('GSE'):
                geo_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
            elif accession.startswith('GDS'):
                geo_url = f"https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc={accession}"
            else:
                geo_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"

            results.append({
                'source': 'GEO',
                'dataset_id': f"geo_{accession.lower()}",
                'accession': accession,
                'description': f"{title}. {summary[:200]}",
                'disease': query,
                'data_type': gds_type or 'RNA-seq',
                'n_samples': n_samples,
                'species': taxon,
                'tissue': '',
                'access_type': 'open',
                'download_url': geo_url,
                'download_instructions': f"1. Go to {geo_url}. 2. Download supplementary files or series matrix. 3. For processed DE results, check the paper's supplementary materials.",
                'paper_doi': '',
                'paper_citation': '',
                'match_score': 0,
            })

        return results

    except Exception as e:
        log.warning(f"  GEO search failed: {e}")
        return []


# ============================================================
# ArrayExpress Search (EBI)
# ============================================================

def search_arrayexpress(query, max_results=10):
    """Search EMBL-EBI ArrayExpress/BioStudies for datasets."""
    import requests

    log.info(f"Searching ArrayExpress for: {query}")
    try:
        resp = requests.get("https://www.ebi.ac.uk/biostudies/api/v1/search", params={
            'query': query,
            'type': 'study',
            'pageSize': max_results,
        }, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        hits = data.get('hits', [])
        total = data.get('totalHits', 0)
        log.info(f"  Found {total} ArrayExpress results ({len(hits)} retrieved)")

        results = []
        for hit in hits:
            accession = hit.get('accession', '')
            title = hit.get('title', '')
            description = hit.get('content', '') or title

            results.append({
                'source': 'ArrayExpress',
                'dataset_id': f"ae_{accession.lower()}",
                'accession': accession,
                'description': f"{title}. {str(description)[:200]}",
                'disease': query,
                'data_type': 'RNA-seq',
                'n_samples': '',
                'species': '',
                'tissue': '',
                'access_type': 'open',
                'download_url': f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}",
                'download_instructions': f"1. Go to ArrayExpress page for {accession}. 2. Download processed data files.",
                'paper_doi': '',
                'paper_citation': '',
                'match_score': 0,
            })

        return results

    except Exception as e:
        log.warning(f"  ArrayExpress search failed: {e}")
        return []


# ============================================================
# Display and Update
# ============================================================

def display_results(results, verbose=False):
    """Print search results in a readable format."""
    if not results:
        print("\n  No datasets found.\n")
        return

    print(f"\n{'='*70}")
    print(f"  Found {len(results)} datasets")
    print(f"{'='*70}\n")

    for i, r in enumerate(results, 1):
        access_tag = "[OPEN]" if r.get('access_type') == 'open' else "[LOGIN REQUIRED]"
        print(f"  {i}. {r['accession']} {access_tag}")
        print(f"     Source: {r['source']} | Type: {r['data_type']} | Samples: {r['n_samples']}")
        print(f"     {r['description'][:120]}")
        if r.get('paper_citation'):
            print(f"     Paper: {r['paper_citation']}")
        if r.get('download_url'):
            print(f"     URL: {r['download_url']}")
        if verbose and r.get('download_instructions'):
            print(f"     How to get it: {r['download_instructions']}")
        print()


def update_registry(new_results, registry_path=REGISTRY_PATH):
    """Add new results to the registry, skipping duplicates."""
    if not new_results:
        return 0

    existing = pd.read_csv(registry_path) if os.path.exists(registry_path) else pd.DataFrame()
    existing_accessions = set(existing['accession'].tolist()) if 'accession' in existing.columns else set()

    new_rows = []
    for r in new_results:
        if r['accession'] not in existing_accessions:
            new_rows.append({
                'dataset_id': r.get('dataset_id', ''),
                'disease': r.get('disease', ''),
                'data_type': r.get('data_type', ''),
                'source': r.get('source', ''),
                'accession': r.get('accession', ''),
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
                'last_verified': datetime.now().strftime('%Y-%m-%d'),
            })

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        combined = pd.concat([existing, new_df], ignore_index=True)
        combined.to_csv(registry_path, index=False)
        log.info(f"Added {len(new_rows)} new datasets to registry (total: {len(combined)})")
    else:
        log.info("No new datasets to add (all already in registry)")

    return len(new_rows)


# ============================================================
# Main: Combined search
# ============================================================

def find_datasets(query, registry_only=False, max_results=20, verbose=False):
    """
    Search all sources for datasets matching the query.
    Returns combined, deduplicated results.
    """
    all_results = []

    # Always search local registry first
    registry_results = search_registry(query)
    all_results.extend(registry_results)

    if not registry_only:
        # Search GEO
        geo_results = search_geo(query, max_results=max_results)
        all_results.extend(geo_results)

        # Search ArrayExpress
        ae_results = search_arrayexpress(query, max_results=max_results // 2)
        all_results.extend(ae_results)

    # Deduplicate by accession
    seen = set()
    deduped = []
    for r in all_results:
        acc = r.get('accession', '')
        if acc and acc not in seen:
            seen.add(acc)
            deduped.append(r)
        elif not acc:
            deduped.append(r)

    # Sort: registry first (they're curated), then by sample count
    def sort_key(r):
        is_registry = 1 if r['source'] == 'registry' else 0
        try:
            n = int(str(r.get('n_samples', '0')).replace('~', '').replace(',', ''))
        except (ValueError, TypeError):
            n = 0
        return (is_registry, n)

    deduped.sort(key=sort_key, reverse=True)

    return deduped


def main():
    parser = argparse.ArgumentParser(
        description="TruthSeq Dataset Search: Find publicly available expression datasets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python3 dataset_search.py --query "autism brain RNA-seq"
    python3 dataset_search.py --query "CRISPR perturbation screen"
    python3 dataset_search.py --query "breast cancer differential expression"
    python3 dataset_search.py --query "schizophrenia" --update-registry
    python3 dataset_search.py --query "Parkinson" --registry-only
        """
    )
    parser.add_argument('--query', required=True, help='Search terms (disease, tissue, method, etc.)')
    parser.add_argument('--registry-only', action='store_true', help='Search local registry only, no API calls')
    parser.add_argument('--update-registry', action='store_true', help='Add new GEO/AE results to the local registry')
    parser.add_argument('--max-results', type=int, default=20, help='Max results per source (default: 20)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Show download instructions')
    parser.add_argument('--json', action='store_true', help='Output as JSON instead of text')

    args = parser.parse_args()

    results = find_datasets(
        args.query,
        registry_only=args.registry_only,
        max_results=args.max_results,
        verbose=args.verbose,
    )

    if args.json:
        print(json.dumps(results, indent=2, default=str))
    else:
        display_results(results, verbose=args.verbose)

    if args.update_registry and not args.registry_only:
        n_added = update_registry(results)
        if n_added > 0:
            print(f"  Updated registry: +{n_added} new datasets")

    return 0


if __name__ == '__main__':
    sys.exit(main())
