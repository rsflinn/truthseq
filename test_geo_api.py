#!/usr/bin/env python3
"""
Quick test: verify GEO E-utilities API works from this machine.
Run: python3 test_geo_api.py
"""

import sys

def test_geo():
    import requests

    print("Testing GEO E-utilities API...")
    try:
        resp = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={
                'db': 'gds',
                'term': '(autism spectrum disorder) AND ("expression profiling by high throughput sequencing"[DataSet Type])',
                'retmax': 5,
                'retmode': 'json',
                'sort': 'relevance',
            },
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        count = int(data.get('esearchresult', {}).get('count', 0))
        ids = data.get('esearchresult', {}).get('idlist', [])
        print(f"  GEO search: {count} total results, {len(ids)} returned")
        if count > 0:
            print("  GEO API: OK")
        else:
            print("  GEO API: responded but zero results (unexpected)")
    except Exception as e:
        print(f"  GEO API: FAILED - {e}")
        return False

    print("\nTesting ArrayExpress API...")
    try:
        resp = requests.get(
            "https://www.ebi.ac.uk/biostudies/api/v1/search",
            params={'query': 'autism', 'type': 'study', 'pageSize': 3},
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        total = data.get('totalHits', 0)
        print(f"  ArrayExpress search: {total} total results")
        if total > 0:
            print("  ArrayExpress API: OK")
        else:
            print("  ArrayExpress API: responded but zero results (unexpected)")
    except Exception as e:
        print(f"  ArrayExpress API: FAILED - {e}")
        return False

    print("\nTesting dataset_search.py full search...")
    try:
        from dataset_search import find_datasets
        results = find_datasets("autism brain RNA-seq", max_results=5)
        print(f"  Combined search returned {len(results)} results")
        for r in results[:3]:
            print(f"    {r['accession']} ({r['source']}): {r['description'][:80]}")
        print("  Full search: OK")
    except Exception as e:
        print(f"  Full search: FAILED - {e}")
        return False

    print("\nAll tests passed.")
    return True


if __name__ == '__main__':
    success = test_geo()
    sys.exit(0 if success else 1)
