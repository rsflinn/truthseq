"""
TruthSeq Disease Expression Lookup
====================================

Finds disease-relevant differential expression data for Tier 2 validation.

Priority order:
  1. User-supplied file (--disease-expr)
  2. Local dataset registry match (dataset_registry.csv)
  3. Online search (GEO + ArrayExpress) → shows user what's available

If no local data exists, this module tells the user exactly what datasets
are available and how to download them. It does NOT try to auto-download
and auto-process arbitrary datasets — that's fragile and error-prone.

Returns a standardized DataFrame with columns:
  gene, log2fc, padj, cell_type, tissue, disease, source
"""

import os
import json
import logging
import tempfile

import pandas as pd

log = logging.getLogger(__name__)

# Standard column schema for all disease expression data
REQUIRED_COLUMNS = ['gene', 'log2fc', 'padj']
OPTIONAL_COLUMNS = ['cell_type', 'tissue', 'disease', 'source']


# ============================================================
# Column Standardization
# ============================================================

def standardize_columns(df):
    """Map common column name variants to TruthSeq standard names."""
    col_map = {}
    for col in df.columns:
        cl = col.lower().strip()
        if cl in ('gene', 'gene_name', 'gene_symbol', 'symbol', 'genename', 'external_gene_name'):
            col_map[col] = 'gene'
        elif cl in ('log2fc', 'log2foldchange', 'log2_fold_change', 'logfc', 'lfc', 'avg_log2fc'):
            col_map[col] = 'log2fc'
        elif cl in ('padj', 'p_val_adj', 'fdr', 'adj.p.val', 'q_value', 'qvalue', 'padjust', 'p.adjust'):
            col_map[col] = 'padj'
        elif cl in ('cell_type', 'celltype', 'cluster', 'cell_type_name', 'ident'):
            col_map[col] = 'cell_type'
        elif cl in ('tissue', 'tissue_type', 'organ'):
            col_map[col] = 'tissue'
        elif cl in ('disease', 'condition', 'phenotype'):
            col_map[col] = 'disease'
        elif cl in ('source', 'dataset_source', 'study', 'dataset'):
            col_map[col] = 'source'

    df = df.rename(columns=col_map)

    # Verify required columns
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns after standardization: {missing}. "
                         f"Available: {list(df.columns)}")

    # Add optional columns with defaults
    for col in OPTIONAL_COLUMNS:
        if col not in df.columns:
            df[col] = ''

    # Clean
    df['gene'] = df['gene'].astype(str).str.strip()
    df['log2fc'] = pd.to_numeric(df['log2fc'], errors='coerce')
    df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
    df = df.dropna(subset=['gene', 'log2fc', 'padj'])

    return df


# ============================================================
# User-Supplied File Loader
# ============================================================

def load_user_file(filepath):
    """Load a user-supplied disease expression file (TSV, CSV, or Parquet)."""
    log.info(f"Loading user-supplied disease expression file: {filepath}")

    if filepath.endswith('.parquet'):
        df = pd.read_parquet(filepath)
    elif filepath.endswith('.tsv') or filepath.endswith('.txt'):
        df = pd.read_csv(filepath, sep='\t')
    else:
        df = pd.read_csv(filepath)

    df = standardize_columns(df)
    log.info(f"  Loaded: {len(df)} entries, {df['gene'].nunique()} unique genes")
    return df


# ============================================================
# Disease Synonym Map (shared across functions)
# ============================================================

DISEASE_ALIASES = {
    'autism': ['autism', 'autistic', 'asd', 'autism spectrum'],
    'schizophrenia': ['schizophrenia', 'schizophren', 'scz'],
    'alzheimer': ['alzheimer', 'alzheimers', "alzheimer's", 'neurodegenerat'],
    'parkinson': ['parkinson', 'parkinsons', "parkinson's"],
    'bipolar': ['bipolar', 'manic depressive', 'mood disorder', 'mania'],
    'epilepsy': ['epilepsy', 'epileptic', 'seizure'],
    'depression': ['depression', 'depressive', 'mdd', 'major depressive'],
    'cancer': ['cancer', 'tumor', 'carcinoma', 'neoplasm', 'oncolog'],
    'diabetes': ['diabetes', 'diabetic', 'insulin resistance', 'glycemic'],
    'huntington': ['huntington', 'huntingtons', "huntington's"],
    'als': ['amyotrophic lateral sclerosis', 'motor neuron disease'],
    'fibromyalgia': ['fibromyalgia', 'fibromyalgic'],
    'crohn': ['crohn', 'crohns', "crohn's", 'inflammatory bowel'],
    'lupus': ['lupus', 'systemic lupus', 'sle'],
    'multiple sclerosis': ['multiple sclerosis', 'ms', 'demyelinat'],
}


def _get_search_terms(disease_keyword):
    """Get expanded search terms for a disease keyword."""
    keyword_lower = disease_keyword.lower().strip()
    for key, aliases in DISEASE_ALIASES.items():
        if keyword_lower in aliases or key in keyword_lower:
            return aliases
    return [keyword_lower]


# ============================================================
# Relevance Validation
# ============================================================

def check_dataset_relevance(description, disease_keyword):
    """
    Check whether a dataset's description is actually about the queried disease.

    Returns a dict with:
      - relevant (bool): True if description mentions the disease or its synonyms
      - confidence (str): 'high', 'low', or 'mismatch'
      - reason (str): human-readable explanation
      - detected_diseases (list): other diseases detected in the description
    """
    desc_lower = str(description).lower()
    search_terms = _get_search_terms(disease_keyword)

    # Check if the queried disease appears in the description
    query_hits = [t for t in search_terms if t in desc_lower]

    # Check what OTHER diseases appear in the description
    # Use word boundary checking to avoid substring matches (e.g. "als" in "reveals")
    import re
    other_diseases = []
    for disease_name, aliases in DISEASE_ALIASES.items():
        # Skip the queried disease itself
        if disease_keyword.lower() in aliases or disease_name in disease_keyword.lower():
            continue
        for alias in aliases:
            if len(alias) <= 3:
                # Short terms need word boundary matching
                if re.search(r'\b' + re.escape(alias) + r'\b', desc_lower):
                    other_diseases.append(disease_name)
                    break
            else:
                if alias in desc_lower:
                    other_diseases.append(disease_name)
                    break

    if query_hits:
        return {
            'relevant': True,
            'confidence': 'high',
            'reason': f"Description mentions: {', '.join(query_hits)}",
            'detected_diseases': other_diseases,
        }

    if other_diseases and not query_hits:
        return {
            'relevant': False,
            'confidence': 'mismatch',
            'reason': (
                f"Description does NOT mention {disease_keyword}. "
                f"Actual subject appears to be: {', '.join(other_diseases)}"
            ),
            'detected_diseases': other_diseases,
        }

    # No disease terms found at all — could be a general/healthy dataset
    return {
        'relevant': True,  # give benefit of the doubt for non-disease datasets
        'confidence': 'low',
        'reason': "No specific disease terms found in description (may be a control/healthy dataset)",
        'detected_diseases': [],
    }


# ============================================================
# Registry Search
# ============================================================

def search_registry(disease_keyword, registry_path=None):
    """
    Search the dataset registry for matching datasets.
    Returns a list of matching entries, sorted by relevance.
    Each entry includes a 'relevance_check' field with validation results.
    """
    if registry_path is None:
        registry_path = os.path.join(os.path.dirname(__file__), 'dataset_registry.csv')

    if not os.path.exists(registry_path):
        log.warning(f"Dataset registry not found: {registry_path}")
        return []

    df = pd.read_csv(registry_path)
    search_terms = _get_search_terms(disease_keyword)

    matches = []
    for idx, row in df.iterrows():
        searchable = ' '.join(str(v).lower() for v in row.values)
        score = sum(1 for term in search_terms if term in searchable)
        if score > 0:
            entry = row.to_dict()
            # Run relevance check against the description
            entry['relevance_check'] = check_dataset_relevance(
                entry.get('description', ''), disease_keyword
            )
            matches.append((score, entry))

    matches.sort(key=lambda x: x[0], reverse=True)
    return [m[1] for m in matches]


def display_registry_matches(matches, disease_keyword):
    """Print registry matches with download instructions and relevance warnings."""
    if not matches:
        return

    relevant = [m for m in matches if m.get('relevance_check', {}).get('relevant', True)]
    flagged = [m for m in matches if not m.get('relevance_check', {}).get('relevant', True)]

    log.info(f"")
    log.info(f"  Found {len(matches)} datasets in registry matching '{disease_keyword}':")
    if flagged:
        log.info(f"  ({len(flagged)} flagged as potentially mistagged — see warnings below)")
    log.info(f"")

    for i, m in enumerate(matches, 1):
        access_tag = "[OPEN]" if m.get('access_type') == 'open' else "[LOGIN REQUIRED]"
        relevance = m.get('relevance_check', {})

        # Show warning for mistagged datasets
        if not relevance.get('relevant', True):
            log.warning(f"  {i}. {m.get('accession', '')} {access_tag} — RELEVANCE WARNING")
            log.warning(f"     {relevance.get('reason', 'Description may not match queried disease')}")
            log.warning(f"     Tagged as: {m.get('disease', 'unknown')}")
            log.warning(f"     Actual description: {str(m.get('description', ''))[:120]}")
            log.warning(f"     This dataset will be SKIPPED for auto-analysis. Use --disease-expr to include it manually.")
            log.info(f"")
        else:
            log.info(f"  {i}. {m.get('accession', '')} {access_tag} — {m.get('paper_citation', '')}")
            log.info(f"     {m.get('data_type', '')}, {m.get('n_samples', '?')} samples, {m.get('species', '')}")
            log.info(f"     {str(m.get('description', ''))[:100]}")
            log.info(f"     URL: {m.get('download_url', '')}")
            log.info(f"")


# ============================================================
# Online Dataset Discovery
# ============================================================

def search_online(disease_keyword):
    """
    Search GEO and ArrayExpress for matching datasets.
    Returns results but does NOT download anything.
    """
    try:
        from dataset_search import search_geo, search_arrayexpress
    except ImportError:
        log.info("  dataset_search.py not found. Online search unavailable.")
        return []

    results = []

    geo_results = search_geo(disease_keyword, max_results=10)
    results.extend(geo_results)

    ae_results = search_arrayexpress(disease_keyword, max_results=5)
    results.extend(ae_results)

    return results


# ============================================================
# Main lookup function
# ============================================================

def find_disease_expression(disease=None, disease_expr_file=None,
                             catalog_path=None, registry_path=None,
                             cache_dir=None):
    """
    Main entry point: find disease expression data.

    Priority:
      1. User-supplied file (--disease-expr)
      2. Local registry → guide user to download
      3. Online search → show user what exists

    Returns:
      (DataFrame, source_description) or (None, reason_string)
    """
    # Priority 1: User-supplied file
    if disease_expr_file and os.path.exists(disease_expr_file):
        try:
            df = load_user_file(disease_expr_file)
            return df, f"User-supplied file: {disease_expr_file}"
        except Exception as e:
            log.error(f"Could not load user file {disease_expr_file}: {e}")
            return None, f"Error loading {disease_expr_file}: {e}"

    if not disease:
        log.info("No --disease or --disease-expr specified. Tier 2 validation skipped.")
        return None, "No disease context specified"

    log.info(f"Searching for '{disease}' expression data...")

    # Priority 2: Check local registry for known datasets
    log.info(f"  Step 1: Checking dataset registry...")
    matches = search_registry(disease, registry_path=registry_path)

    if matches:
        display_registry_matches(matches, disease)

        # Check if any match has a local file already downloaded
        # (look for parquet/csv/tsv files in common locations)
        # SKIP datasets that failed the relevance check
        script_dir = os.path.dirname(os.path.abspath(__file__))
        for m in matches:
            relevance = m.get('relevance_check', {})
            if not relevance.get('relevant', True):
                log.info(f"  Skipping {m.get('dataset_id', '?')} — {relevance.get('reason', 'failed relevance check')}")
                continue

            dataset_id = m.get('dataset_id', '')
            for ext in ['.parquet', '.tsv', '.csv']:
                for search_dir in [script_dir, os.path.join(script_dir, 'data'),
                                   os.path.join(script_dir, 'Inputs')]:
                    candidate = os.path.join(search_dir, f"{dataset_id}{ext}")
                    if os.path.exists(candidate):
                        log.info(f"  Found local file: {candidate}")
                        log.info(f"  Relevance: {relevance.get('confidence', '?')} — {relevance.get('reason', '')}")
                        try:
                            df = load_user_file(candidate)
                            return df, f"Registry match (local): {m.get('paper_citation', dataset_id)}"
                        except Exception as e:
                            log.warning(f"  Could not load {candidate}: {e}")

        # No local file found — tell user how to get the data
        log.info(f"  No local data file found for '{disease}'.")
        log.info(f"  To use Tier 2 validation, download one of the datasets above")
        log.info(f"  and re-run with: --disease-expr <your_file.tsv>")
        log.info(f"  See format_spec.md for the expected file format.")
        log.info(f"")

        # Return None but with helpful context
        best_match = matches[0]
        return None, (
            f"Registry has {len(matches)} matching datasets for '{disease}' but none are "
            f"downloaded locally. Best match: {best_match.get('paper_citation', best_match.get('accession', ''))} "
            f"({best_match.get('access_type', 'unknown')} access). "
            f"Download and re-run with --disease-expr <file>."
        )

    # Priority 3: Online search
    log.info(f"  Step 2: Searching online databases (GEO, ArrayExpress)...")
    online_results = search_online(disease)

    if online_results:
        log.info(f"  Found {len(online_results)} datasets online:")
        for i, r in enumerate(online_results[:5], 1):
            log.info(f"    {i}. {r['accession']} — {str(r['description'])[:80]}")
            log.info(f"       URL: {r['download_url']}")
        log.info(f"")
        log.info(f"  To use these, download the DE results and run with --disease-expr <file>.")
        log.info(f"  Run 'python3 dataset_search.py --query \"{disease}\" --verbose' for full details.")
        log.info(f"  Run with --update-registry to save these results for future searches.")
        return None, f"Found {len(online_results)} datasets online for '{disease}' but none downloaded locally."

    # Nothing found anywhere
    log.info(f"  No matching datasets found for '{disease}'.")
    log.info(f"  Try searching manually: python3 dataset_search.py --query \"{disease}\"")
    return None, f"No datasets found for '{disease}'"
