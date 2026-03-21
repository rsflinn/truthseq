#!/usr/bin/env python3
"""
TruthSeq Setup: Download and prepare the Replogle Perturb-seq reference data.
=============================================================================

This is the first thing a new user runs after cloning the repo. It:
  1. Downloads the Replogle genome-wide Perturb-seq pseudo-bulk h5ad from Figshare
  2. Processes it into parquet files that TruthSeq uses for validation
  3. Optionally builds a gene ID mapping for Open Targets (Tier 3)

The Replogle atlas contains ~11,000 single-gene knockdowns in human K562 cells,
measuring expression changes across ~8,000 genes per knockdown. This is the
core reference data: for any gene-gene regulatory claim, TruthSeq checks
whether knocking down the upstream gene actually changed expression of the
downstream gene in real human cells.

After setup, run:
    python3 truthseq_validate.py --claims your_claims.csv \\
        --replogle replogle_knockdown_effects.parquet \\
        --replogle-stats replogle_knockdown_stats.parquet

Usage:
    python3 setup.py                    # Full setup
    python3 setup.py --skip-download    # Re-process existing h5ad
    python3 setup.py --status           # Check what's already set up

Requirements:
    pip3 install scanpy anndata pandas pyarrow numpy scipy requests
"""

import os
import sys
import argparse
import logging
import hashlib
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

# Figshare download URL for the Replogle pseudo-bulk h5ad
# Source: Replogle et al. 2022, "Mapping information-rich genotype-phenotype
# landscapes with genome-scale Perturb-seq" (Cell)
FIGSHARE_URL = "https://figshare.com/ndownloader/files/35773219"
H5AD_FILENAME = "K562_gwps_normalized_bulk_01.h5ad"
H5AD_SIZE_MB = 2700  # ~2.7 GB

# Output filenames
EFFECTS_PARQUET = "replogle_knockdown_effects.parquet"
STATS_PARQUET = "replogle_knockdown_stats.parquet"
GENE_MAP_FILE = "gene_id_mapping.tsv"


def check_status(work_dir):
    """Report what data files exist and their stats."""
    print("\n=== TruthSeq Data Status ===\n")

    files = {
        H5AD_FILENAME: "Replogle h5ad (raw source, ~2.7 GB)",
        EFFECTS_PARQUET: "Knockdown effects parquet (Tier 1 data)",
        STATS_PARQUET: "Knockdown distribution stats (for percentile calculations)",
        GENE_MAP_FILE: "Gene ID mapping (for Open Targets, Tier 3)",
    }

    all_ready = True
    for fname, desc in files.items():
        fpath = os.path.join(work_dir, fname)
        if os.path.exists(fpath):
            size_mb = os.path.getsize(fpath) / 1024 / 1024
            print(f"  [OK] {fname} ({size_mb:.1f} MB) — {desc}")
        else:
            print(f"  [--] {fname} — {desc}")
            if fname in (EFFECTS_PARQUET, STATS_PARQUET):
                all_ready = False

    print()
    if all_ready:
        print("Ready to validate! Run:")
        print(f"  python3 truthseq_validate.py --claims your_claims.csv \\")
        print(f"      --replogle {EFFECTS_PARQUET} \\")
        print(f"      --replogle-stats {STATS_PARQUET}")
    else:
        print("Setup needed. Run: python3 setup.py")

    print()
    return all_ready


def download_h5ad(work_dir):
    """Download the Replogle h5ad from Figshare."""
    import requests

    output_path = os.path.join(work_dir, H5AD_FILENAME)

    if os.path.exists(output_path):
        size_mb = os.path.getsize(output_path) / 1024 / 1024
        if size_mb > 100:  # Sanity check: file should be >100 MB
            log.info(f"h5ad file already exists: {output_path} ({size_mb:.0f} MB)")
            return output_path
        else:
            log.warning(f"h5ad file exists but is suspiciously small ({size_mb:.1f} MB). Re-downloading.")

    print()
    print("=" * 70)
    print("Step 1: Download Replogle Perturb-seq data from Figshare")
    print("=" * 70)
    print(f"  Source: Replogle et al. 2022 (Cell)")
    print(f"  File size: ~{H5AD_SIZE_MB / 1000:.1f} GB")
    print(f"  URL: {FIGSHARE_URL}")
    print()

    log.info("Starting download (this may take 10-30 minutes depending on your connection)...")

    try:
        resp = requests.get(FIGSHARE_URL, stream=True, timeout=30)
        resp.raise_for_status()

        total_size = int(resp.headers.get('content-length', 0))
        downloaded = 0
        start_time = time.time()

        with open(output_path, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):  # 1 MB chunks
                f.write(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    pct = downloaded / total_size * 100
                    elapsed = time.time() - start_time
                    speed_mbps = (downloaded / 1024 / 1024) / max(elapsed, 0.1)
                    remaining = (total_size - downloaded) / (downloaded / max(elapsed, 0.1)) if downloaded > 0 else 0
                    print(f"\r  {pct:.0f}% ({downloaded/1024/1024:.0f}/{total_size/1024/1024:.0f} MB) "
                          f"— {speed_mbps:.1f} MB/s, ~{remaining/60:.0f} min remaining", end='', flush=True)

        print()
        size_mb = os.path.getsize(output_path) / 1024 / 1024
        elapsed = time.time() - start_time
        log.info(f"Download complete: {size_mb:.0f} MB in {elapsed/60:.1f} minutes")
        return output_path

    except Exception as e:
        log.error(f"Download failed: {e}")
        log.error("")
        log.error("You can download manually instead:")
        log.error(f"  1. Go to: {FIGSHARE_URL}")
        log.error(f"  2. Save as: {output_path}")
        log.error(f"  3. Re-run: python3 setup.py --skip-download")
        if os.path.exists(output_path):
            os.remove(output_path)  # Clean up partial download
        return None


def process_h5ad(h5ad_path, work_dir):
    """
    Process the h5ad into parquet files.
    Extracts knockdown effects and per-knockdown distribution statistics.
    """
    import re
    import scanpy as sc
    from scipy import stats as sp_stats

    print()
    print("=" * 70)
    print("Step 2: Process h5ad into TruthSeq format")
    print("=" * 70)
    print()

    effects_path = os.path.join(work_dir, EFFECTS_PARQUET)
    stats_path = os.path.join(work_dir, STATS_PARQUET)

    import numpy as np
    import pandas as pd

    log.info(f"Loading {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)
    log.info(f"  Shape: {adata.shape[0]} observations x {adata.shape[1]} genes")
    log.info(f"  Obs columns: {list(adata.obs.columns)}")

    # ---- Parse gene names from obs index ----
    # Format: NUMBER_GENENAME_PERTURBATION_ENSEMBLID
    # e.g., 3291_GATA1_P1P2_ENSG00000102145
    raw_index = adata.obs_names.values.astype(str)

    def clean_gene_name(name):
        name = re.sub(r'^\d+_', '', name)           # Strip leading number
        name = re.sub(r'_ENSG\d+', '', name)         # Strip trailing Ensembl ID
        name = re.sub(r'_P\d+(P\d+)?', '', name)     # Strip perturbation suffix
        return name

    gene_labels = np.array([clean_gene_name(name) for name in raw_index])
    unique_genes = set(gene_labels)
    log.info(f"  Unique knockdown genes (after name cleaning): {len(unique_genes)}")

    # ---- Map affected gene columns (Ensembl IDs -> symbols) ----
    if 'gene_name' in adata.var.columns:
        col_names = adata.var['gene_name'].values.astype(str)
        log.info(f"  Using gene_name column from adata.var for affected genes")
    else:
        col_names = adata.var_names.values.astype(str)
        log.info(f"  Using var_names index for affected genes (may be Ensembl IDs)")

    # ---- Build expression matrix ----
    log.info("Building expression matrix...")

    expr_data = adata.X
    if hasattr(expr_data, 'toarray'):
        expr_data = expr_data.toarray()

    expr_matrix = pd.DataFrame(expr_data, index=gene_labels, columns=col_names)

    # Handle duplicate column names (multiple Ensembl IDs -> same gene symbol)
    if expr_matrix.columns.duplicated().any():
        n_dups = expr_matrix.columns.duplicated().sum()
        log.info(f"  Averaging {n_dups} duplicate gene symbols...")
        expr_matrix = expr_matrix.T.groupby(level=0).mean().T

    log.info(f"  Expression matrix: {expr_matrix.shape}")

    # ---- Detect if data needs Z-scoring ----
    sample_vals = expr_matrix.iloc[:100].values.flatten()
    sample_vals = sample_vals[~np.isnan(sample_vals)]
    grand_median = float(np.median(np.abs(sample_vals)))
    log.info(f"  Grand median |value|: {grand_median:.4f}")

    needs_zscore = grand_median < 0.5
    if needs_zscore:
        log.info("  Data appears to be fold-changes (not Z-scores). Will Z-score per knockdown.")
    else:
        log.info("  Data appears to already be Z-scored.")

    # ---- Process each knockdown ----
    log.info("Processing knockdowns...")
    unique_kd_genes = sorted(set(gene_labels))
    all_pairs = []
    all_stats = []
    z_threshold = 1.0

    for i, kd_gene in enumerate(unique_kd_genes):
        if (i + 1) % 500 == 0:
            log.info(f"  {i+1}/{len(unique_kd_genes)} knockdowns processed...")

        mask = gene_labels == kd_gene
        kd_expr = expr_matrix.loc[mask]

        if len(kd_expr) == 0:
            continue

        # Average across replicates
        mean_effects = kd_expr.mean(axis=0)

        if needs_zscore:
            kd_mean = float(mean_effects.mean())
            kd_std = float(mean_effects.std())
            if kd_std > 0:
                z_effects = (mean_effects - kd_mean) / kd_std
            else:
                z_effects = mean_effects * 0
        else:
            z_effects = mean_effects

        abs_z = z_effects.abs()

        # Store distribution stats for this knockdown
        stats_row = {
            'knocked_down_gene': kd_gene,
            'n_genes_tested': len(z_effects),
            'n_replicates': int(mask.sum()),
            'median_abs_z': float(abs_z.median()),
            'mean_abs_z': float(abs_z.mean()),
            'std_abs_z': float(abs_z.std()),
            'max_abs_z': float(abs_z.max()),
            'n_sig': int((abs_z > z_threshold).sum()),
        }

        # Add quantile breakpoints
        for q in [5, 10, 25, 50, 75, 80, 85, 90, 95, 97, 99]:
            stats_row[f'q{q:02d}'] = float(np.percentile(abs_z.values, q))

        all_stats.append(stats_row)

        # Store significant pairs (|Z| > threshold)
        sig_mask = abs_z > z_threshold
        sig_genes = z_effects[sig_mask]

        for affected_gene, z_score in sig_genes.items():
            if affected_gene == kd_gene:
                continue  # Skip self-effects
            all_pairs.append({
                'knocked_down_gene': kd_gene,
                'affected_gene': affected_gene,
                'z_score': round(float(z_score), 4),
                'cell_line': 'K562',
            })

    log.info(f"  Total knockdowns: {len(all_stats)}")
    log.info(f"  Total significant pairs (|Z|>{z_threshold}): {len(all_pairs):,}")

    # ---- Save parquets ----
    effects_df = pd.DataFrame(all_pairs)
    effects_df.to_parquet(effects_path, index=False)
    log.info(f"  Saved: {effects_path} ({os.path.getsize(effects_path)/1024/1024:.1f} MB)")

    stats_df = pd.DataFrame(all_stats)
    stats_df.to_parquet(stats_path, index=False)
    log.info(f"  Saved: {stats_path} ({os.path.getsize(stats_path)/1024/1024:.1f} MB)")

    return effects_path, stats_path


def build_gene_map(work_dir):
    """Build gene symbol -> Ensembl ID mapping for Open Targets queries."""
    gene_map_path = os.path.join(work_dir, GENE_MAP_FILE)

    if os.path.exists(gene_map_path):
        import pandas as pd
        existing = pd.read_csv(gene_map_path, sep='\t')
        log.info(f"Gene map already exists: {len(existing)} genes")
        return gene_map_path

    print()
    print("=" * 70)
    print("Step 3: Build gene ID mapping (for Open Targets, Tier 3)")
    print("=" * 70)

    # Try to extract from the h5ad if it's still loaded, or from the parquet
    h5ad_path = os.path.join(work_dir, H5AD_FILENAME)
    effects_path = os.path.join(work_dir, EFFECTS_PARQUET)

    all_genes = set()

    if os.path.exists(effects_path):
        import pandas as pd
        df = pd.read_parquet(effects_path)
        all_genes.update(df['knocked_down_gene'].unique())
        all_genes.update(df['affected_gene'].unique())

    if not all_genes:
        log.warning("No gene list available. Skipping gene map build.")
        return None

    log.info(f"Querying Ensembl BioMart for {len(all_genes)} gene symbols...")

    try:
        import requests

        # BioMart query for human gene symbols -> Ensembl IDs
        biomart_url = "http://www.ensembl.org/biomart/martservice"
        xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1"
       uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="hsapiens_gene_ensembl" interface="default">
        <Attribute name="hgnc_symbol"/>
        <Attribute name="ensembl_gene_id"/>
    </Dataset>
</Query>"""

        resp = requests.get(biomart_url, params={'query': xml_query}, timeout=120)

        if resp.status_code == 200 and len(resp.text) > 100:
            import io
            import pandas as pd
            biomart_df = pd.read_csv(io.StringIO(resp.text), sep='\t')
            biomart_df.columns = ['symbol', 'ensembl_id']
            biomart_df = biomart_df[biomart_df['symbol'].isin(all_genes)]
            biomart_df = biomart_df.drop_duplicates(subset='symbol', keep='first')
            biomart_df.to_csv(gene_map_path, sep='\t', index=False)
            log.info(f"  Saved: {gene_map_path} ({len(biomart_df)} genes mapped)")
            return gene_map_path
        else:
            log.warning(f"  BioMart query failed (status {resp.status_code}). Gene map not built.")
            log.info("  Tier 3 (Open Targets) will still work but may be slower.")
            return None

    except Exception as e:
        log.warning(f"  Could not build gene map: {e}")
        log.info("  This is optional. Tier 1 validation works without it.")
        return None


def main():
    import pandas as pd
    import numpy as np

    parser = argparse.ArgumentParser(
        description="TruthSeq Setup: Download and prepare Replogle Perturb-seq reference data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
After setup, validate your claims:
    python3 truthseq_validate.py --claims your_claims.csv \\
        --replogle replogle_knockdown_effects.parquet \\
        --replogle-stats replogle_knockdown_stats.parquet

Optional: add disease tissue expression for Tier 2:
    python3 truthseq_validate.py --claims your_claims.csv \\
        --replogle replogle_knockdown_effects.parquet \\
        --replogle-stats replogle_knockdown_stats.parquet \\
        --disease-expr your_disease_de_results.tsv

See format_spec.md for disease expression file format.
        """
    )
    parser.add_argument('--status', action='store_true',
                        help='Check what data files are already set up')
    parser.add_argument('--skip-download', action='store_true',
                        help='Skip h5ad download (re-process existing file)')
    parser.add_argument('--skip-gene-map', action='store_true',
                        help='Skip gene ID mapping (Tier 3 optional)')
    parser.add_argument('--dir', default='.',
                        help='Working directory (default: current)')

    args = parser.parse_args()
    work_dir = os.path.abspath(args.dir)

    if args.status:
        check_status(work_dir)
        return 0

    print()
    print("=" * 70)
    print("  TruthSeq Setup")
    print("  Preparing Replogle Perturb-seq reference data")
    print("=" * 70)
    print()
    print("  This will download ~2.7 GB of Perturb-seq data from Figshare")
    print("  and process it into parquet files for fast validation.")
    print()
    print(f"  Working directory: {work_dir}")
    print()

    # Step 1: Download
    h5ad_path = os.path.join(work_dir, H5AD_FILENAME)
    if args.skip_download:
        if not os.path.exists(h5ad_path):
            log.error(f"h5ad not found at {h5ad_path}")
            log.error("Run without --skip-download to download it first.")
            return 1
        log.info(f"Skipping download. Using existing: {h5ad_path}")
    else:
        h5ad_path = download_h5ad(work_dir)
        if h5ad_path is None:
            return 1

    # Step 2: Process
    effects_path = os.path.join(work_dir, EFFECTS_PARQUET)
    stats_path = os.path.join(work_dir, STATS_PARQUET)

    if os.path.exists(effects_path) and os.path.exists(stats_path):
        log.info(f"Parquet files already exist. Skipping processing.")
        log.info(f"  (Delete {EFFECTS_PARQUET} and {STATS_PARQUET} to force reprocessing)")
    else:
        process_h5ad(h5ad_path, work_dir)

    # Step 3: Gene map (optional)
    if not args.skip_gene_map:
        build_gene_map(work_dir)

    # Summary
    print()
    print("=" * 70)
    print("  Setup Complete!")
    print("=" * 70)
    check_status(work_dir)

    return 0


if __name__ == '__main__':
    sys.exit(main())
