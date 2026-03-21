"""
TruthSeq v2: Download and process Replogle Perturb-seq data
=============================================================
Run this on your local Mac (not in the Cowork VM).

v2 CHANGE: In addition to significant gene-gene pairs (|Z|>1), this script
now saves per-knockdown distribution statistics. This lets TruthSeq compute
accurate percentile ranks for ANY downstream gene — including those with
modest effects that fall below the |Z|>1 threshold.

This solves the v1 problem where known-true regulatory relationships
(MYT1L→MEF2C, GATA1→HBB) scored WEAK because their effects were filtered out.

TWO APPROACHES (same as v1):

  APPROACH A (recommended): Figshare pseudo-bulk h5ad.
    - Download ~2-3 GB h5ad file manually from Figshare
    - Produces BOTH the hits file AND full distribution stats
    - Most accurate percentile calculations

  APPROACH B (fallback): Harmonizome API.
    - No large download needed
    - Only gets pre-filtered significant associations
    - Distribution stats are estimated from the data (less accurate)
    - Still much better than v1 because it stores the estimates

Requirements:
    pip3 install scanpy anndata pandas pyarrow requests tqdm numpy scipy

Outputs:
    replogle_knockdown_effects.parquet   — gene-gene pairs with |Z|>1
    replogle_knockdown_stats.parquet     — NEW: per-knockdown distribution stats
"""

import os
import sys
import argparse
import logging
import time

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)


# ============================================================
# APPROACH A: Figshare pseudo-bulk h5ad (recommended)
# ============================================================

def download_figshare_files(output_dir):
    """Prompt user to download h5ad from Figshare."""
    os.makedirs(output_dir, exist_ok=True)

    expected_file = os.path.join(output_dir, "K562_gwps_normalized_bulk_01.h5ad")

    if not os.path.exists(expected_file):
        log.info("=" * 70)
        log.info("MANUAL DOWNLOAD REQUIRED")
        log.info("=" * 70)
        log.info("")
        log.info("Visit this URL in your browser:")
        log.info("https://plus.figshare.com/articles/dataset/20029387")
        log.info("")
        log.info("Download: K562_gwps_normalized_bulk_01.h5ad")
        log.info(f"Save it to: {output_dir}")
        log.info("")
        log.info("File size is approximately 2-3 GB.")
        log.info("After downloading, re-run this script.")
        log.info("=" * 70)
        return False

    return True


def process_h5ad_v2(input_dir, output_dir):
    """
    Convert Replogle pseudo-bulk h5ad into:
      1. Gene-gene effect lookup (|Z|>1 pairs) — same as v1
      2. Per-knockdown distribution stats — NEW in v2

    The distribution stats file stores quantile breakpoints for each
    knockdown, enabling accurate percentile calculations for any
    downstream gene without storing millions of pairs.
    """
    import anndata as ad

    log.info("Loading K562 genome-wide pseudo-bulk data...")
    h5ad_path = os.path.join(input_dir, "K562_gwps_normalized_bulk_01.h5ad")
    adata = ad.read_h5ad(h5ad_path)

    log.info(f"  Shape: {adata.shape} (perturbations x genes)")
    log.info(f"  Obs columns: {list(adata.obs.columns)}")

    # Find the gene column
    gene_col = None
    for candidate in ['gene', 'target_gene', 'perturbation', 'gene_name',
                       'knocked_down', 'sgRNA_group', 'gene_id']:
        if candidate in adata.obs.columns:
            gene_col = candidate
            break

    if gene_col is None:
        log.error(f"Could not find gene column. Available: {list(adata.obs.columns)}")
        return False

    log.info(f"  Using '{gene_col}' as the knocked-down gene identifier")

    # Extract expression matrix
    log.info("Extracting expression matrix...")
    if hasattr(adata.X, 'toarray'):
        expr_matrix = pd.DataFrame(
            adata.X.toarray(), index=adata.obs.index, columns=adata.var_names
        )
    else:
        expr_matrix = pd.DataFrame(
            adata.X, index=adata.obs.index, columns=adata.var_names
        )

    expr_matrix[gene_col] = adata.obs[gene_col].values

    # Process each knockdown
    log.info("Building gene-gene effects AND per-knockdown stats...")
    unique_genes = expr_matrix[gene_col].unique()
    n_genes = len(unique_genes)
    log.info(f"  Total perturbations: {n_genes}")

    effect_records = []  # gene-gene pairs with |Z|>1
    stats_records = []   # per-knockdown distribution stats

    # Quantile breakpoints to store (enough for accurate interpolation)
    quantile_points = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50,
                       0.60, 0.70, 0.80, 0.90, 0.95, 0.99]

    for i, kd_gene in enumerate(unique_genes):
        if i % 500 == 0:
            log.info(f"  Processing {i}/{n_genes} perturbations...")

        kd_rows = expr_matrix[expr_matrix[gene_col] == kd_gene]
        mean_effects = kd_rows.drop(columns=[gene_col]).mean(axis=0)

        # --- Stats for this knockdown (ALL genes) ---
        abs_z = mean_effects.abs()
        quantiles = np.quantile(abs_z.values, quantile_points)

        stats_row = {
            'knocked_down_gene': kd_gene,
            'n_genes_tested': len(mean_effects),
            'mean_abs_z': round(float(abs_z.mean()), 4),
            'std_abs_z': round(float(abs_z.std()), 4),
            'median_abs_z': round(float(abs_z.median()), 4),
        }
        for q, val in zip(quantile_points, quantiles):
            stats_row[f'q{int(q*100):02d}'] = round(float(val), 4)

        stats_records.append(stats_row)

        # --- Significant effects (|Z|>1) ---
        significant = mean_effects[mean_effects.abs() > 1.0]
        for affected_gene, z_score in significant.items():
            effect_records.append({
                'knocked_down_gene': kd_gene,
                'affected_gene': affected_gene,
                'z_score': round(float(z_score), 4),
                'cell_line': 'K562'
            })

    # Save effects file
    effects_df = pd.DataFrame(effect_records)
    effects_path = os.path.join(output_dir, 'replogle_knockdown_effects.parquet')
    effects_df.to_parquet(effects_path, index=False)
    log.info(f"  Effects file: {len(effects_df):,} pairs (|Z|>1)")
    log.info(f"    Unique knockdowns: {effects_df['knocked_down_gene'].nunique():,}")
    log.info(f"    Saved to {effects_path}")

    # Save stats file
    stats_df = pd.DataFrame(stats_records)
    stats_path = os.path.join(output_dir, 'replogle_knockdown_stats.parquet')
    stats_df.to_parquet(stats_path, index=False)
    log.info(f"  Stats file: {len(stats_df)} knockdowns with full distribution data")
    log.info(f"    Saved to {stats_path}")

    return True


# ============================================================
# APPROACH B: Harmonizome API (fallback)
# ============================================================

def process_via_harmonizome_v2(output_dir):
    """
    Query Harmonizome API for Replogle perturbation signatures.

    v2 change: Also estimates per-knockdown distribution stats from
    the available data. Since Harmonizome only returns significant
    associations, the stats are estimated (less accurate than Approach A).
    """
    import requests
    import re

    BASE = "https://maayanlab.cloud/Harmonizome/api/1.0"
    DATASET = "Replogle et al., Cell, 2022 K562 Genome-wide Perturb-seq Gene Perturbation Signatures"

    log.info("Querying Harmonizome API for Replogle K562 signatures...")

    resp = requests.get(f"{BASE}/dataset/{DATASET}", timeout=30)
    if resp.status_code != 200:
        log.error(f"Failed to fetch dataset: {resp.status_code}")
        return False

    data = resp.json()
    gene_sets = data.get('geneSets', [])
    log.info(f"  Found {len(gene_sets)} perturbation signatures")
    log.info(f"  Estimated time: {len(gene_sets) * 0.15 / 60:.0f} minutes")

    def extract_gene_from_name(name):
        cleaned = re.sub(r'_P\d*(P\d*)?$', '', name)
        parts = cleaned.split('_', 1)
        return parts[1] if len(parts) >= 2 else name

    effect_records = []
    stats_records = []

    # Approximate total genes tested per knockdown in the full Replogle dataset
    APPROX_TOTAL_GENES = 8000

    for i, gs in enumerate(gene_sets):
        if i % 100 == 0:
            log.info(f"  Fetching {i}/{len(gene_sets)}...")

        name = gs.get('name', '')
        kd_gene = extract_gene_from_name(name)
        gs_href = gs.get('href', '')

        if not gs_href:
            continue

        try:
            gs_resp = requests.get(
                f"https://maayanlab.cloud/Harmonizome{gs_href}", timeout=15
            )
            if gs_resp.status_code != 200:
                continue

            gs_data = gs_resp.json()
            kd_z_scores = []

            for assoc in gs_data.get('associations', []):
                gene_info = assoc.get('gene', {})
                affected_gene = gene_info.get('symbol', '')
                value = assoc.get('standardizedValue', 0)

                if affected_gene and abs(value) > 0:
                    effect_records.append({
                        'knocked_down_gene': kd_gene,
                        'affected_gene': affected_gene,
                        'z_score': round(float(value), 4),
                        'cell_line': 'K562'
                    })
                    kd_z_scores.append(abs(float(value)))

            # Estimate distribution stats
            # We only see the significant tail. Estimate the full distribution
            # by assuming the remaining genes have |Z| ~ half-normal(0, 0.5)
            n_sig = len(kd_z_scores)
            n_nonsig = max(0, APPROX_TOTAL_GENES - n_sig)

            if kd_z_scores:
                # Generate estimated non-significant portion
                np.random.seed(hash(kd_gene) % 2**31)
                nonsig_z = np.abs(np.random.normal(0, 0.5, n_nonsig))
                full_dist = np.concatenate([np.array(kd_z_scores), nonsig_z])

                quantile_points = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50,
                                   0.60, 0.70, 0.80, 0.90, 0.95, 0.99]
                quantiles = np.quantile(full_dist, quantile_points)

                stats_row = {
                    'knocked_down_gene': kd_gene,
                    'n_genes_tested': APPROX_TOTAL_GENES,
                    'mean_abs_z': round(float(np.mean(full_dist)), 4),
                    'std_abs_z': round(float(np.std(full_dist)), 4),
                    'median_abs_z': round(float(np.median(full_dist)), 4),
                    '_estimated': True,  # Flag that these are estimated, not exact
                }
                for q, val in zip(quantile_points, quantiles):
                    stats_row[f'q{int(q*100):02d}'] = round(float(val), 4)

                stats_records.append(stats_row)

            time.sleep(0.1)

        except Exception as e:
            log.warning(f"  Error at {i}: {e}")
            continue

    # Save effects
    effects_df = pd.DataFrame(effect_records)
    effects_path = os.path.join(output_dir, 'replogle_knockdown_effects.parquet')
    effects_df.to_parquet(effects_path, index=False)
    log.info(f"  Effects: {len(effects_df):,} pairs from {effects_df['knocked_down_gene'].nunique():,} knockdowns")

    # Save stats
    stats_df = pd.DataFrame(stats_records)
    stats_path = os.path.join(output_dir, 'replogle_knockdown_stats.parquet')
    stats_df.to_parquet(stats_path, index=False)
    log.info(f"  Stats: {len(stats_df)} knockdowns (estimated distributions)")
    log.info(f"  NOTE: Distribution stats are ESTIMATED from Harmonizome filtered data.")
    log.info(f"  For accurate stats, use Approach A (Figshare h5ad).")

    return True


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="TruthSeq v2: Download and process Replogle Perturb-seq data"
    )
    parser.add_argument('--approach', choices=['figshare', 'harmonizome'],
                        default='figshare')
    parser.add_argument('--input-dir', default='./raw_data',
                        help='Directory for h5ad files (figshare approach)')
    parser.add_argument('--output-dir', default='.',
                        help='Output directory for parquet files')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    if args.approach == 'figshare':
        log.info("=== Approach A: Figshare h5ad (full distribution) ===")
        if download_figshare_files(args.input_dir):
            success = process_h5ad_v2(args.input_dir, args.output_dir)
        else:
            log.info("")
            log.info("Use --approach harmonizome to skip the large download.")
            success = False
    else:
        log.info("=== Approach B: Harmonizome API (estimated distributions) ===")
        success = process_via_harmonizome_v2(args.output_dir)

    if success:
        log.info("")
        log.info("SUCCESS! Output files:")
        log.info(f"  {args.output_dir}/replogle_knockdown_effects.parquet")
        log.info(f"  {args.output_dir}/replogle_knockdown_stats.parquet")

    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())
