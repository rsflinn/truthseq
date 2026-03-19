#!/usr/bin/env python3
"""
TruthSeq: Validate computational gene regulatory claims against experimental data.
==================================================================================

Takes a CSV of gene-gene regulatory predictions and checks each one against:
  1. Direct perturbation data (Replogle Perturb-seq atlas)
  2. Disease tissue expression (PsychENCODE ASD differential expression)
  3. Genetic associations (Open Targets Platform API)

Outputs a confidence grade for each claim: VALIDATED, PARTIALLY_SUPPORTED,
CONTRADICTED, WEAK, UNTESTABLE, or CELL_TYPE_CAVEAT.

Usage:
    python truthseq_validate.py \
        --claims claims.csv \
        --replogle replogle_knockdown_effects.parquet \
        --psychencode psychencode_asd_de.parquet \
        --gene-map gene_id_mapping.tsv \
        --output truthseq_report
"""

import os
import sys
import argparse
import logging
import json
from datetime import datetime

import pandas as pd
import numpy as np
from scipy import stats

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

RANDOM_SEED = 42
NULL_SAMPLE_SIZE = 500
PERCENTILE_VALIDATED = 90
PERCENTILE_PARTIAL = 50


# ============================================================
# Data Loading
# ============================================================

def load_claims(path):
    df = pd.read_csv(path)
    required = ['upstream_gene', 'downstream_gene', 'predicted_direction']
    missing = [c for c in required if c not in df.columns]
    if missing:
        log.error(f"Claims file missing required columns: {missing}")
        log.error(f"Found columns: {list(df.columns)}")
        sys.exit(1)

    df['predicted_direction'] = df['predicted_direction'].str.upper().str.strip()
    if 'cell_type_context' not in df.columns:
        df['cell_type_context'] = ''
    if 'source' not in df.columns:
        df['source'] = ''

    log.info(f"Loaded {len(df)} claims from {path}")
    return df


def load_replogle(path):
    if path and os.path.exists(path):
        df = pd.read_parquet(path)
        log.info(f"Loaded Replogle data: {len(df):,} gene-gene pairs")
        log.info(f"  Unique knockdowns: {df['knocked_down_gene'].nunique():,}")
        return df
    log.warning(f"Replogle data not found at {path} — Tier 1 validation disabled")
    return None


def load_psychencode(path):
    if path and os.path.exists(path):
        df = pd.read_parquet(path)
        log.info(f"Loaded PsychENCODE data: {len(df)} gene-cell_type entries")
        return df
    log.warning(f"PsychENCODE data not found at {path} — Tier 2 validation disabled")
    return None


def load_gene_map(path):
    if path and os.path.exists(path):
        df = pd.read_csv(path, sep='\t')
        mapping = dict(zip(df['symbol'], df['ensembl_id']))
        log.info(f"Loaded gene mapping: {len(mapping)} genes")
        return mapping
    log.warning(f"Gene map not found at {path} — Tier 3 validation disabled")
    return {}


# ============================================================
# Tier 1: Perturbation Lookup
# ============================================================

def validate_perturbation(claims, replogle_df):
    """
    For each claim, look up whether the upstream gene was knocked down
    in the Replogle dataset and what happened to the downstream gene.

    Returns a dict keyed by claim index with perturbation evidence.
    """
    if replogle_df is None:
        return {i: _empty_perturb_result("NO_DATA") for i in claims.index}

    results = {}
    np.random.seed(RANDOM_SEED)

    all_affected_genes = replogle_df['affected_gene'].unique()
    available_kd_genes = set(replogle_df['knocked_down_gene'].unique())

    for idx, row in claims.iterrows():
        upstream = row['upstream_gene']
        downstream = row['downstream_gene']
        predicted_dir = row['predicted_direction']

        if upstream not in available_kd_genes:
            results[idx] = _empty_perturb_result("UPSTREAM_NOT_TESTED")
            continue

        kd_data = replogle_df[replogle_df['knocked_down_gene'] == upstream]
        target_hit = kd_data[kd_data['affected_gene'] == downstream]

        if len(target_hit) == 0:
            # Upstream was knocked down but downstream didn't make the
            # significance cutoff (|Z| > 1 in preprocessing). This means
            # the effect was small — the downstream gene was NOT notably
            # affected by this knockdown.
            results[idx] = {
                'perturb_status': 'WEAK_OR_ABSENT',
                'perturb_z_score': 0.0,
                'perturb_direction_match': None,
                'perturb_percentile': 0.0,
                'perturb_cell_line': 'K562',
                'perturb_null_mean': None,
                'perturb_null_std': None,
                'perturb_note': (
                    f"Knockdown of {upstream} did not produce a notable "
                    f"expression change in {downstream} (below |Z|>1 threshold). "
                    f"This gene was not more affected than typical genes."
                )
            }
            continue

        observed_z = float(target_hit.iloc[0]['z_score'])
        cell_line = target_hit.iloc[0].get('cell_line', 'K562')

        # Determine observed direction
        # In CRISPRi data: negative Z means gene went DOWN after knockdown
        observed_dir = "DOWN" if observed_z < 0 else "UP"
        direction_match = (observed_dir == predicted_dir)

        # Compute null distribution: random genes' responses to same knockdown
        all_z_for_this_kd = kd_data['z_score'].values
        if len(all_z_for_this_kd) > NULL_SAMPLE_SIZE:
            null_sample = np.random.choice(all_z_for_this_kd, NULL_SAMPLE_SIZE, replace=False)
        else:
            null_sample = all_z_for_this_kd

        null_abs = np.abs(null_sample)
        observed_abs = abs(observed_z)
        percentile = float(stats.percentileofscore(null_abs, observed_abs))

        results[idx] = {
            'perturb_status': 'DATA_FOUND',
            'perturb_z_score': round(observed_z, 4),
            'perturb_direction_match': direction_match,
            'perturb_percentile': round(percentile, 1),
            'perturb_cell_line': cell_line,
            'perturb_null_mean': round(float(np.mean(null_abs)), 4),
            'perturb_null_std': round(float(np.std(null_abs)), 4),
            'perturb_note': (
                f"Knockdown of {upstream} changed {downstream} expression "
                f"(Z={observed_z:.2f}, direction={'matches' if direction_match else 'OPPOSES'} prediction, "
                f"percentile={percentile:.0f}% vs random genes after same knockdown). "
                f"Cell line: {cell_line}."
            )
        }

    return results


def _empty_perturb_result(status):
    return {
        'perturb_status': status,
        'perturb_z_score': None,
        'perturb_direction_match': None,
        'perturb_percentile': None,
        'perturb_cell_line': None,
        'perturb_null_mean': None,
        'perturb_null_std': None,
        'perturb_note': f"No perturbation data available ({status})"
    }


# ============================================================
# Tier 2: Disease Tissue Expression
# ============================================================

def validate_disease_expression(claims, psychencode_df):
    """
    Check whether downstream genes are actually dysregulated in ASD brain tissue.
    """
    if psychencode_df is None:
        return {i: _empty_de_result("NO_DATA") for i in claims.index}

    results = {}

    for idx, row in claims.iterrows():
        downstream = row['downstream_gene']
        cell_context = row.get('cell_type_context', '')

        gene_data = psychencode_df[psychencode_df['gene'] == downstream]

        if len(gene_data) == 0:
            results[idx] = _empty_de_result("GENE_NOT_IN_DATASET")
            continue

        # If user specified a cell type, try to match it
        best_match = None
        cell_type_matched = False

        if cell_context:
            context_lower = cell_context.lower().replace(' ', '_')
            for _, grow in gene_data.iterrows():
                ct = str(grow['cell_type']).lower().replace(' ', '_')
                if context_lower in ct or ct in context_lower:
                    best_match = grow
                    cell_type_matched = True
                    break

        # If no cell type match, use the most significant result
        if best_match is None:
            gene_data_sorted = gene_data.sort_values('padj')
            best_match = gene_data_sorted.iloc[0]
            cell_type_matched = (cell_context == '')

        padj = float(best_match['padj'])
        log2fc = float(best_match['log2fc'])
        ct = best_match['cell_type']
        source = best_match.get('dataset_source', 'PsychENCODE')

        is_significant = padj < 0.05
        direction = "DOWN" if log2fc < 0 else "UP"

        results[idx] = {
            'de_status': 'SIGNIFICANT' if is_significant else 'NOT_SIGNIFICANT',
            'de_log2fc': round(log2fc, 4),
            'de_padj': round(padj, 6),
            'de_cell_type': ct,
            'de_cell_type_matched': cell_type_matched,
            'de_direction': direction,
            'de_source': source,
            'de_note': (
                f"{downstream} is {'significantly' if is_significant else 'NOT significantly'} "
                f"dysregulated in ASD brain ({ct}: log2FC={log2fc:.3f}, padj={padj:.4f}). "
                f"{'Cell type matches claim context.' if cell_type_matched else f'NOTE: Data from {ct}, not {cell_context}.'} "
                f"Source: {source}."
            )
        }

    return results


def _empty_de_result(status):
    return {
        'de_status': status,
        'de_log2fc': None,
        'de_padj': None,
        'de_cell_type': None,
        'de_cell_type_matched': None,
        'de_direction': None,
        'de_source': None,
        'de_note': f"No disease tissue expression data available ({status})"
    }


# ============================================================
# Tier 3: Genetic Association (Open Targets API)
# ============================================================

def validate_genetic_association(claims, gene_map):
    """
    Query Open Targets for disease associations for each gene.
    """
    if not gene_map:
        return {i: _empty_ot_result("NO_GENE_MAP") for i in claims.index}

    import requests

    OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
    ASD_EFO = "MONDO_0005258"  # autism spectrum disorder

    results = {}
    all_genes = set(claims['upstream_gene'].tolist() + claims['downstream_gene'].tolist())
    gene_scores = {}

    for gene in all_genes:
        ensembl_id = gene_map.get(gene)
        if not ensembl_id:
            gene_scores[gene] = {'score': None, 'status': 'NO_ENSEMBL_ID'}
            continue

        query = """
        query($ensemblId: String!, $efoId: String!) {
          disease(efoId: $efoId) {
            associatedTargets(page: {size: 1, index: 0},
                             aggregationFilters: [{name: "pathwayTypes", path: "target.id", value: $ensemblId}]) {
              rows {
                target { approvedSymbol }
                score
                datatypeScores { id score }
              }
            }
          }
          target(ensemblId: $ensemblId) {
            approvedSymbol
            associatedDiseases(page: {size: 5, index: 0}) {
              rows {
                disease { id name }
                score
              }
            }
          }
        }
        """

        try:
            resp = requests.post(OT_API, json={
                'query': query,
                'variables': {'ensemblId': ensembl_id, 'efoId': ASD_EFO}
            }, timeout=10)

            if resp.status_code == 200:
                data = resp.json().get('data', {})

                # Check for direct ASD association
                asd_rows = (data.get('disease', {}) or {}).get('associatedTargets', {}).get('rows', [])
                asd_score = asd_rows[0]['score'] if asd_rows else 0.0

                # Get top disease associations for context
                target_data = data.get('target', {}) or {}
                top_diseases = []
                for drow in (target_data.get('associatedDiseases', {}) or {}).get('rows', []):
                    top_diseases.append({
                        'disease': drow['disease']['name'],
                        'score': drow['score']
                    })

                gene_scores[gene] = {
                    'score': round(asd_score, 4),
                    'status': 'FOUND',
                    'top_diseases': top_diseases[:3]
                }
            else:
                gene_scores[gene] = {'score': None, 'status': f'API_ERROR_{resp.status_code}'}

        except Exception as e:
            gene_scores[gene] = {'score': None, 'status': f'ERROR: {str(e)[:50]}'}

    # Map back to claims
    for idx, row in claims.iterrows():
        up_info = gene_scores.get(row['upstream_gene'], {'score': None, 'status': 'UNKNOWN'})
        down_info = gene_scores.get(row['downstream_gene'], {'score': None, 'status': 'UNKNOWN'})

        results[idx] = {
            'ot_upstream_asd_score': up_info.get('score'),
            'ot_downstream_asd_score': down_info.get('score'),
            'ot_upstream_status': up_info.get('status'),
            'ot_downstream_status': down_info.get('status'),
            'ot_note': (
                f"Open Targets ASD association: {row['upstream_gene']}="
                f"{up_info.get('score', 'N/A')}, {row['downstream_gene']}="
                f"{down_info.get('score', 'N/A')}"
            )
        }

    return results


def _empty_ot_result(status):
    return {
        'ot_upstream_asd_score': None,
        'ot_downstream_asd_score': None,
        'ot_upstream_status': status,
        'ot_downstream_status': status,
        'ot_note': f"Genetic association lookup unavailable ({status})"
    }


# ============================================================
# Confidence Grading
# ============================================================

def assign_confidence_grades(claims, perturb_results, de_results, ot_results):
    """
    Combine evidence from all three tiers into a confidence grade for each claim.
    """
    grades = {}

    for idx in claims.index:
        pr = perturb_results.get(idx, _empty_perturb_result("MISSING"))
        dr = de_results.get(idx, _empty_de_result("MISSING"))
        ot = ot_results.get(idx, _empty_ot_result("MISSING"))

        perturb_status = pr['perturb_status']
        percentile = pr.get('perturb_percentile') or 0
        direction_match = pr.get('perturb_direction_match')
        de_sig = dr.get('de_status') == 'SIGNIFICANT'
        de_ct_match = dr.get('de_cell_type_matched', False)

        # Grade assignment logic
        if perturb_status == 'DATA_FOUND':
            if direction_match is False:
                grade = 'CONTRADICTED'
                reason = (
                    "Perturbation data shows a significant effect in the OPPOSITE "
                    "direction to the prediction."
                )
            elif percentile >= PERCENTILE_VALIDATED and direction_match and de_sig:
                grade = 'VALIDATED'
                reason = (
                    f"Perturbation confirms predicted direction (top {100-percentile:.0f}% "
                    f"vs null), and disease tissue expression is consistent."
                )
            elif percentile >= PERCENTILE_VALIDATED and direction_match:
                grade = 'PARTIALLY_SUPPORTED'
                reason = (
                    f"Perturbation confirms predicted direction (top {100-percentile:.0f}% "
                    f"vs null), but disease tissue data is missing or non-significant."
                )
            elif percentile >= PERCENTILE_PARTIAL and direction_match:
                grade = 'PARTIALLY_SUPPORTED'
                reason = (
                    f"Perturbation effect is real but modest ({percentile:.0f}th percentile "
                    f"vs null). Direction matches prediction."
                )
            else:
                grade = 'WEAK'
                reason = (
                    f"Perturbation data exists but the downstream gene is not more "
                    f"affected than random genes ({percentile:.0f}th percentile vs null)."
                )

            # Cell type caveat overlay
            if grade in ('VALIDATED', 'PARTIALLY_SUPPORTED') and not de_ct_match and de_sig:
                grade = 'CELL_TYPE_CAVEAT'
                reason += (
                    f" NOTE: Perturbation data is from {pr.get('perturb_cell_line', 'K562')}, "
                    f"disease tissue data from {dr.get('de_cell_type', 'unknown')} — "
                    f"neither matches the claimed context."
                )

        elif perturb_status == 'WEAK_OR_ABSENT':
            if de_sig:
                grade = 'PARTIALLY_SUPPORTED'
                reason = (
                    "No notable perturbation effect, but the downstream gene IS "
                    "dysregulated in disease tissue. The regulatory relationship "
                    "may exist through indirect mechanisms not captured by single-gene knockdown."
                )
            else:
                grade = 'WEAK'
                reason = (
                    "Knockdown of the upstream gene did not notably affect the downstream "
                    "gene, and the downstream gene is not significantly dysregulated in "
                    "disease tissue."
                )

        elif perturb_status in ('UPSTREAM_NOT_TESTED', 'NO_DATA'):
            if de_sig:
                grade = 'PARTIALLY_SUPPORTED'
                reason = (
                    "No perturbation data available for the upstream gene. However, "
                    "the downstream gene IS dysregulated in disease tissue, which is "
                    "consistent with the claim (but cannot confirm the specific regulatory link)."
                )
            else:
                grade = 'UNTESTABLE'
                reason = (
                    "No perturbation data for the upstream gene and no significant "
                    "disease tissue dysregulation for the downstream gene. This claim "
                    "cannot be validated or refuted with available data."
                )
        else:
            grade = 'UNTESTABLE'
            reason = f"Insufficient data to evaluate ({perturb_status})."

        grades[idx] = {
            'confidence_grade': grade,
            'grade_reason': reason
        }

    return grades


# ============================================================
# Output Generation
# ============================================================

def build_results_table(claims, perturb_results, de_results, ot_results, grades):
    rows = []
    for idx in claims.index:
        claim = claims.loc[idx]
        pr = perturb_results.get(idx, {})
        dr = de_results.get(idx, {})
        ot = ot_results.get(idx, {})
        gr = grades.get(idx, {})

        rows.append({
            'upstream_gene': claim['upstream_gene'],
            'downstream_gene': claim['downstream_gene'],
            'predicted_direction': claim['predicted_direction'],
            'cell_type_context': claim.get('cell_type_context', ''),
            'source': claim.get('source', ''),
            # Tier 1
            'perturb_status': pr.get('perturb_status'),
            'perturb_z_score': pr.get('perturb_z_score'),
            'perturb_direction_match': pr.get('perturb_direction_match'),
            'perturb_percentile': pr.get('perturb_percentile'),
            'perturb_cell_line': pr.get('perturb_cell_line'),
            'perturb_null_mean': pr.get('perturb_null_mean'),
            'perturb_null_std': pr.get('perturb_null_std'),
            # Tier 2
            'de_status': dr.get('de_status'),
            'de_log2fc': dr.get('de_log2fc'),
            'de_padj': dr.get('de_padj'),
            'de_cell_type': dr.get('de_cell_type'),
            'de_cell_type_matched': dr.get('de_cell_type_matched'),
            # Tier 3
            'ot_upstream_asd_score': ot.get('ot_upstream_asd_score'),
            'ot_downstream_asd_score': ot.get('ot_downstream_asd_score'),
            # Grade
            'confidence_grade': gr.get('confidence_grade'),
            'grade_reason': gr.get('grade_reason'),
            # Evidence narratives
            'perturb_evidence': pr.get('perturb_note', ''),
            'de_evidence': dr.get('de_note', ''),
            'ot_evidence': ot.get('ot_note', ''),
        })

    return pd.DataFrame(rows)


def compute_base_rate(results_df, replogle_df):
    """
    What fraction of random gene pairs would score VALIDATED or PARTIALLY_SUPPORTED?
    This is the base rate comparison that makes TruthSeq useful.
    """
    if replogle_df is None:
        return None

    np.random.seed(RANDOM_SEED + 1)
    n_claims = len(results_df)
    n_simulations = 1000

    all_kd_genes = replogle_df['knocked_down_gene'].unique()
    all_affected_genes = replogle_df['affected_gene'].unique()

    positive_counts = []

    for _ in range(n_simulations):
        random_ups = np.random.choice(all_kd_genes, n_claims, replace=True)
        random_downs = np.random.choice(all_affected_genes, n_claims, replace=True)

        count = 0
        for up, down in zip(random_ups, random_downs):
            hit = replogle_df[
                (replogle_df['knocked_down_gene'] == up) &
                (replogle_df['affected_gene'] == down)
            ]
            if len(hit) > 0:
                z = abs(float(hit.iloc[0]['z_score']))
                kd_all = replogle_df[replogle_df['knocked_down_gene'] == up]['z_score'].abs()
                pct = stats.percentileofscore(kd_all.values, z)
                if pct >= PERCENTILE_PARTIAL:
                    count += 1
        positive_counts.append(count)

    return {
        'null_mean': round(np.mean(positive_counts), 2),
        'null_std': round(np.std(positive_counts), 2),
        'observed': int(results_df['confidence_grade'].isin(
            ['VALIDATED', 'PARTIALLY_SUPPORTED']).sum()),
        'n_claims': n_claims,
        'n_simulations': n_simulations
    }


def generate_summary_report(results_df, base_rate, output_dir):
    """Write the human-readable Markdown summary report."""
    grade_counts = results_df['confidence_grade'].value_counts()
    n = len(results_df)

    lines = [
        f"# TruthSeq Validation Report",
        f"",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"",
        f"## Overall Scorecard",
        f"",
        f"Claims analyzed: {n}",
        f"",
    ]

    grade_order = ['VALIDATED', 'PARTIALLY_SUPPORTED', 'CELL_TYPE_CAVEAT',
                   'WEAK', 'CONTRADICTED', 'UNTESTABLE']
    for g in grade_order:
        c = grade_counts.get(g, 0)
        pct = (c / n * 100) if n > 0 else 0
        lines.append(f"  {g}: {c}  ({pct:.1f}%)")

    lines.append("")

    # Base rate comparison
    if base_rate:
        obs = base_rate['observed']
        null_mean = base_rate['null_mean']
        null_std = base_rate['null_std']

        if null_std > 0:
            z_vs_null = (obs - null_mean) / null_std
            p_val = 1 - stats.norm.cdf(z_vs_null)
            lines.append(f"## Base Rate Comparison")
            lines.append(f"")
            lines.append(
                f"In a random set of {base_rate['n_claims']} gene pairs drawn from the same "
                f"perturbation database, the expected number scoring VALIDATED or PARTIALLY_SUPPORTED "
                f"is {null_mean} +/- {null_std}. This claim set scores {obs}, "
                f"which is {'above' if obs > null_mean else 'at or below'} the null expectation "
                f"(Z={z_vs_null:.2f}, p={p_val:.3f})."
            )
        else:
            lines.append(f"## Base Rate Comparison")
            lines.append(f"")
            lines.append(f"Null distribution had zero variance — likely too few perturbation matches to compute a reliable base rate.")

        lines.append("")

    # Per-claim details
    lines.append("## Per-Claim Evidence")
    lines.append("")

    for _, row in results_df.iterrows():
        lines.append(f"### {row['upstream_gene']} -> {row['downstream_gene']} "
                      f"(predicted: {row['predicted_direction']})")
        lines.append(f"")
        lines.append(f"**Grade: {row['confidence_grade']}**")
        lines.append(f"")
        lines.append(f"Reason: {row['grade_reason']}")
        lines.append(f"")
        if row.get('perturb_evidence'):
            lines.append(f"Perturbation: {row['perturb_evidence']}")
            lines.append(f"")
        if row.get('de_evidence'):
            lines.append(f"Disease tissue: {row['de_evidence']}")
            lines.append(f"")
        if row.get('ot_evidence'):
            lines.append(f"Genetic association: {row['ot_evidence']}")
            lines.append(f"")
        lines.append("---")
        lines.append("")

    # Warnings
    lines.append("## Warnings and Limitations")
    lines.append("")

    has_k562 = any(results_df['perturb_cell_line'] == 'K562')
    if has_k562:
        lines.append(
            "**Cell type mismatch**: All Tier 1 perturbation data is from K562 cells "
            "(chronic myeloid leukemia line). Gene regulatory relationships can differ "
            "across cell types. Claims about neuronal gene regulation should be validated "
            "in neuronal perturbation datasets when available."
        )
        lines.append("")

    lines.append(
        "**Perturbation =/= regulation**: A knockdown effect shows that gene X's "
        "expression affects gene Y, but does not prove direct transcriptional regulation. "
        "The effect could be indirect, mediated through intermediate genes or pathways."
    )
    lines.append("")

    lines.append(
        "**Observational =/= causal**: Disease tissue differential expression (Tier 2) "
        "shows correlation with disease state, not causation. A gene can be dysregulated "
        "in ASD brain as a downstream consequence rather than a causal driver."
    )

    report_path = os.path.join(output_dir, "truthseq_summary.md")
    with open(report_path, 'w') as f:
        f.write('\n'.join(lines))
    log.info(f"Summary report written to {report_path}")


def generate_heatmap(results_df, output_dir):
    """Create a confidence heatmap visualization."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        log.warning("matplotlib not available — skipping heatmap generation")
        return

    grade_colors = {
        'VALIDATED': '#2ecc71',
        'PARTIALLY_SUPPORTED': '#f1c40f',
        'CELL_TYPE_CAVEAT': '#e67e22',
        'WEAK': '#bdc3c7',
        'CONTRADICTED': '#e74c3c',
        'UNTESTABLE': '#ecf0f1',
    }

    fig, ax = plt.subplots(figsize=(12, max(4, len(results_df) * 0.6)))

    labels = [f"{r['upstream_gene']} -> {r['downstream_gene']}" for _, r in results_df.iterrows()]
    tiers = ['Perturbation\n(Tier 1)', 'Disease Tissue\n(Tier 2)', 'Overall\nGrade']

    for i, (_, row) in enumerate(results_df.iterrows()):
        # Tier 1 color
        ps = row.get('perturb_status', '')
        if ps == 'DATA_FOUND':
            pct = row.get('perturb_percentile', 0) or 0
            dm = row.get('perturb_direction_match')
            if dm is False:
                c1 = '#e74c3c'
            elif pct >= 90:
                c1 = '#2ecc71'
            elif pct >= 50:
                c1 = '#f1c40f'
            else:
                c1 = '#bdc3c7'
        elif ps == 'WEAK_OR_ABSENT':
            c1 = '#bdc3c7'
        else:
            c1 = '#ecf0f1'

        # Tier 2 color
        ds = row.get('de_status', '')
        if ds == 'SIGNIFICANT':
            c2 = '#2ecc71'
        elif ds == 'NOT_SIGNIFICANT':
            c2 = '#bdc3c7'
        else:
            c2 = '#ecf0f1'

        # Overall grade color
        grade = row.get('confidence_grade', 'UNTESTABLE')
        c3 = grade_colors.get(grade, '#ecf0f1')

        for j, color in enumerate([c1, c2, c3]):
            rect = plt.Rectangle((j, i), 0.9, 0.8, facecolor=color, edgecolor='white', linewidth=2)
            ax.add_patch(rect)

    ax.set_xlim(-0.1, 3)
    ax.set_ylim(-0.1, len(results_df))
    ax.set_xticks([0.45, 1.45, 2.45])
    ax.set_xticklabels(tiers, fontsize=10)
    ax.set_yticks([i + 0.4 for i in range(len(results_df))])
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_title("TruthSeq Confidence Heatmap", fontsize=14, fontweight='bold', pad=15)

    legend_patches = [mpatches.Patch(color=c, label=g) for g, c in grade_colors.items()]
    ax.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(1.02, 1),
              fontsize=8, frameon=True)

    plt.tight_layout()
    heatmap_path = os.path.join(output_dir, "truthseq_heatmap.png")
    plt.savefig(heatmap_path, dpi=150, bbox_inches='tight')
    plt.close()
    log.info(f"Heatmap saved to {heatmap_path}")


# ============================================================
# Main
# ============================================================

def create_demo_claims(output_path):
    """Create the demo claims file if none provided."""
    demo = pd.DataFrame([
        {"upstream_gene": "MYT1L", "downstream_gene": "MEF2C", "predicted_direction": "DOWN",
         "cell_type_context": "neuron", "source": "GRN inference"},
        {"upstream_gene": "TCF4", "downstream_gene": "MEF2C", "predicted_direction": "DOWN",
         "cell_type_context": "neuron", "source": "GRN inference"},
        {"upstream_gene": "MYT1L", "downstream_gene": "SCN2A", "predicted_direction": "DOWN",
         "cell_type_context": "excitatory_neuron", "source": "co-expression"},
        {"upstream_gene": "EP300", "downstream_gene": "MYT1L", "predicted_direction": "DOWN",
         "cell_type_context": "neuron", "source": "chromatin regulation"},
        {"upstream_gene": "FOXP1", "downstream_gene": "TCF4", "predicted_direction": "DOWN",
         "cell_type_context": "neuron", "source": "TF binding prediction"},
        {"upstream_gene": "MEF2C", "downstream_gene": "KCNA2", "predicted_direction": "DOWN",
         "cell_type_context": "neuron", "source": "target gene analysis"},
        {"upstream_gene": "MEF2C", "downstream_gene": "GRIN2B", "predicted_direction": "DOWN",
         "cell_type_context": "excitatory_neuron", "source": "target gene analysis"},
        {"upstream_gene": "MECP2", "downstream_gene": "CDKL5", "predicted_direction": "DOWN",
         "cell_type_context": "neuron", "source": "literature"},
    ])
    demo.to_csv(output_path, index=False)
    log.info(f"Created demo claims file: {output_path}")
    return output_path


def main():
    parser = argparse.ArgumentParser(description="TruthSeq: Validate gene regulatory claims")
    parser.add_argument('--claims', default=None, help='CSV file with claims (or omit for demo)')
    parser.add_argument('--replogle', default='replogle_knockdown_effects.parquet')
    parser.add_argument('--psychencode', default='psychencode_asd_de.parquet')
    parser.add_argument('--gene-map', default='gene_id_mapping.tsv')
    parser.add_argument('--output', default='truthseq_report', help='Output directory')
    parser.add_argument('--skip-ot', action='store_true', help='Skip Open Targets API queries')
    parser.add_argument('--skip-base-rate', action='store_true',
                        help='Skip base rate simulation (faster)')

    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)

    # Load claims
    if args.claims and os.path.exists(args.claims):
        claims = load_claims(args.claims)
    else:
        log.info("No claims file provided — using demo ORC gene claims")
        demo_path = os.path.join(args.output, "demo_claims.csv")
        create_demo_claims(demo_path)
        claims = load_claims(demo_path)

    # Load reference data
    replogle_df = load_replogle(args.replogle)
    psychencode_df = load_psychencode(args.psychencode)
    gene_map = load_gene_map(args.gene_map) if not args.skip_ot else {}

    # Run three tiers of validation
    log.info("")
    log.info("=== Tier 1: Perturbation Lookup ===")
    perturb_results = validate_perturbation(claims, replogle_df)

    log.info("")
    log.info("=== Tier 2: Disease Tissue Expression ===")
    de_results = validate_disease_expression(claims, psychencode_df)

    log.info("")
    log.info("=== Tier 3: Genetic Association ===")
    if args.skip_ot:
        log.info("Skipped (--skip-ot flag)")
        ot_results = {i: _empty_ot_result("SKIPPED") for i in claims.index}
    else:
        ot_results = validate_genetic_association(claims, gene_map)

    # Assign grades
    log.info("")
    log.info("=== Confidence Grading ===")
    grades = assign_confidence_grades(claims, perturb_results, de_results, ot_results)

    # Build results table
    results_df = build_results_table(claims, perturb_results, de_results, ot_results, grades)

    # Save CSV
    csv_path = os.path.join(args.output, "truthseq_results.csv")
    results_df.to_csv(csv_path, index=False)
    log.info(f"Results table saved to {csv_path}")

    # Base rate comparison
    base_rate = None
    if not args.skip_base_rate and replogle_df is not None:
        log.info("")
        log.info("=== Base Rate Simulation (1000 random gene sets) ===")
        base_rate = compute_base_rate(results_df, replogle_df)
        if base_rate:
            log.info(f"  Null expectation: {base_rate['null_mean']} +/- {base_rate['null_std']}")
            log.info(f"  Observed: {base_rate['observed']}")

    # Generate reports
    log.info("")
    log.info("=== Generating Reports ===")
    generate_summary_report(results_df, base_rate, args.output)
    generate_heatmap(results_df, args.output)

    # Print summary to console
    print("")
    print("=" * 60)
    print("TruthSeq Validation Complete")
    print("=" * 60)
    grade_counts = results_df['confidence_grade'].value_counts()
    for g in ['VALIDATED', 'PARTIALLY_SUPPORTED', 'CELL_TYPE_CAVEAT',
              'WEAK', 'CONTRADICTED', 'UNTESTABLE']:
        c = grade_counts.get(g, 0)
        print(f"  {g:25s} {c}")
    print(f"  {'TOTAL':25s} {len(results_df)}")
    print("")
    print(f"Output files in: {args.output}/")
    print(f"  truthseq_results.csv    — full evidence table")
    print(f"  truthseq_summary.md     — human-readable report")
    print(f"  truthseq_heatmap.png    — confidence visualization")
    print("")

    return 0


if __name__ == '__main__':
    sys.exit(main())
