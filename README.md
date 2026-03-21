# TruthSeq

Validate computational gene regulatory claims against real experimental data.

You predicted Gene X regulates Gene Y? TruthSeq checks what actually happened when Gene X was knocked down in human cells. You can also layer on disease tissue expression and genetic association data to build a multi-evidence confidence grade for each claim.

## What it does

TruthSeq takes a CSV of gene-gene regulatory predictions and validates each one against up to three independent data sources:

**Tier 1 — Perturbation data (core).** The Replogle genome-wide Perturb-seq atlas: ~11,000 single-gene CRISPR knockdowns in human K562 cells, measuring expression changes across ~8,000 genes. If your upstream gene was knocked down, TruthSeq reports what happened to your downstream gene — the Z-score, direction, and percentile rank versus all other genes affected by the same knockdown.

**Tier 2 — Disease tissue expression (optional).** Supply a differential expression dataset from disease tissue, or let TruthSeq search for one. It checks whether your downstream gene is actually dysregulated in that disease context. This adds biological relevance but can't prove causation.

**Tier 3 — Genetic association (optional).** Queries Open Targets to check whether your genes have known genetic associations with the disease of interest.

Each claim gets a confidence grade: VALIDATED, PARTIALLY_SUPPORTED, WEAK, CONTRADICTED, or UNTESTABLE.

## Quick start

```bash
# 1. Install dependencies
pip3 install scanpy anndata pandas pyarrow numpy scipy requests

# 2. Run setup (downloads ~2.7 GB Perturb-seq data from Figshare)
python3 setup.py

# 3. Validate your claims (Tier 1 only)
python3 truthseq_validate.py \
    --claims your_claims.csv \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet \
    --output my_results
```

## Adding disease context (Tier 2)

TruthSeq includes a dataset registry and search tool to help you find expression data for your disease of interest.

```bash
# Search for available datasets
python3 dataset_search.py --query "autism brain RNA-seq" --verbose

# Or just specify the disease and TruthSeq will check the registry
python3 truthseq_validate.py \
    --claims claims.csv \
    --disease "autism" \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet
```

If TruthSeq finds matching datasets in the registry but you don't have one downloaded locally, it tells you exactly what's available and how to get it. Once you have a disease expression file, supply it directly:

```bash
python3 truthseq_validate.py \
    --claims claims.csv \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet \
    --disease-expr my_disease_de.tsv \
    --output results
```

See `format_spec.md` for the expected file format.

## Finding datasets

The dataset search tool queries GEO (NCBI) and ArrayExpress (EBI) for publicly available expression datasets:

```bash
# Search broadly
python3 dataset_search.py --query "Parkinson disease brain RNA-seq"

# Search with full download instructions
python3 dataset_search.py --query "breast cancer differential expression" --verbose

# Add new results to the local registry
python3 dataset_search.py --query "schizophrenia single-cell" --update-registry

# Search local registry only (no API calls)
python3 dataset_search.py --query "perturbation" --registry-only
```

The registry (`dataset_registry.csv`) ships with curated entries for major disease areas and perturbation datasets. You can add your own entries or update it automatically by running searches with `--update-registry`.

## Claims file format

CSV with these columns:

| Column | Required | Description |
|--------|----------|-------------|
| upstream_gene | yes | Gene symbol of the predicted regulator |
| downstream_gene | yes | Gene symbol of the predicted target |
| predicted_direction | yes | UP or DOWN (predicted effect of upstream on downstream) |
| cell_type_context | no | Cell type where the regulation is predicted to occur |
| source | no | Where the prediction came from |

Example:
```csv
upstream_gene,downstream_gene,predicted_direction,cell_type_context,source
MEF2C,SCN2A,DOWN,excitatory_neuron,GRN_inference
TCF4,KCNA1,DOWN,inhibitory_neuron,co-expression
GATA1,HBB,DOWN,,published_literature
```

## How grading works

TruthSeq combines evidence across tiers:

- **VALIDATED**: Perturbation confirms direction (top 10% effect vs null) AND disease tissue expression is consistent.
- **PARTIALLY_SUPPORTED**: Perturbation confirms direction but effect is modest, OR perturbation data unavailable but disease tissue evidence supports the claim.
- **WEAK**: Perturbation was tested but downstream gene wasn't notably affected, or evidence is insufficient.
- **CONTRADICTED**: Perturbation shows a significant effect in the opposite direction from the prediction.
- **UNTESTABLE**: Upstream gene wasn't in the Perturb-seq atlas and no other evidence available.

## Important caveats

The Replogle atlas uses K562 cells (chronic myeloid leukemia line). Regulatory relationships specific to other cell types (neurons, immune cells, etc.) may not show up. A WEAK grade doesn't mean the relationship is wrong — it means it wasn't detectable in this system. Disease tissue expression (Tier 2) can help fill this gap for cell-type-specific claims.

## Contributing datasets

If you know of a publicly available expression dataset that should be in the registry, add a row to `dataset_registry.csv` with the accession, description, and download instructions, then submit a pull request.

## Data sources

- Replogle et al. 2022. "Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq." *Cell* 185, 5689-5710. [Figshare](https://figshare.com/articles/dataset/Replogle_GWPS/19968745)
- Open Targets Platform: platform.opentargets.org
- NCBI GEO: ncbi.nlm.nih.gov/geo
- EMBL-EBI ArrayExpress: ebi.ac.uk/biostudies/arrayexpress
