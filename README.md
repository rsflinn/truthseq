# TruthSeq

Most findings in biology don't survive replication. TruthSeq is designed to help you determine if your data is a true signal that can be found in human tissue, or just a computational artifact or hallucination.

A 2016 Nature survey of more than 1,500 scientists found that over 70% had tried and failed to reproduce another researcher's results. Computational biology tools are powerful enough to find patterns in any dataset, but that pattern could be noise or an artifact of the specific dataset that wouldn't appear in a living organism.

This problem is only going to increase now that more genomic datasets are published online and made open source. Combined with powerful AI tools that can run computational analyses in minutes that would have taken humans months or years to calculate, the potential for false positives will continue to plague scientific research.

You can run a perfectly valid analysis, get a statistically significant result, and still be wrong. TruthSeq aims to provide that check.

## What it does

TruthSeq validates claims about which genes control other genes — a question that comes up across cancer research, neuroscience, immunology and most other areas of genomics. If your computational analysis produced a finding about gene regulation, TruthSeq checks it against real experimental data from human cells.

You give it a list of predictions and it checks each one against up to three layers of independent evidence:

**Tier 1 — Lab knockdown data (core).** The Replogle Perturb-seq atlas: researchers at the Broad Institute knocked out nearly every human gene one by one (~11,000 genes) in human cells and measured what happened to the rest of the genome. If you claim Gene X controls Gene Y, TruthSeq looks up what actually happened to Gene Y when Gene X was disabled. This directly tests cause and effect.

**Tier 2 — Disease tissue expression (optional but recommended).** Is the target gene actually changed in the disease or tissue you're studying? TruthSeq searches a registry of publicly available gene expression datasets and two major repositories (GEO and ArrayExpress) for published data from real patient samples. This adds biological context specific to your research area.

**Tier 3 — Genetic association (optional).** Do genetic variants near these genes associate with the disease? Queries Open Targets, a public database linking genes to diseases through population genetics.

Each prediction gets a grade: VALIDATED, PARTIALLY_SUPPORTED, WEAK, CONTRADICTED, or UNTESTABLE.

### Why the tiers matter

Tier 1 is the strongest evidence — it directly tests whether knocking out one gene changes another. But it comes from one cell type (K562, a blood-derived cell line), because that's the only cell type where this experiment has been done at genome scale. A relationship that's real in brain cells or immune cells might not show up in K562.

That's why Tier 2 exists. If you're studying Parkinson's disease, TruthSeq can search for published gene expression data from brain tissue and check whether your target genes are actually disrupted there. Tier 1 tells you whether the cause-and-effect relationship is real in any human cell. Tier 2 tells you whether it's relevant to the disease or tissue you care about. Used together, they're a much stronger filter than either one alone.

## Who this is for

This tool was built by a citizen scientist who spent months chasing computational patterns that looked compelling on paper but didn't hold up under scrutiny. If you're someone who:

- Uses AI tools or bioinformatics software to explore genetic data
- Found something interesting in a computational analysis and wants to know if it's real
- Is new to genomics and wants a sanity check before going too far down a rabbit hole
- Is a working researcher who wants a fast, independent check on gene regulatory predictions

TruthSeq is for you. No biology degree required. If you can run a few terminal commands and make a spreadsheet, you can use it.

### What it can and can't test

TruthSeq tests whether one gene controls another gene's expression. That's a specific question, but it's central to a wide range of research. Cancer biologists studying whether a tumor suppressor regulates downstream targets, neuroscientists modeling transcription factor networks, immunologists tracing cytokine signaling pathways — all of these involve gene-to-gene regulatory claims that TruthSeq can check.

It can't test broader biological questions like "does this mutation cause drug resistance?" or "does this pathway drive inflammation?" Those involve protein function, cell behavior and other mechanisms beyond gene expression. But if the upstream step in your model is "Gene X controls Gene Y," TruthSeq can tell you whether that step holds up.

## What you need

- Python 3.8 or later (type `python3 --version` in your terminal to check)
- About 500 MB of free disk space for the experimental dataset
- A CSV file of gene regulatory predictions to test (or use the included example to try it out)

## Setup

```bash
# Install the Python libraries TruthSeq needs
pip3 install scanpy anndata pandas pyarrow numpy scipy requests

# Download the experimental dataset (~357 MB from Figshare, a public data repository)
python3 setup.py
```

The setup script downloads the Replogle Perturb-seq atlas from Figshare (~357 MB). This is real lab data, not a model or simulation. It takes a few minutes depending on your connection.

## Your first validation

TruthSeq reads a simple CSV where each row is one prediction: "I think this gene regulates that gene in this direction." Here's what the file looks like:

```csv
upstream_gene,downstream_gene,predicted_direction,source
SLC30A1,MT2A,DOWN,known_biology
GATA1,TYROBP,DOWN,known_biology
CHMP6,SOD2,DOWN,my_analysis
```

Quick translation:
- **upstream_gene** — the gene you think is the regulator (the one doing the controlling)
- **downstream_gene** — the gene you think is being controlled
- **predicted_direction** — when the regulator is active, does the target go UP or DOWN?
- **source** — where the prediction came from (optional, just for your notes)

Run it:

```bash
python3 truthseq_validate.py \
    --claims your_claims.csv \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet \
    --output my_results
```

## Try the example first

The included `example_claims.csv` has 11 test predictions: some textbook-true relationships, some deliberately wrong ones, and some where the direction is flipped. It's designed to show you what each grade means in practice.

```bash
python3 truthseq_validate.py \
    --claims example_claims.csv \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet \
    --output example_results
```

What to expect: SLC30A1 → MT2A (zinc transporter knockout triggers metallothionein response) and GATA1 → TYROBP (blood cell master regulator derepresses myeloid genes) should score VALIDATED. The wrong-direction controls should come back CONTRADICTED. BRCA1 and RB1 claims will score WEAK — real genes, but their effects are too diffuse to stand out in this dataset. TP53 isn't in the knockdown data at all, so it returns UNTESTABLE.

## How the grades work

- **VALIDATED** — The lab data confirms it. When the regulator was knocked out, the target gene changed significantly in the predicted direction.
- **PARTIALLY_SUPPORTED** — The direction matches but the effect was modest, or the knockdown data wasn't available but other evidence lines up.
- **WEAK** — The regulator was tested but the target gene didn't respond much. This doesn't mean the relationship is definitely wrong (see caveats below).
- **CONTRADICTED** — The lab data shows the opposite of what was predicted.
- **UNTESTABLE** — The regulator gene wasn't in the experimental dataset.

## Adding disease context

The basic validation (Tier 1) tells you what happened in a lab dish. For most research, you'll also want Tier 2: is the target gene actually disrupted in the disease or tissue you're studying?

```bash
# Search for publicly available datasets for your disease of interest
python3 dataset_search.py --query "Parkinson disease brain RNA-seq" --verbose

# Or let TruthSeq search its built-in registry
python3 truthseq_validate.py \
    --claims claims.csv \
    --disease "breast cancer" \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet
```

TruthSeq maintains a registry of publicly available gene expression datasets and searches two major repositories (GEO and ArrayExpress) for more. If it finds data you don't have locally, it tells you where to get it. A weekly automated scan adds newly published datasets to the registry.

You can also supply a disease expression file directly:

```bash
python3 truthseq_validate.py \
    --claims claims.csv \
    --replogle replogle_knockdown_effects.parquet \
    --replogle-stats replogle_knockdown_stats.parquet \
    --disease-expr my_disease_data.tsv \
    --output results
```

See `format_spec.md` for the expected file format. TruthSeq is flexible with column names — `log2fc`, `logFC`, `log2FoldChange` all work.

## Important caveats

The Tier 1 lab dataset uses K562 cells, a blood-derived cell line chosen because it's the only cell type where genome-scale knockdown experiments have been done. No equivalent dataset exists for brain cells, immune cells or any other tissue — this is the only one of its kind on the planet. Gene regulation varies by cell type, so a relationship that's real in brain cells might not show up in K562. A WEAK grade means "not detectable in this system," not "this is wrong." Tier 2 disease tissue data helps compensate by checking whether the genes are actually disrupted in the tissue relevant to your work.

TruthSeq can't prove that Gene X controls Gene Y in a living organism. It can tell you whether disabling Gene X changed Gene Y's expression in a controlled lab setting, and whether Gene Y is disrupted in real disease tissue. The point is to catch findings that are clearly unsupported before you invest months following them.

## Finding and adding datasets

```bash
# Search public repositories
python3 dataset_search.py --query "schizophrenia single-cell"

# Search with download instructions
python3 dataset_search.py --query "breast cancer differential expression" --verbose

# Add results to the local registry
python3 dataset_search.py --query "epilepsy RNA-seq" --update-registry

# Search local registry only (no API calls)
python3 dataset_search.py --query "perturbation" --registry-only
```

The registry (`dataset_registry.csv`) ships with seed entries and grows automatically. You can add entries by hand or submit a pull request.

## Claims file format (reference)

| Column | Required | Description |
|--------|----------|-------------|
| upstream_gene | yes | Gene symbol of the predicted regulator |
| downstream_gene | yes | Gene symbol of the predicted target |
| predicted_direction | yes | UP or DOWN |
| cell_type_context | no | Cell type where the regulation is predicted |
| source | no | Where the prediction came from |

## Data sources

- **Replogle et al. 2022.** "Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq." *Cell* 185, 5689-5710. [Figshare data](https://figshare.com/articles/dataset/Replogle_GWPS/19968745)
- **Open Targets Platform**: [platform.opentargets.org](https://platform.opentargets.org)
- **NCBI GEO**: [ncbi.nlm.nih.gov/geo](https://ncbi.nlm.nih.gov/geo)
- **ArrayExpress**: [ebi.ac.uk/biostudies/arrayexpress](https://ebi.ac.uk/biostudies/arrayexpress)

## Contributing

Found a public gene expression dataset that belongs in the registry? Add a row to `dataset_registry.csv` and submit a pull request. Bugs and feature ideas go in Issues.
