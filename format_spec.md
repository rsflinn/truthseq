# TruthSeq Disease Expression File Format

TruthSeq accepts disease expression data in TSV, CSV, or Parquet format. The data should contain gene-level differential expression results comparing disease tissue to control.

## Required columns

| Column | Type | Description |
|---|---|---|
| gene | string | Gene symbol (e.g., FOXP1, SCN2A) |
| log2fc | float | Log2 fold change (disease vs. control). Negative = downregulated in disease |
| padj | float | Adjusted p-value (FDR-corrected) |

## Optional columns

| Column | Type | Description |
|---|---|---|
| cell_type | string | Cell type or tissue compartment (e.g., Excitatory_Neuron_L2_3, Astrocyte) |
| tissue | string | Tissue of origin (e.g., cortex, breast, liver) |
| disease | string | Disease name (e.g., autism, breast cancer) |
| source | string | Dataset source (e.g., Velmeshev2019, TCGA) |

## Example

```tsv
gene	log2fc	padj	cell_type	tissue	disease	source
SCN2A	-0.35	0.0001	Excitatory_Neuron_L2_3	cortex	autism	Velmeshev2019
KCNA2	-0.42	0.0002	Inhibitory_Neuron_PV	cortex	autism	Velmeshev2019
MEF2C	-0.25	0.003	Excitatory_Neuron_L2_3	cortex	autism	Velmeshev2019
GFAP	0.45	0.0001	Astrocyte	cortex	autism	Velmeshev2019
ACTB	0.02	0.85	Excitatory_Neuron_L2_3	cortex	autism	Velmeshev2019
```

## Column name flexibility

TruthSeq recognizes common variations in column naming. These are all equivalent:

- Gene: `gene`, `gene_symbol`, `gene_name`, `symbol`, `GeneName`
- Log2FC: `log2fc`, `log2FoldChange`, `logFC`, `log2_fc`, `lfc`
- Adjusted p-value: `padj`, `p_adj`, `fdr`, `adj_pvalue`, `qvalue`
- Cell type: `cell_type`, `celltype`, `cluster`, `cell_class`

## How to create a disease expression file

Most published single-cell or bulk RNA-seq differential expression results can be converted to this format. Common sources:

- Supplementary tables from papers (usually Excel, downloadable from journal websites)
- GEO Series Matrix files processed through DESeq2 or edgeR
- Expression Atlas experiment downloads
- Your own RNA-seq analysis pipeline output
