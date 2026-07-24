# GWAS Confounders and Interpretation Rules

GWAS does NOT measure gene expression or pathway activity.
It measures **genetic association with phenotype** through linkage disequilibrium (LD).

Enrichment reflects:
- proximity to causal variants
- tissue relevance
- developmental timing
- gene density

NOT necessarily pathway activation.

---

## Major GWAS Confounders

### Gene Density Bias
Regions with many genes (MHC, subtelomeres) produce many hits.

### LD Blocks
One causal SNP implicates many nearby genes.

### Regulatory Variants
GWAS often points to enhancers, not coding genes.

### Developmental Effects
Genes active only in embryogenesis can drive adult phenotypes.

### Pleiotropy
One gene affects multiple phenotypes.

---

## Pathway Enrichment in GWAS Means

"Genes near variants related to this phenotype are enriched in this pathway."

NOT:
"This pathway is active in adult tissue."

---

## How to Tell If a GWAS Pathway Is Real

Ask:
1. Is the pathway expressed in the **relevant tissue**?
2. Are the genes **regulatory hubs** or just neighbors?
3. Does it match known disease biology?

Immune pathways in autoimmune GWAS → real  
Cell cycle in height GWAS → often developmental  
Synapse pathways in schizophrenia → real  
Ribosome pathways → usually gene-density artifacts

---

## Validation

- eQTL colocalization
- Tissue-specific expression
- CRISPR in relevant cell types
- Mendelian disease overlap
