# ATAC-seq Confounders

ATAC-seq measures chromatin accessibility, not transcription.

---

## Major Confounders

### Global Chromatin Opening
Stress, apoptosis, or differentiation  
Causes many promoters to look “activated”

### Mitochondrial Reads
High MT reads → dying cells

### TF Motif Over-calling
Motif enrichment ≠ TF activity

---

## Interpretation Rules

Open chromatin:
- does NOT mean transcription increased
- may reflect chromatin remodeling or stress

TF motifs:
- require validation with ChIP-seq or CUT&RUN

---

## Validation

- RNA-seq correlation
- ChIP/CUT&RUN for TF binding
- Histone marks (H3K27ac, H3K4me3)
