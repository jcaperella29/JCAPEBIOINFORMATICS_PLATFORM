# Epigenetic vs Transcriptional Priors (and Overlap Rules)

Different assays measure different layers of biology:

- RNA-seq: transcript abundance (output)
- ATAC-seq: chromatin accessibility (potential)
- Histone marks (ChIP/CUT&RUN): regulatory state
- TF motifs: binding potential, not activity
- Proteomics/phospho: execution layer

Interpretation must reflect the layer.

---

## Key Principle

Chromatin changes indicate **regulatory permission**.
RNA changes indicate **transcriptional output**.

You can have:
- chromatin opens without RNA change (poised)
- RNA change without chromatin change (post-transcriptional or already-accessible promoters)
- both change together (strong mechanistic confidence)

---

## Strongest Evidence Hierarchy (Same Direction)

Highest confidence for pathway/program involvement when:

1) TF motif enrichment AND
2) increased accessibility at relevant elements AND
3) RNA increases for downstream targets AND
4) protein activity supports it (if measurable)

Confidence drops if only one layer changes.

---

## Common Epigenetic Patterns That Mislead

### Poised/Primed Enhancers
ATAC opens, H3K27ac absent, RNA unchanged.
Meaning: readiness, not activation.

### Broad Promoter Opening
Global accessibility increase is common with stress/apoptosis.
Meaning: general chromatin disruption.

### Motif Enrichment Over-calling
Motif enrichment ≠ TF activation.
TF activity usually depends on phosphorylation, localization, cofactors.

---

## Disambiguation Rules by Assay

### ATAC-only pathway enrichment
Treat as "regulatory potential" until:
- matched RNA supports it
- TF binding is confirmed (CUT&RUN/ChIP)
- direction fits phenotype

### RNA-only pathway enrichment
Possible explanations:
- post-transcriptional regulation already permits transcription
- signaling pathway activity causes transcription without major chromatin change
Validate with:
- phospho markers
- TF localization
- target gene kinetics

---

## Overlapping Pathways: What They Share and How to Separate

Some programs overlap heavily in gene sets and can be misclassified.

### Growth vs Translation vs Cell Cycle
- Translation/ribosome can track growth rate without proliferation.
- Cell cycle indicates division, not just anabolism.
- mTOR can increase translation without increasing division.

Rule:
- Translation only + no E2F/cell cycle = anabolism/hypertrophy more likely
- Translation + E2F/cell cycle = proliferation/growth axis
- Cell cycle + DDR/p53 = toxicity/fitness confounder unless phenotype implies DNA damage

### Inflammation vs IFN vs Stress
Many genes sit near the boundary (chemokines, STATs).
Use the dedicated disambiguation doc to separate.

---

## Validation Menu (Cross-layer)

- ATAC motif → CUT&RUN/ChIP for TF binding
- ATAC + RNA → timecourse to confirm directionality
- Suspected signaling (mTOR, NFκB, IFN) → phospho markers
- Regulatory priming → enhancer marks (H3K27ac/H3K4me1)
