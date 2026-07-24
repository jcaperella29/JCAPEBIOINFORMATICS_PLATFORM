# RNA-seq Confounders and Interpretation Rules

RNA-seq measures RNA abundance, not protein activity. Many enriched pathways arise from **cell state, stress, or composition**, not the mechanism of interest.

---

## High-Frequency RNA-seq Confounders

### Cell Cycle / Proliferation
Genes: MKI67, TOP2A, CCNB1, CDK1, MCMs  
Cause: growth rate changes, selection, toxicity  
Meaning: often reflects fitness rather than phenotype

### Ribosome / Translation
Genes: RPL*, RPS*, EIF*, EEF*  
Cause: growth, stress, RNA quality, library complexity  
Meaning: capacity for protein synthesis, not necessarily pathway activation

### Mitochondrial (MT- genes)
Cause: cell death, stress, low RNA quality  
Meaning: often QC or viability signal

### Interferon / ISG
Genes: IFIT*, ISG15, OAS*, MX1  
Cause: viral mimic, contamination, dissociation stress, innate immune response  
Meaning: only real if phenotype is immune/viral

### Immediate Early Genes
Genes: FOS, JUN, EGR1  
Cause: handling, dissociation, acute stimulation  
Meaning: often artifact

---

## How to Tell If a Signal Is Real

Ask:

1. Is it **supported by many genes**?
2. Is it **consistent with phenotype**?
3. Is it **expected in this cell type**?
4. Does it form a **coherent biological program**?

Isolated enrichment with no phenotype tie is usually reactive.

---

## RNA-seq Does Not Tell You

- Kinase activity
- Phosphorylation
- Protein turnover
- Pathway flux

These require follow-up assays.

---

## RNA-seq Validation Playbook

| Program | Validation |
|--------|------------|
| IFN | ISG qPCR, pSTAT1 |
| NFκB | Cytokine ELISA, RELA nuclear localization |
| mTOR | pS6, p4EBP1 |
| UPR | ATF4, CHOP |
| ECM | Collagen IF, hydroxyproline |

RNA-seq suggests hypotheses; experiments confirm them.
