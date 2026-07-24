# Single-Cell RNA-seq Confounders

scRNA-seq measures cell-level expression but is extremely sensitive to stress, dissociation, and composition shifts.

---

## Major Confounders

### Dissociation Stress
Genes: FOS, JUN, EGR1, HSPs  
Cause: enzymatic digestion, mechanical stress  
Meaning: not biology of interest

### Interferon Response
Genes: IFIT, ISG15, OAS, MX  
Cause: dissociation, dying cells, ambient RNA  
Meaning: often artifact unless immune context

### Mitochondrial Inflation
Cause: dying cells  
Meaning: low-quality or apoptotic cells

### Cell Type Proportion Shifts
If fibroblasts expand → ECM genes appear enriched  
If immune cells expand → cytokines appear enriched  
Meaning: composition change, not regulation

---

## How to Tell Regulation vs Composition

- Are marker genes of a cell type driving the signal?
- Does the same program appear **within** a cell type?

If only cell type markers change → composition  
If pathway genes change inside one cell type → regulation

---

## Validation

- Compare per-cell type DE, not pooled
- Validate with flow cytometry or cell fraction counts
