# Perturb-seq (CRISPR-scRNA) Confounders

Perturb-seq is biased by **fitness and survival effects**.

---

## Major Confounders

### Cell Cycle & Growth
Knockouts that slow growth appear depleted  
Pathways enriched reflect **selection**, not phenotype

### p53 / DDR
Many sgRNAs induce DNA damage  
p53, apoptosis, and cell cycle arrest are common artifacts

### Stress & IFN
Cas9, viral vectors, and transduction trigger IFN

---

## Interpretation Rules

If a perturbation:
- reduces cell number
- induces p53
- causes cell cycle arrest

then enrichment reflects **viability**, not target pathway.

---

## Validation

- Use non-targeting controls
- Measure sgRNA dropout
- Compare to CRISPRi/a where possible
