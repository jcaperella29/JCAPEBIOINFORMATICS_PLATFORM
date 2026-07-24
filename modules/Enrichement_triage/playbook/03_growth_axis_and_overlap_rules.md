# Growth Axis, Proliferation, and Overlap Rules

Many enriched pathways fall on a shared **growth axis** that can be biologically correct in very different phenotypes.
These include PI3K–AKT–mTOR, MYC, E2F, ribosome biogenesis, translation, and cell cycle.

These pathways should NOT be auto-labeled as artifacts. Their meaning depends on the phenotype.

---

## Core Growth Axis Components

| Pathway / Program | What it reflects |
|------------------|------------------|
| PI3K–AKT–mTOR | Anabolic signaling, nutrient sensing, protein synthesis |
| MYC | Transcriptional amplification, biomass accumulation |
| E2F / Cell cycle | Proliferation, S-phase entry |
| Ribosome biogenesis | Capacity for protein synthesis |
| Translation / elongation | Actual protein output |
| Nucleotide synthesis | Growth support |

---

## When Growth Programs Are Likely **Causal Drivers**

Growth programs are likely real drivers when the phenotype involves:

- Cancer or tumor growth
- Height or developmental growth
- Muscle hypertrophy
- Tissue regeneration
- Hyperplasia or increased cell mass
- Metabolic/anabolic activation

In these contexts, PI3K, mTOR, MYC, ribosome, and cell cycle programs are often the **mechanism**, not noise.

---

## When Growth Programs Are Often **Confounders**

Growth programs are often reactive or confounding when the phenotype is:

- Inflammation without proliferation
- Immune activation
- Differentiation without expansion
- Fibrosis without cell division
- Toxicity or stress
- CRISPR or drug screens affecting viability

In these cases, growth signatures often reflect:
- Cells dying or stopping division
- Selection for fitter clones
- Global stress effects

They are downstream, not the primary biology.

---

## How to Disambiguate Growth vs Confounder

Ask:

1. **Does the phenotype imply growth?**
   - If yes → growth programs are plausible drivers.
   - If no → treat growth with skepticism.

2. **Do multiple growth layers agree?**
   - True growth: mTOR + ribosome + cell cycle + nucleotide synthesis
   - Artifact: only one layer (e.g., ribosome alone)

3. **Do growth genes appear as leading drivers?**
   - MYC, CCND1, MKI67, RPL/RPS families, EIFs

---

## Validation Experiments for Growth Programs

| Hypothesis | Experiment | Readout |
|----------|------------|--------|
| mTOR-driven growth | Rapamycin / Torin | pS6, p4EBP1 |
| Increased translation | SUnSET assay | Puromycin incorporation |
| Proliferation | EdU / Ki67 | Cell cycle entry |
| Hypertrophy | Cell size, protein content | Mass, protein/cell |
| Ribosome biogenesis | rRNA synthesis | qPCR / nucleolar markers |

Growth should be suppressed by pathway inhibition if it is causal.

---

## Key Rule

**Growth pathways can be correct in both cancer and height/muscle phenotypes.  
They must be judged relative to phenotype, not discarded.**
