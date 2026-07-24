# Tissue and Cell-Type Prior Map

Certain pathways are biologically plausible only in specific tissues.

These priors help separate signal from nonsense.

---

## Muscle (Skeletal / Cardiac)

Plausible programs:
- mTOR / translation
- Mitochondrial biogenesis
- Sarcomere assembly
- Calcium handling
- ECM remodeling (injury, hypertrophy)

Red flags:
- Synapse signaling
- Antigen presentation
- Neurotransmitter pathways

---

## Tumor / Cancer Cells

Plausible programs:
- Cell cycle
- DNA damage / p53
- PI3K–AKT–mTOR
- Hypoxia / angiogenesis
- EMT / ECM
- Metabolic rewiring

Red flags:
- Neuronal synapse
- Muscle contraction
- Cilia unless specific tumor type

---

## Immune Cells

Plausible programs:
- Cytokines (IL, TNF)
- IFN / antiviral
- Antigen presentation
- T cell receptor, BCR
- Metabolic switching (glycolysis vs OXPHOS)

Red flags:
- Sarcomere
- ECM (unless stromal)
- Neurotransmission

---

## Brain / Neurons

Plausible programs:
- Synaptic transmission
- Ion channels
- Neurotransmitter metabolism
- Vesicle trafficking
- Mitochondria

Red flags:
- Collagen / ECM
- Immune IFN (unless neuroinflammation)
- Cell cycle (unless tumor or development)

---

## Rule

If a pathway is enriched but violates tissue plausibility, treat it as:
- composition change
- developmental artifact
- mapping noise
until proven otherwise.
