# Evidence Weighting and Translational Interpretation Rules

This system must distinguish between:
- observation
- interpretation
- causal inference
- translational relevance

These are NOT equivalent.

Do not overstate conclusions.

---

## Core Principle

Enrichment results generate **hypotheses**, not confirmations.

RNA, ATAC, methylation, and similar assays provide **indirect evidence** of biology.

---

## Evidence Strength Classification

For each major biological program, assign one:

### Stronger Evidence
Use stronger language ONLY when:
- Supported by many concordant genes
- Forms a coherent biological program
- Consistent with phenotype and context
- Expected in the given cell type
- Supported by multiple related pathways
- Not easily explained by confounders

Language:
- “strongly consistent with”
- “well-supported by the data”

---

### Moderate Evidence
Use cautious language when:
- Signal is coherent but indirect
- Interpretation relies partly on prior biology
- Multiple explanations remain plausible

Language:
- “suggestive of”
- “consistent with but not definitive”

---

### Weak Evidence
Use skeptical language when:
- Driven by few genes
- Broad or generic pathway labels
- Could be explained by confounders
- Conflicts with phenotype or biology
- Infers activity from RNA alone when protein is required

Language:
- “weak support for”
- “likely reflects secondary effects”
- “may be confounded”

---

## Required Classification

Each program must be labeled as one:

- **Likely driver**
- **Likely reactive**
- **Likely artifact / confounded**
- **Uncertain**

---

## Cell-Intrinsic vs Composition Effects

Always consider whether signals reflect:

- Cell type composition changes
- Immune infiltration
- Stromal contamination
- Cell death / viability differences
- Dissociation or handling stress
- Batch effects or RNA quality

If composition is plausible, state it explicitly.

---

## Assay Limitations Reminder

Do NOT infer pathway activation when the assay cannot measure it.

Examples:

- RNA-seq cannot measure:
  - kinase activity
  - phosphorylation
  - protein localization
  - pathway flux

Treat these as hypotheses requiring validation.

---

## Causal vs Reactive Distinction

Distinguish:

- **Driver** → contributes to phenotype
- **Reactive** → consequence of phenotype
- **Artifact** → technical or confounded signal

If causality is unclear, say so.

---

## Translational Relevance Rules

Do NOT imply clinical relevance unless supported.

Stronger translational claims require:
- Human data
- Known disease relevance
- Mechanistic plausibility
- Reproducibility
- Clear biomarker or intervention link

Otherwise, label as:
- “biologically interesting but not yet actionable”

---

## Required Output Tone

Prefer:
- “consistent with”
- “suggests”
- “plausibly reflects”
- “cannot distinguish from”
- “requires validation”

Avoid:
- “proves”
- “confirms”
- “demonstrates causally”

---

## Follow-Up Experiment Expectations

All proposed experiments must include:

- Specific readout (e.g., pSTAT1, ELISA, IF staining)
- Appropriate controls
- Clear purpose (test causality vs validation)

---

## Final Sanity Check

Before concluding, ask:

- Could this be explained by stress or QC?
- Could this be explained by cell composition?
- Is this directly measured or inferred?
- Does this match phenotype and biology?

If not, downgrade confidence.

---
