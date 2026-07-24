# DNA Methylation (5mC) Assay Confounders & Interpretation Rules

This playbook applies when the assay is **DNA methylation** (5mC; e.g., WGBS, RRBS, EPIC/450K arrays, targeted bisulfite panels, cfDNA methylation).  
Goal: distinguish **true regulatory methylation shifts** from **cell mixture, technical artifacts, and epigenetic memory**.

---

## 0) First principles (must obey)
1. **Methylation is not expression.** A methylation change is evidence of altered regulatory state or cell composition, not a guaranteed transcriptional effect.
2. **Cell-type composition is the #1 confounder** in bulk methylation. Treat it as guilty until proven otherwise.
3. Interpret changes by **genomic context**: promoter CpG islands vs shores/shelves, enhancers, gene bodies, repeats, imprinted loci.
4. Use **directionality carefully**:
   - Promoter CpG island **hypermethylation** often correlates with **repression**, but not always.
   - Gene body methylation can correlate with **active transcription** or splicing regulation, and is not simply “silencing.”
5. Many signals reflect **epigenetic memory** (prior differentiation, inflammation, aging, smoking) more than the acute perturbation.

---

## 1) Ask these metadata questions before interpreting
- Assay type: **WGBS / RRBS / array / targeted / cfDNA?**
- Sample type: bulk tissue? sorted cells? single-cell methylation? cfDNA?
- Case/control structure: paired? batch? time series?
- Covariates available: age, sex, ancestry, smoking, BMI, meds, inflammation markers, tumor purity, cell counts.
- Processing: bisulfite conversion method, library prep kit, array processing center.

If any of this is unknown, say so and limit claims.

---

## 2) Dominant confounders (bulk methylation)
### A) Cell-type mixture / purity (most important)
**Problem:** different cell types have distinct baseline methylomes; shifts in proportions look like differential methylation.

**Red flags**
- DMRs align strongly to known cell identity markers.
- Many changes map to immune-lineage loci in tissue experiments.
- Tumor samples show widespread hypomethylation + focal hypermethylation patterns but purity is unknown.

**Mitigations**
- Do deconvolution (reference-based if possible; otherwise reference-free).
- Include cell-type fractions/purity as covariates.
- Validate in sorted populations or single-cell methylation.

**Follow-ups**
- Flow cytometry / FACS counts.
- scRNA-seq or CITE-seq to estimate composition; compare to methylation signal.
- Marker-based qPCR / immunohistochemistry.

### B) Age / smoking / inflammation / metabolic state
**Problem:** large, reproducible methylation shifts occur with age (epigenetic clocks), smoking, chronic inflammation, obesity, diabetes, stress hormones.

**Red flags**
- Strong enrichment for known smoking-associated CpGs (AHRR, etc.) or clock-like patterns.
- Effects correlate with systemic inflammatory markers (CRP, IL6) or treatment duration.

**Mitigations**
- Include covariates; stratify; sensitivity analysis.
- Compare against known methylation signature sets (smoking, age, inflammation).

### C) Sex / X-inactivation / imprinting
**Problem:** sex chromosomes and imprinting regions generate large methylation differences unrelated to the main phenotype.

**Rules**
- Handle sex as a covariate.
- Treat X/Y signals separately; don’t let them dominate the narrative.
- Imprinted loci are special: interpret cautiously and confirm.

---

## 3) Technical & preprocessing confounders (assay-specific)

### A) Bisulfite conversion efficiency (WGBS/RRBS/targeted)
**Problem:** incomplete conversion inflates apparent methylation; over-conversion/harsh conditions can bias coverage.

**Red flags**
- Global methylation shifts without plausible biology.
- Low conversion controls / weird non-CpG methylation levels (if tracked).

**Mitigations**
- Use conversion spike-ins/controls when available.
- QC: mapping rate, duplication, coverage uniformity, CpG coverage distribution.

### B) Coverage depth & missingness (sequencing methylation)
**Problem:** low depth increases noise; RRBS has biased representation (CpG-rich regions).

**Rules**
- Don’t over-interpret single CpGs at low coverage.
- Prefer region-level (DMR) calls with consistent direction and adequate depth.

### C) Array-specific biases (EPIC/450K)
**Common issues**
- Probe cross-reactivity and SNPs at probe sites
- Dye bias / batch effects / position effects
- Sex chromosome probes can dominate

**Mitigations**
- Standard normalization and probe filtering (cross-reactive, SNP-affected probes).
- Include batch/plate/position covariates.
- Check PCA/UMAP of methylation beta/M-values for batch structure.

### D) PCR bias / GC bias / mapping artifacts
**Problem:** biases around repeats, segmental duplications, and mappability.

**Rules**
- Be cautious with DMRs dominated by low-mappability regions or repeats.
- Prefer well-mapped regions for mechanistic claims.

### E) Hydroxymethylation (5hmC) confounding
**Problem:** standard bisulfite cannot distinguish 5mC vs 5hmC (tissue-dependent, especially brain).

**Rules**
- If tissue likely has high 5hmC, explicitly mention ambiguity.
- If interpretation hinges on 5mC vs 5hmC, recommend oxBS/TAB-seq or similar.

---

## 4) Genomic context interpretation rules (use these heuristics)
When reporting “drivers,” always specify context:

### Promoters / CpG islands
- **Hypermethylation at promoter CpG islands** can suggest repression of the associated gene.
- But: if promoter is already unmethylated across conditions, small changes may not matter.

### Enhancers / distal regulatory
- Methylation changes at enhancers can indicate altered cell state or transcription factor occupancy.
- Strongly consider cell composition and activation state.

### Gene bodies
- Not simple repression; can correlate with active transcription and isoform regulation.
- Avoid claiming “silencing” from gene body hypermethylation alone.

### Repeats / global methylation
- Global hypomethylation can reflect proliferation, tumor biology, stress, or technical bias.
- Treat global shifts as supportive context, not a single-cause explanation.

---

## 5) Causality rules (do not overclaim)
You may propose mechanisms, but you must label them as **hypotheses** unless supported by orthogonal evidence.

### Stronger evidence requires:
- concordant gene expression (RNA-seq) changes
- chromatin accessibility (ATAC) or histone marks (ChIP-seq/CUT&Tag)
- TF occupancy evidence
- perturbation experiments (DNMT/TET inhibition or CRISPR epigenetic editing)

---

## 6) What “good results” look like (sanity checks)
Prefer interpretations where:
- DMRs cluster in coherent pathways *and* are not entirely explained by cell mixture.
- Region-level signals replicate across cohorts/batches.
- Effects are consistent with known biology of the phenotype *and* sample type.

If results are dominated by batch/cell-type signals, say so.

---

## 7) Follow-up experiment menu (recommend with readouts + controls)
Pick 3–6 high-value follow-ups:

### Validate methylation
- Targeted bisulfite amplicon sequencing on top DMRs (technical validation).
- Replicate in an independent cohort/batch.

### Disentangle cell composition
- FACS sorting of key cell types + methylation profiling.
- Deconvolution + sensitivity analysis including estimated fractions.

### Link to function
- RNA-seq for expression correlation.
- ATAC-seq / CUT&Tag for regulatory state.
- Reporter assays for enhancer candidates (if feasible).
- CRISPR-dCas9 DNMT3A/TET1 editing at top loci (mechanistic).

### If cfDNA
- Confirm tissue-of-origin signatures.
- Control for fragmentation and total cfDNA load.

---

## 8) Output style requirements
In the final answer:
- Always name the assay as DNA methylation.
- Call out the **top 2–4 likely confounders** explicitly.
- Separate **likely drivers** vs **likely reactive** vs **likely artifacts/confounders**.
- Tie any mechanistic claims to genomic context (promoter/enhancer/gene body).
- Provide follow-ups that test *both* biology and confounding explanations.
