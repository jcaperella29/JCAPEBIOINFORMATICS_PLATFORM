Patch 3: report cleanup + duplicate-claim control

Replace these files in your repo:

  schemas.py
  summarizer.py

What changed:

1. Claim de-duplication
   - normalize_claims() now keeps one primary claim per program.
   - If GPT emits a weaker conflicting claim for the same program, it is preserved as an explicit alternative explanation.
   - Handling/FOS/JUN/IEG caveats are relabeled as:
     Alternative explanation: handling-stress / immediate-early gene artifact (<PROGRAM>)

2. Program table ordering
   - PDF program table now sorts by claim priority first:
     Likely driver > Likely reactive > Uncertain > Likely artifact/confounded
     Stronger > Moderate > Weak
   - Deterministic program score is used as the tie-breaker.

3. Validation cleanup
   - Duplicate validation sections with the same heading are suppressed.
   - Alternative explanations keep their alternative heading, so they remain visible without looking like repeated program sections.

4. Synthetic positive-control notice
   - If any term contains "synthetic positive control", the PDF adds a visible demo-input notice.
   - This helps avoid accidentally treating synthetic test reports as biological evidence.

Expected result on the synthetic airway dexamethasone CSV:

  GLUCOCORTICOID_GR_RESPONSE should appear above ECM_FIBROSIS in the ranked program table.
  Duplicate MAPK claims should either collapse or the weaker FOS/JUN artifact claim should be labeled as an alternative explanation.
