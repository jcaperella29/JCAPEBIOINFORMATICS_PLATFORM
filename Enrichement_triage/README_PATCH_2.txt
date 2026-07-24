Patch 2: ontology + PubMed grounding

Files to replace in your repo:
- program_summarizer.py
- pubmed_client.py

What changed:
1. Adds GLUCOCORTICOID_GR_RESPONSE so dexamethasone/glucocorticoid terms and canonical GR genes do not fall into OTHER.
2. Adds MAPK_STRESS_KINASE_ATTENUATION so DUSP/MAPK/p38/JNK terms get their own deterministic bucket.
3. Improves PubMed fallback queries for dexamethasone + airway smooth muscle + glucocorticoid contexts.

After applying, rerun the synthetic dexamethasone positive-control CSV. The ranked program table should now show GLUCOCORTICOID_GR_RESPONSE as the top driver instead of putting those terms under OTHER.
