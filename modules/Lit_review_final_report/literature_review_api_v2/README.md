# JCAP Literature Review API v2

This rebuild keeps the original role of the literature/final-report stage but changes the evidence logic.

## Sources

- PubMed
- Europe PMC
- OpenAlex

Semantic Scholar is intentionally omitted.

## Design principle

**Retrieval is deliberately liberal. Evidence classification is strict.**

The search layer now runs both assay-aware and assay-free queries. It can also search nearby biological consequences from the internal playbook. A broad retrieval hit does **not** earn a stronger evidence label by itself.

## Biological consequence playbook

`config/biology_playbook.yaml` is the editable, curated biological knowledge file. `litreview/biology_playbook.py` only loads, validates, and traverses it. This keeps biological curation separate from classifier code and makes every relationship auditable.

Relationships include:

- `DIRECT_COMPONENT`
- `DIRECT_REGULATOR`
- `CATALYZES`
- `ACTIVATES`
- `INHIBITS`
- `UPSTREAM_OF`
- `DOWNSTREAM_OF`
- `CAUSES_PHENOTYPE`
- `CONTRIBUTES_TO_PHENOTYPE`
- `MODULATES`
- `REQUIRED_FOR`
- `ASSOCIATED_WITH` (context/retrieval only; never enough by itself for direct evidence)

Each edge can include direction, context and caveats. The seed graph is intentionally conservative and is meant to grow as JCAP encounters validated biology. Each edge can carry aliases, direction, context, caveats, and curator notes. See `PLAYBOOK_CURATION.md`.

Examples of the intended reasoning:

- cell-cycle activation -> cell proliferation
- proliferation -> tumor growth, but with an explicit warning that proliferation alone is not proof of cancer causation
- enzyme abundance -> enzyme activity, with the caveat that expression does not guarantee activity
- enzyme activity -> enzyme-dependent reaction -> enzyme-dependent pathway, subject to substrate/cofactor/pathway-control limitations

## Direct evidence rule

A paper can be classified as `direct` only when:

1. It uses the same assay type or an explicitly equivalent assay; and
2. It examines either:
   - the same biological target/question,
   - the same named pathway/phenotype, or
   - an endpoint connected by an allowed **<=2-edge** direct causal/mechanistic path in the biological playbook; and
3. It is primary experimental literature rather than a review.

Context mismatch does not automatically erase direct status, but it adds an explicit caution.

A <=3-edge pathway can support `mechanistic_support`, but cannot become `direct` unless the direct rule is independently satisfied.

Other evidence buckets are:

- `mechanistic_support`
- `contextual_support`
- `background`
- `exclude`

The API returns the playbook path and its relationship labels so the final report can show *why* a paper was considered biologically close rather than relying on an opaque similarity score.

## Ranking vs. evidence

Ranking is intentionally separate from evidence class. A paper can be retrieved and rank reasonably well because it is biologically related while still being only contextual or background evidence. The classifier must not upgrade it merely because retrieval similarity is high.

## Run

```bash
pip install -r requirements.txt
uvicorn litreview.app:app --host 0.0.0.0 --port 8000
```

Then POST `example_request.json` to `/review`.

## Next pieces to add

1. Fetch richer PubMed abstracts when only summary metadata are returned.
2. Move assay equivalence into editable configuration rather than hard-coded aliases.
3. Expand `config/biology_playbook.yaml` with curated pathway-to-pathway, enzyme-to-reaction, regulator-to-target and phenotype relationships encountered in validated JCAP runs.
4. Add structured evidence extraction from full abstracts.
5. Add final PI-facing synthesis from upstream quantitative/network/triage data plus classified evidence.
6. Emit JSON, Markdown and PDF outputs for Nextflow handoff.
7. Require the report generator to preserve evidence class, causal path and cautions verbatim rather than allowing an LLM to silently upgrade evidence strength.
