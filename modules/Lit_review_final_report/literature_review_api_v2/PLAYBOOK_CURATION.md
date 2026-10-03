# Curating the biological consequence playbook

Edit `config/biology_playbook.yaml`. The classifier code should normally **not** need to change when adding biology.

## What belongs in an edge

Each entry has a `source`, `relation`, and `target`. Optional fields are:

- `aliases`: alternate wording for the source term used when matching pipeline output
- `direction`: usually `increase`, `decrease`, or another concise directional description
- `context`: where the relationship is expected to hold
- `caveat`: conditions that prevent the relationship from being treated as universal
- `evidence_note`: short curator note explaining why the edge is appropriate

Example:

```yaml
- source: "example kinase activity"
  aliases: ["EXK activity"]
  relation: "ACTIVATES"
  target: "example downstream pathway"
  direction: "increase"
  context: "stimulated mammalian cells"
  caveat: "Requires an intact downstream signaling complex."
```

## Choosing the relation

Use a **direct relation** only when the target is an immediate component, regulator, reaction, or well-established immediate downstream consequence. These relations may participate in a short path used for direct-evidence classification, but a paper still needs the same/equivalent assay and primary experimental data.

Use `CONTRIBUTES_TO_PHENOTYPE`, `MODULATES`, or `REQUIRED_FOR` when the biology is real but conditional enough that it should normally support only `mechanistic_support`.

Use `ASSOCIATED_WITH` when the relationship is correlational, epidemiologic, broad, or otherwise unsuitable for causal inference. `ASSOCIATED_WITH` is useful for liberal retrieval/context but cannot create `direct` or `mechanistic_support` status by itself.

## Curation guardrails

1. Encode the **measured consequence**, not the most dramatic downstream disease interpretation.
2. Do not turn gene expression into protein activity unless that connection is itself established and appropriately caveated.
3. Do not turn enzyme abundance into pathway flux without a caveat about substrate/cofactors/pathway control.
4. Record lineage/tissue/stimulus/genotype dependence in `context` when known.
5. Prefer several short defensible edges over one sweeping causal edge.
6. When uncertain, use a weaker relation or leave the edge out pending review.

## Validation

Run:

```bash
python validate_playbook.py
```

The loader rejects missing required fields and relationship names that are not declared in `relation_policy`.
