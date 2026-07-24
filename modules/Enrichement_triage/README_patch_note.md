# PubMed relevance-label patch

This patch keeps PubMed retrieval permissive, but adds honesty labels so weak hits are not treated as direct evidence.

## Replace these files

- `pubmed_client.py`
- `reasoner.py`

## Behavior

- PubMed still returns papers when possible.
- Each paper gets:
  - `relevance_label`: `direct_support_candidate`, `general_support`, or `weak_background`
  - `relevance_score`
  - `relevance_reason`
  - `literature_use`
- The PubMed block also gets:
  - `retrieval_quality`
  - `retrieval_quality_reason`
  - `literature_use_guidance`
  - `relevance_counts`

## Why

This avoids the two bad extremes:

1. too strict -> no papers -> the app feels empty
2. too permissive -> random papers look like evidence

The app now keeps the papers but tells GPT and users how strongly to treat them.
