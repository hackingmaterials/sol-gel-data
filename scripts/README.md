# scripts/

Regenerates every numbered table and figure in the manuscript and the Supplementary
Information.

## Setup

From the repository root:

```bash
unzip sol_gel_dataset.jsonl.zip
pip install -r requirements.txt
jupyter notebook scripts/figures_tables.ipynb
```

Run the cells in order — later cells reuse helper functions and data loaded by earlier ones.

## Contents

- `figures_tables.ipynb` — one section per manuscript item. Table 2, Figure 2 and Figure 3
  are built from `sol_gel_dataset.jsonl`; Tables 3–5, Tables S2–S6 and Figures S1–S4 from
  the validation labels below.
- `validation/` — the human and LLM annotations for the 200-paper validation set:
  paper-level classification labels (`classification_human.jsonl`,
  `classification_LLM.jsonl`), the 286 records extracted by Gemini 3.0 Flash
  (`validation_LLM.jsonl`) and the 279 ground-truth records they are matched against
  (`validation_human.jsonl`, line-aligned with the LLM file), and
  `validation_paper_level_metrics_human.jsonl`, which categorizes each paper to account for
  experimental-variant errors. Papers are identified by DOI only; no paper text is
  redistributed.

Figure 1, Table 1 and Table S1 have no computational step. Statistics quoted only in the
running text are covered by `tutorial.ipynb` and
`post-processing/semantic_inconsistency.py`.
