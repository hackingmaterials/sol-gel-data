# Sol-Gel Synthesis Dataset

A text-mined dataset of inorganic sol-gel synthesis recipes with phase purity outcomes, extracted from the scientific literature using Gemini 3.0 Flash.

`sol_gel_dataset.jsonl` is the final post-processed dataset. See `tutorial.ipynb` for a walkthrough of the data structure and example analyses.

## Tutorial

To run `tutorial.ipynb`, unzip the dataset and install its dependencies:

```bash
unzip sol_gel_dataset.jsonl.zip
pip install -r requirements.txt
```


## Raw Data

The raw Gemini 3.0 Flash data extraction before any post-processing is provided in `post-processing/raw_dataset.jsonl`. If you want to apply your own normalization scheme (different chemical parser, atmosphere labels, characterization method mappings, etc.), you can re-run the pipeline from this file. 

## Post-Processing Pipeline

### 1. Set up environment

Install dependencies and unzip the raw dataset (the pipeline runs from inside `post-processing`):

```bash
pip install -r post-processing/requirements.txt
cd post-processing
unzip raw_dataset.jsonl.zip
```

### 2. Schema validation

Validates each raw extraction against a Pydantic schema, strips unexpected fields, and logs failures by error type. 

```bash
python schema_checker.py
# input:  raw_dataset.jsonl
# output: cleaned_raw_dataset.jsonl
```

### 3. Normalization

Restructures the dataset and parses chemical entities  using text2chem and pubchempy (`reagent_mapping.json`). Depends on the custom entity mappings in `atmosphere_manual.json`,`cmt_manual.json`. Edit these files to customize atmosphere labels, characterization method names, or reagent lookups before running.

```bash
python normalization.py
# input:  cleaned_raw_dataset.jsonl
# output: sol_gel_dataset.jsonl
```
