import json
from typing import List, Literal, Annotated
from collections import Counter
from pydantic import BaseModel, Field, ValidationError, ConfigDict
from tqdm import tqdm

# --- CONFIGURATION ---
INPUT_FILE = "raw_dataset.jsonl" 
CLEANED_RAW_OUT = "cleaned_raw_dataset.jsonl"
FAILED_LOG = "failed_schema_checker.jsonl"

# --- SCHEMA DEFINITION ---
RatioTriple = Annotated[List[str], Field(min_length=3, max_length=3)]

class BaseStrictModel(BaseModel):
    # 'ignore' tells Pydantic to strip extra fields during model_validate()
    model_config = ConfigDict(extra='ignore') 

class TargetSchema(BaseStrictModel):
    name: str
    form: Literal["powder", "thin film", "other"]

class ReagentSchema(BaseStrictModel):
    name: str
    role: str

class PhasePuritySchema(BaseStrictModel):
    impurity_phase: List[str]
    phase_purity_details: str
    classification: Literal["Pure", "Impure", "Insufficient info"]

class OperationSchema(BaseStrictModel):
    operation_type: Literal["mixing", "heating", "deposition", "pressing", "grinding", "other"]
    operation_details: str = ""
    operation_time: str = ""
    temperature_value: str = ""
    pH: str = ""
    ratios: List[RatioTriple]
    atmosphere: str = ""

class RecipeSchema(BaseStrictModel):
    doi: str
    target: TargetSchema
    metal_precursors: List[str]
    reagents: List[ReagentSchema]
    phase_purity: PhasePuritySchema
    characterization_methods: List[str]
    operations: List[OperationSchema]

def _classify_validation_error(e: ValidationError) -> str:
    """Returns the most specific failure category for a ValidationError."""
    for err in e.errors():
        if err["type"] == "missing":
            return "missing_field"
        if err["type"] == "literal_error":
            return "bad_literal"
    return "type_error"


def run_checker():
    stats = {
        "total": 0,
        "perfect": 0,
        "fixed": 0,
        "fail_missing_field": 0,
        "fail_bad_literal": 0,
        "fail_type_error": 0,
        "fail_malformed_json": 0,
        "warn_empty_doi": 0,
        "warn_empty_target_name": 0,
    }

    extra_field_counts = Counter()  
    allowed_keys = set(RecipeSchema.model_fields.keys())

    with open(INPUT_FILE, "r", encoding="utf-8") as f, \
         open(CLEANED_RAW_OUT, "w", encoding="utf-8") as v_f, \
         open(FAILED_LOG, "w", encoding="utf-8") as i_f:

        for line in tqdm(f, desc="Checking & Fixing"):
            if not line.strip():
                continue
            stats["total"] += 1
            obj = None

            try:
                obj = json.loads(line)

                # Record which unexpected top-level keys are present
                extra_found = set(obj.keys()) - allowed_keys
                for key in extra_found:
                    extra_field_counts[key] += 1

                # Validate + strip extra fields
                valid_recipe = RecipeSchema.model_validate(obj)

                # Semantic empty-field warnings — record passes but is flagged
                if not valid_recipe.doi.strip():
                    stats["warn_empty_doi"] += 1
                if not valid_recipe.target.name.strip():
                    stats["warn_empty_target_name"] += 1

                v_f.write(valid_recipe.model_dump_json() + "\n")

                if extra_found:
                    stats["fixed"] += 1
                else:
                    stats["perfect"] += 1

            except json.JSONDecodeError as e:
                # Malformed JSON 
                stats["fail_malformed_json"] += 1
                i_f.write(json.dumps({"doi": None, "error_type": "malformed_json", "error": str(e)}) + "\n")

            except ValidationError as e:
                # Pydantic failure categorised by error type
                category = _classify_validation_error(e)
                stats[f"fail_{category}"] += 1
                i_f.write(json.dumps({
                    "doi": obj.get("doi") if obj else None,
                    "error_type": category,
                    "error": str(e)
                }) + "\n")

    # --- Report ---
    total_fails = stats["fail_missing_field"] + stats["fail_bad_literal"] + \
                  stats["fail_type_error"] + stats["fail_malformed_json"]

    print(f"\n{'='*45}")
    print(f"  VALIDATION REPORT")
    print(f"{'='*45}")
    print(f"  Total processed:        {stats['total']:>7,}")
    print(f"  Perfect (no extras):    {stats['perfect']:>7,}")
    print(f"  Fixed (extras stripped):{stats['fixed']:>7,}")
    print(f"  Critical failures:      {total_fails:>7,}")
    print(f"    - Missing field:       {stats['fail_missing_field']:>6,}")
    print(f"    - Bad literal value:   {stats['fail_bad_literal']:>6,}")
    print(f"    - Type error:          {stats['fail_type_error']:>6,}")
    print(f"    - Malformed JSON:      {stats['fail_malformed_json']:>6,}")
    print(f"\n  Semantic warnings (passed validation):")
    print(f"    - Empty DOI:           {stats['warn_empty_doi']:>6,}")
    print(f"    - Empty target name:   {stats['warn_empty_target_name']:>6,}")

    if extra_field_counts:
        print(f"\n  Stripped extra fields:")
        for key, count in extra_field_counts.most_common():
            print(f"    - '{key}': {count:,}")
    print(f"{'='*45}\n")


if __name__ == "__main__":
    run_checker()