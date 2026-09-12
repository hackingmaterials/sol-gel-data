import json
from collections import Counter
from tqdm import tqdm

# --- CONFIGURATION ---
# Runs on the normalized dataset, i.e. the output of normalization.py. The
# checks below compare parsed chemistry (material_formula / composition), not
# raw strings, so they cannot run inside schema_checker.py — that script sees
# the dataset before text2chem parsing.
INPUT_FILE = "sol_gel_dataset.jsonl"
FLAGGED_LOG = "semantic_inconsistencies.jsonl"

# C, H, N and O are never treated as unsourced: they are the elements of
# carbonates (BaCO3, SrCO3 ...) and of the organic residues left by the
# chelating agents and alkoxides used in sol-gel synthesis, none of which need
# a source in the target or in the metal precursors to be a plausible impurity.
IGNORED_ELEMENTS = {"C", "H", "N", "O"}


def material_elements(material):
    """Every element appearing in a parsed material entity's composition."""
    elements = set()
    for compound in material.get("composition", []):
        elements.update(compound.get("elements", {}).keys())
    return elements


def check_record(record):
    """
    Returns the list of semantic inconsistencies found in one normalized
    record. An empty list means the record is self-consistent.

    None of these are schema violations — every flagged record is valid
    against the schema and keeps its other fields (precursors, operations,
    conditions), which is why they are reported rather than removed.

    Checks 2 and 3 compare two fields of the same record. Checks 1 and 4 depend
    on the text2chem parse and flag parsing artifacts and legitimate chemistry
    alongside genuine inconsistencies; see the comments on each.
    """
    issues = []

    protocol = record.get("protocol", {})
    phase_purity = record.get("phase_purity", {})
    target = protocol.get("target", {})
    impurities = phase_purity.get("impurity_phase", [])
    classification = phase_purity.get("classification", "")

    # 1. The target is listed as one of its own impurity phases. Compared on the
    #    parsed formula, so a phase descriptor cannot hide the match: target
    #    'TiO2' and impurity 'rutile TiO2' both parse to 'TiO2'.
    #    The comparison also matches records that are correct as written:
    #    composite targets reduced to one component ('TiO2:SiO2:Fe3O4' against
    #    impurity 'Fe3O4'), and polymorphs of the target ('Al2O3' against
    #    'θ-Al2O3'), which are real impurity phases. A verbatim repeat
    #    ('La2Zr2O7' -> ['La2Zr2O7']) is the unambiguous case.
    target_formula = target.get("material_formula", "").strip()
    if target_formula and any(
        imp.get("material_formula", "").strip() == target_formula for imp in impurities
    ):
        issues.append("target_listed_as_impurity")

    # 2/3. The classification and the impurity list contradict each other.
    if classification == "Pure" and impurities:
        issues.append("pure_with_impurity_phase")
    if classification == "Impure" and not impurities:
        issues.append("impure_without_impurity_phase")

    # 4. An impurity phase contains an element that is in neither the target nor
    #    any metal precursor, so the record shows no source for it.
    #    An incomplete parse on either side produces the same flag: the element
    #    is present in the raw string but was not parsed out (composite target
    #    'Ni:ZrO2:Sm2O3', precursor 'zirconium propoxide'), the string yields a
    #    placeholder pseudo-element ('(Ba,Sr)2TiSi2O8' -> 'M'), or the string
    #    did not parse at all ('diatomite', which is SiO2). Where the parse is
    #    complete, the element may still come from a source these fields do not
    #    record: an NaOH mineralizer, a quartz substrate, an alumina crucible,
    #    zirconia milling media.
    if impurities:
        available = material_elements(target)
        for precursor in protocol.get("metal_precursors", []):
            available |= material_elements(precursor)

        unsourced = set()
        for imp in impurities:
            unsourced |= material_elements(imp) - IGNORED_ELEMENTS - available
        if unsourced:
            issues.append("unsourced_element_in_impurity:" + ",".join(sorted(unsourced)))

    return issues


def run_checker():
    stats = Counter()
    total = 0

    with open(INPUT_FILE, "r", encoding="utf-8") as f, \
         open(FLAGGED_LOG, "w", encoding="utf-8") as out_f:

        for line in tqdm(f, desc="Checking consistency"):
            if not line.strip():
                continue
            total += 1

            record = json.loads(line)
            issues = check_record(record)
            if not issues:
                continue

            stats["flagged"] += 1
            for issue in issues:
                stats[issue.split(":")[0]] += 1

            protocol = record.get("protocol", {})
            phase_purity = record.get("phase_purity", {})
            out_f.write(json.dumps({
                "doi": protocol.get("doi", ""),
                "target": protocol.get("target", {}).get("material_string", ""),
                "classification": phase_purity.get("classification", ""),
                "impurity_phase": [
                    imp.get("material_string", "") for imp in phase_purity.get("impurity_phase", [])
                ],
                "issues": issues,
            }) + "\n")

    # --- Report ---
    def pct(n):
        return f"{(100 * n / total):.2f}%" if total else "0.00%"

    print(f"\n{'='*62}")
    print(f"  SEMANTIC CONSISTENCY REPORT")
    print(f"{'='*62}")
    print(f"  Total records:                          {total:>7,}")
    print(f"  Records with >= 1 inconsistency:        {stats['flagged']:>7,}  ({pct(stats['flagged'])})")
    print(f"\n  By check (a record may hit more than one):")
    print(f"    - Target listed as its own impurity:  {stats['target_listed_as_impurity']:>7,}  ({pct(stats['target_listed_as_impurity'])})")
    print(f"    - Labeled pure, impurity phase listed:{stats['pure_with_impurity_phase']:>7,}  ({pct(stats['pure_with_impurity_phase'])})")
    print(f"    - Labeled impure, no impurity listed: {stats['impure_without_impurity_phase']:>7,}  ({pct(stats['impure_without_impurity_phase'])})")
    print(f"    - Impurity has unsourced element:     {stats['unsourced_element_in_impurity']:>7,}  ({pct(stats['unsourced_element_in_impurity'])})")
    print(f"\n  Flagged records are NOT removed from {INPUT_FILE}; they are")
    print(f"  listed in {FLAGGED_LOG} so they can be excluded if needed.")
    print(f"{'='*62}\n")


if __name__ == "__main__":
    run_checker()
