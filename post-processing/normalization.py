import json
import math
import functools
import signal
import re
from tqdm import tqdm

# --- CONFIGURATION ---
INPUT_FILE = "cleaned_raw_dataset.jsonl"
OUTPUT_FILE = "sol_gel_dataset.jsonl"
CMT_MANUAL_FILE = "cmt_manual.json"
ATMO_MANUAL_FILE = "atmosphere_manual.json"
REAGENT_MAPPING_FILE = "reagent_mapping.json"

# --- text2chem Setup ---
try:
    from text2chem.regex_parser import RegExParser
    from text2chem.parser_pipeline import ParserPipelineBuilder
    from text2chem.preprocessing_tools.additives_processing import AdditivesProcessing
    from text2chem.preprocessing_tools.chemical_name_processing import ChemicalNameProcessing
    from text2chem.preprocessing_tools.phase_processing import PhaseProcessing
    from text2chem.preprocessing_tools.mixture_processing import MixtureProcessing
    from text2chem.postprocessing_tools.substitute_additives import SubstituteAdditives

    mp = ParserPipelineBuilder() \
        .add_preprocessing(AdditivesProcessing) \
        .add_preprocessing(ChemicalNameProcessing) \
        .add_preprocessing(PhaseProcessing) \
        .add_preprocessing(MixtureProcessing)\
        .add_postprocessing(SubstituteAdditives)\
        .set_regex_parser(RegExParser)\
        .build()
except ImportError:
    print("text2chem not installed. Chemistry parsing disabled.")
    mp = None

try:
    from pymatgen.core import Composition
    PYMATGEN_AVAILABLE = True
except ImportError:
    print("pymatgen not installed. Oxidation-state prediction disabled.")
    PYMATGEN_AVAILABLE = False

# --- MAPPINGS ---
# Fold every dash variant (hyphen, non-breaking hyphen, figure dash, en dash,
# em dash, horizontal bar, minus sign) to a plain ASCII hyphen, and collapse
# any run of whitespace (incl. non-breaking spaces) to a single space.
_DASH_RE = re.compile('[‐‑‒–—―−]')
_WS_RE = re.compile(r'\s+')

def normalize_method_string(s):
    """Canonical form for matching characterization strings against the
    manual map: lowercase, dash variants folded to '-', whitespace collapsed."""
    s = s.lower().strip()
    s = _DASH_RE.sub('-', s)
    s = _WS_RE.sub(' ', s).strip()
    return s

def load_manual_map(manual_filepath, label="normalization"):
    """
    Builds a flat raw-string -> canonical lookup from a curated manual mapping
    file (e.g. cmt_manual.json, atmosphere_manual.json). Keys are normalized
    with normalize_method_string; keys starting with '_' (e.g. '_comment') are
    metadata and skipped.

    Strings not present in the map are left unchanged / dropped downstream.
    """
    lookup = {}
    try:
        with open(manual_filepath, 'r', encoding='utf-8') as f:
            manual = json.load(f)
        for key, val in manual.items():
            if not key.startswith('_'):
                lookup[normalize_method_string(key)] = val
    except FileNotFoundError:
        print(f"Warning: {manual_filepath} not found. {label} disabled.")

    return lookup

CHAR_LOOKUP = load_manual_map(CMT_MANUAL_FILE, label="Characterization normalization")
ATMO_LOOKUP = load_manual_map(ATMO_MANUAL_FILE, label="Atmosphere normalization")

# Combinatorial dilute-H2 reducing mixtures (5% h2/ar, 30% h2 + 70% n2,
# ar-4%h2, 1 vol.% h2 in ar, ...) are too varied to enumerate in the manual
# map, so classify them by structure: H2 present + an inert carrier present +
# a mixture marker (%, /, +, :, "in", "and").
_FG_H2_RE = re.compile(r'\bh2\b|hydrogen')
_FG_INERT_RE = re.compile(r'\bar\b|\bn2\b|\bhe\b|argon|nitrogen|helium|forming')
_FG_MIX_RE = re.compile(r'[%/+:]|\bin\b|\band\b')

def classify_atmosphere(raw):
    """Normalize a raw atmosphere string to a canonical label, or "" if
    unrecognized (the raw string is preserved separately). Exact manual-map
    matches win; otherwise the dilute-H2 mixture regex fallback applies."""
    if not raw:
        return ""
    s = normalize_method_string(raw)
    if s in ATMO_LOOKUP:
        return ATMO_LOOKUP[s]
    if _FG_H2_RE.search(s) and _FG_INERT_RE.search(s) and _FG_MIX_RE.search(s):
        return "Forming Gas (dilute H2 in N2/Ar)"
    return ""

# --- UTILITIES ---
def parse_numeric_range(raw_str):
    # Text fields (raw_string, unit) use "" when empty; numeric fields
    # (max_value, min_value) use null, since JSON has no NaN and a numeric
    # column should not be polluted with empty strings.
    if not raw_str or not isinstance(raw_str, str):
        return {"raw_string": raw_str if isinstance(raw_str, str) else "", "unit": "", "max_value": None, "min_value": None}
    clean_s = raw_str.lower().strip()
    nums = [float(n) for n in re.findall(r"[-+]?\d*\.\d+|\d+", clean_s)]
    # A unit is only meaningful alongside a magnitude. If no number is parsed
    # (e.g. "room temperature", "overnight", "three days"), report nothing and
    # leave the raw string for downstream handling.
    if not nums:
        return {"raw_string": raw_str, "unit": "", "max_value": None, "min_value": None}
    # The unit must be the token immediately following a digit, not just any
    # trailing word, so "80 °C then cooled" -> "°c" (not "cooled").
    unit_match = re.search(r'\d\s*([a-zA-Z°%]+)', clean_s)
    unit = unit_match.group(1) if unit_match else ""
    if unit:
        unit_map = {"hour": "h", "hours": "h", "hr": "h", "hrs": "h", "sec": "s", "second": "s", "minute": "min", "minutes": "min",
                    "°c": "°C", "°f": "°F", "k": "K"}
        unit = unit_map.get(unit, unit)
    min_val, max_val = min(nums), max(nums)
    return {"raw_string": raw_str, "unit": unit, "max_value": max_val, "min_value": min_val}

def clean_reagent_name(name):
    n = name.lower().strip()
    # Greedily catch all water variants ("deionized water, distilled water, DI water, etc.")
    if "water" in n:
        return "water"
    overrides = {"triethanol amine": "triethanolamine", "absolute ethanol": "ethanol"} # add variants that pubchempy (reagent_mappings.json) fails to capture
    return overrides.get(n, n)

HYDRATE_SEP_RE = re.compile(r'[.\*•×∙⋅]\s*(?=\d*H2O)')

def canonicalize_hydrate_sep(formula):
    """
    Standardizes hydrate separators to the canonical middle-dot '·'.
    Recognizes period, asterisk, bullet (U+2022), multiplication sign
    (U+00D7), bullet operator (U+2219), and dot operator (U+22C5),
    but only when followed by '\\d*H2O' so non-hydrate dots stay intact.

    Example: 'La(NO3)3.6H2O'  -> 'La(NO3)3·6H2O'
             'Ni(NO3)2*6H2O'  -> 'Ni(NO3)2·6H2O'
             'Li1.5La0.5TiO3' -> 'Li1.5La0.5TiO3' (no H2O, preserved)
    """
    if not formula or not isinstance(formula, str):
        return formula
    return HYDRATE_SEP_RE.sub('·', formula).strip()

OXI_TIMEOUT_SECONDS = 30

class _OxiTimeout(Exception):
    pass

def _oxi_timeout_handler(signum, frame):
    raise _OxiTimeout()

@functools.lru_cache(maxsize=4096)
def predict_oxidation_states(formula):
    """
    Returns {element: [ox_state_int, ...]} using pymatgen's best guess.
    Fractional guesses (e.g. Fe in Fe3O4 -> 2.667) collapse to [floor, ceil].

    Memorized — same precursor formulas recur across many recipes, so
    pymatgen only runs once per distinct formula (and once per distinct
    failure too, since failures return {} rather than raising).

    Hydrates: predicts on the anhydrous part only (text before '·'), since
    pymatgen.Composition cannot parse the middle-dot separator.

    Doped / fractional-stoichiometry formulas (e.g. La0.7Sr0.3MnO3) are
    attempted too — the wall-clock timeout below guards against the slow
    combinatorial searches these can trigger.

    Returns {} on any failure (no pymatgen, unparseable formula, no balanced
    assignment, formula contains variables, timeout, etc.).
    """
    if not PYMATGEN_AVAILABLE or not formula or not isinstance(formula, str):
        return {}

    # Hydrates: predict on the anhydrous salt only.
    formula = formula.split('·', 1)[0].strip()
    if not formula:
        return {}

    try:
        # Hard wall-clock timeout — pymatgen's oxi_state_guesses can stall on
        # pathological compositions (organic precursors, rare-element systems,
        # mixed-valence enumerations). Bail out and cache {} via lru_cache so
        # the same formula is never paid for twice.
        signal.signal(signal.SIGALRM, _oxi_timeout_handler)
        signal.alarm(OXI_TIMEOUT_SECONDS)
        try:
            guesses = Composition(formula).oxi_state_guesses()
        finally:
            signal.alarm(0)
        if not guesses:
            return {}
        best = guesses[0]
        result = {}
        for el, state in best.items():
            el_str = str(el)
            if float(state).is_integer():
                result[el_str] = [int(state)]
            else:
                result[el_str] = [math.floor(state), math.ceil(state)]
        return result
    except Exception:
        return {}

def _empty_chem(name, form=""):
    """
    Canonical 'failed/empty' material entity. Returns the SAME key set as a
    successful parse so every material entity has an identical schema
    regardless of outcome — no downstream consumer has to distinguish
    "key absent" from "key present but empty". An empty `composition` ([])
    is the single, unambiguous "not normalized" signal.
    """
    return {
        "material_string": name,
        "material_form": form,
        "material_formula": "",
        "composition": [],
        "additives": [],
        "is_mixture": False,
        "phase": "",
        "oxygen_deficiency": "no",
        "oxidation_states": {},
    }

def parse_chem(target_obj):
    """
    Comprehensive chemical parser:
    - Captures Phase, Formula, Species, and Stoichiometry.
    - Normalizes hydrate dots in material_formula.
    """
    if not mp or not isinstance(target_obj, dict):
        name = target_obj.get('name', '') if isinstance(target_obj, dict) else str(target_obj)
        name = canonicalize_hydrate_sep(name)
        return _empty_chem(name)

    # Canonicalize hydrate separators (*, •, ×, etc.) to '·' before
    # text2chem sees the string — text2chem only recognizes a fixed set
    # of separators, and a non-canonical one causes silent mis-splits.
    name = canonicalize_hydrate_sep(target_obj.get('name', ''))
    target_form = target_obj.get('form', '')
    
    try:
        res_obj = mp.parse(name)
        res_dict = res_obj.to_dict()
        
        raw_od = getattr(res_obj, "oxygen_deficiency", None)
        od_status = "yes" if raw_od else "no"

        phase_str = res_dict.get("phase", "")
        full_composition = []
        all_elements = set()
        
        if hasattr(res_obj, 'composition') and res_obj.composition:
            for compound in res_obj.composition:
                comp_elements = dict(getattr(compound, 'elements', {}))
                comp_species = dict(getattr(compound, 'species', {}))
                
                # Normalize sub-compound formulas if they exist
                raw_comp_formula = getattr(compound, 'formula', "")
                clean_comp_formula = canonicalize_hydrate_sep(raw_comp_formula)
                
                comp_data = {
                    "formula": clean_comp_formula,
                    "amount": getattr(compound, 'amount', "1"),
                    "elements": comp_elements,
                    "species": comp_species
                }
                full_composition.append(comp_data)
                all_elements.update(comp_elements.keys())

        # Extract and normalize the main material formula
        raw_material_formula = getattr(res_obj, "material_formula", "")
        normalized_formula = canonicalize_hydrate_sep(raw_material_formula)

        return {
            "material_string": name,
            "material_form": target_form,
            "material_formula": normalized_formula,
            "composition": full_composition,
            "additives": getattr(res_obj, "additives", []),
            "is_mixture": len(full_composition) > 1,
            "phase": phase_str,
            "oxygen_deficiency": od_status,
            "oxidation_states": predict_oxidation_states(normalized_formula),
        }
    except Exception:
        return _empty_chem(name, target_form)

# --- MAIN ---
def run_normalization():
    print(f"Loading reagent mapping from {REAGENT_MAPPING_FILE}...")
    try:
        with open(REAGENT_MAPPING_FILE, "r", encoding="utf-8") as f:
            reagent_map = json.load(f)
    except FileNotFoundError:
        print(f"Error: {REAGENT_MAPPING_FILE} not found. Please run mapping script first.")
        reagent_map = {}

    with open(INPUT_FILE, "r", encoding="utf-8") as infile, \
         open(OUTPUT_FILE, "w", encoding="utf-8") as outfile:
        
        for line in tqdm(infile, desc="Normalizing Data"):
            recipe = json.loads(line)
            
            recipe['target'] = parse_chem(recipe['target'])
            recipe['metal_precursors'] = [parse_chem({"name": p}) for p in recipe['metal_precursors']]
            
            new_reagents = []
            for r in recipe['reagents']:
                lookup = clean_reagent_name(r['name'])
                
                re_obj = {
                    "reagent_string": r['name'],
                    "reagent_role": r['role'],
                    "pubchem_CID": None,   # numeric -> null when unresolved
                    "iupac_name": "",
                    "reagent_formula": "",
                    "smiles": ""
                }

                if lookup in reagent_map and reagent_map[lookup]:
                    m_data = reagent_map[lookup]
                    re_obj.update({
                        "pubchem_CID": m_data.get('cid'),          # null if missing
                        "iupac_name": m_data.get('iupac_name') or "",
                        "reagent_formula": m_data.get('formula') or "",
                        "smiles": m_data.get('smiles') or ""
                    })
                new_reagents.append(re_obj)
            recipe['reagents'] = new_reagents

            # --- Characterization Methods ---
            raw_methods = recipe.get("characterization_methods", [])
            normalized_list = []
            for m in raw_methods:
                if not m: continue
                clean_m = normalize_method_string(m)
                canon = CHAR_LOOKUP.get(clean_m)
                # Drop strings the manual map doesn't recognize — they are NOT
                # passed through as-is. The unmodified inputs remain available
                # in characterization_string for anyone who needs them.
                if canon is not None:
                    normalized_list.append(canon)
            # Multiple raw strings can collapse to the same canonical method
            # (e.g. "ICP-MS" and "ICP-AES" -> ICP-OES); keep first occurrence.
            normalized_list = list(dict.fromkeys(normalized_list))

            # --- Operations ---
            parsed_ops = []
            for op in recipe["operations"]:
                ph_raw = op["pH"]
                ph_res = parse_numeric_range(ph_raw)
                atmo_raw = op.get("atmosphere", "")
                parsed_ops.append({
                    "operation_type": op["operation_type"],
                    "operation_details": op.get("operation_details", ""),
                    "time": parse_numeric_range(op["operation_time"]),
                    "temperature": parse_numeric_range(op["temperature_value"]),
                    "pH": {"raw_string": ph_raw, "max_value": ph_res["max_value"], "min_value": ph_res["min_value"]},
                    "ratios": [{"type": r[0], "description": r[1], "value": r[2]} for r in op["ratios"]],
                    # Raw string preserved; atmosphere_normalized is "" when unrecognized.
                    "atmosphere": {
                        "raw_string": atmo_raw,
                        "atmosphere_normalized": classify_atmosphere(atmo_raw)
                    }
                })

            # --- Impurity Phases ---
            raw_imp = recipe["phase_purity"]["impurity_phase"]
            parsed_imp = [parse_chem({"name": imp}) for imp in raw_imp] if raw_imp else []

            # --- Assemble restructured output ---
            output = {
                "protocol": {
                    "doi": recipe["doi"],
                    "target": recipe["target"],
                    "metal_precursors": recipe["metal_precursors"],
                    "reagents": recipe["reagents"]
                },
                "operations": parsed_ops,
                "phase_purity": {
                    "impurity_phase": parsed_imp,
                    "phase_purity_details": recipe["phase_purity"].get("phase_purity_details", ""),
                    "classification": recipe["phase_purity"].get("classification", ""),
                    "characterization_methods": {
                        "characterization_string": raw_methods,
                        "characterization_normalized": normalized_list
                    }
                }
            }

            outfile.write(json.dumps(output, ensure_ascii=False) + "\n")

    print(f"\nSaved standardized dataset to: {OUTPUT_FILE}")

if __name__ == "__main__":
    run_normalization()