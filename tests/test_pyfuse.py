import os
import pytest
import logging
import re
import pandas as pd
from datetime import datetime
from os import environ
from pyfuse.annotator import RunPyfuse
from pyfuse.utils.common_utils import config, utils
from pyfuse.utils.resource_manger import ResourceManager
from pandas._testing import assert_frame_equal
pd.options.display.width = 0

logger = logging.getLogger()

DEFAULT_PYTEST_INPUT = "input_27cases.txt"
DEFAULT_PYTEST_TRUTH_OUTPUT = "expected_output_27cases.xlsx"
DEFAULT_PYTEST_INPUT_FORMAT = "star"


"""
This module tests the PyFuse annotation output against expected output from perl version of PyFuse


To run the tests, use the following command:

pytest -vrs --html=pyfuse_pytest_report_Mar22-2026.html --self-contained-html -vrs

or use the following command with optional environment variables:

RESOURCE_PATH="./custom_resources/hg38/" \
INPUT_FORMAT="star"  \
FUSION2VCF=TRUE  \
REFERENCE_GENOME="./custom_resources/hg38/hg38.fa" \
pytest -vrsx --html=pyfuse_pytest_report_Feb26-2024.html --self-contained-html -s

+++++++++++ optional environmentatl variables ++++++++++
TRUTH_OUTPUT : path of perl output file from perl version of pyfuse
INPUT_FORMAT : options: default/star
REFERENCE_GENOME : refence genome file .fa 
OUTPUT_PATH: any output 
RESOURCE_PATH: path to the resources folder
FUSION2VCF: option-TRUE 

"""


def parse_truth_output_file(truth_fus_out_file):

    perl_fus_df = pd.read_csv(truth_fus_out_file, sep='\t', dtype={"Distance_between_breakpoints": object})
    perl_fus_df["5'-3'Gene_Partners"] = perl_fus_df["5'-3'Gene_Partners"].str.replace('~', '-')
    perl_fus_df["Distance_between_breakpoints"] = perl_fus_df['Distance_between_breakpoints'].str.replace('-', 'NA')

    if environ.get('perl') is not None:
        perl_fus_df.rename({
                            "3\'_Frame": 'Frame_3p',
                            "5\'_Frame": 'Frame_5p',
                            'Fusion_Location': 'Fusion_Position',
                            'Fusion_Sequence': 'Fusion-Sequence',
                            }, axis=1, inplace=True)
        #perl_fus_df["5\'co-ordinate"] = perl_fus_df["Chr_5\'"] + ':' + perl_fus_df["Coordinate_5\'"].astype(str)
        #
        # 
        # perl_fus_df["3\'co-ordinate"] = perl_fus_df["Chr_3\'"] + ':' + perl_fus_df["Coordinate_3\'"].astype(str)

    return perl_fus_df.fillna('NA')


def python_out(truth_fus_out_file):

    # Read the Excel file, keep int columns as object to avoid float conversion, and fill NaN with 'NA'
    df = pd.read_excel(truth_fus_out_file, dtype=object)
    df = df.fillna('NA')
    return df

def fix_cols(df, reference_genome=None):

    print(type(df))
    COLS = ["5'-3'Gene_Partners", 'Fusion_Position', 'Frame_Status',
            'Fusion_Annotation', "5'_Exon_Annotation", "3'_Exon_Annotation",
            "5'co-ordinate", "3'co-ordinate", "Frame_5p", "Frame_3p",
            'Distance_between_breakpoints']
    
    if 'Fusion-Sequence' in df.columns and reference_genome is not None:
        df['Fusion_nucleotide_sequence'] = df['Fusion-Sequence']
        COLS.append('Fusion_nucleotide_sequence')

    elif 'Fusion_nucleotide_sequence' in df.columns:
        COLS.append('Fusion_nucleotide_sequence')   

    print('\nColumns in the output:', COLS)

    """
    df = df.apply(pd.to_numeric,
                          downcast='integer',
                          errors='ignore')"""
    #print('\nData types in the output:', df.dtypes)


    #ignoring one of the python cases
    return df[COLS]


def fix_df(df):
    return df.sort_values(by=df.columns.tolist()).reset_index(drop=True)


_INT_LIKE_RE = re.compile(r"^[+-]?\d+(?:\.0+)?$")
_FLOAT_LIKE_RE = re.compile(r"^[+-]?\d+\.\d+$")


def _canon_cell(value):
    if pd.isna(value):
        return 'NA'

    text = str(value).strip()
    if text in {'', 'NA', 'nan', 'None', '<NA>'}:
        return 'NA'

    if _INT_LIKE_RE.match(text):
        try:
            return str(int(float(text)))
        except ValueError:
            return text

    if _FLOAT_LIKE_RE.match(text):
        try:
            return format(float(text), '.15g')
        except ValueError:
            return text

    return text


def _canonicalize_df(df):
    return df.apply(lambda col: col.map(_canon_cell))


def _build_diff_report(expected_df, actual_df, max_examples=5):
    common_cols = [c for c in expected_df.columns if c in actual_df.columns]
    lines = []

    if list(expected_df.columns) != list(actual_df.columns):
        lines.append(f"Column order differs. expected={list(expected_df.columns)} actual={list(actual_df.columns)}")

    if expected_df.shape != actual_df.shape:
        lines.append(f"Shape differs. expected={expected_df.shape} actual={actual_df.shape}")

    expected_norm = _canonicalize_df(expected_df[common_cols])
    actual_norm = _canonicalize_df(actual_df[common_cols])

    for col in common_cols:
        mism_mask = expected_norm[col] != actual_norm[col]
        mismatch_count = int(mism_mask.sum())
        if mismatch_count == 0:
            continue

        lines.append(f"Column '{col}' mismatch count: {mismatch_count}/{len(expected_norm)}")
        mismatch_idx = expected_norm.index[mism_mask].tolist()[:max_examples]
        for idx in mismatch_idx:
            left_raw = expected_df.at[idx, col]
            right_raw = actual_df.at[idx, col]
            left_norm = expected_norm.at[idx, col]
            right_norm = actual_norm.at[idx, col]
            lines.append(
                f"  row={idx} expected={repr(left_raw)} ({type(left_raw).__name__}) -> {repr(left_norm)}; "
                f"actual={repr(right_raw)} ({type(right_raw).__name__}) -> {repr(right_norm)}"
            )

    if not lines:
        lines.append("No value-level mismatches found after canonicalization.")

    return "\n".join(lines)


def get_fusion_output(test_bkpt, truth_output, format_, reference_genome=None):
    param_list = []
    #python_bkpt, res_path

    pyfuse_obj = RunPyfuse(test_bkpt, pd.DataFrame(),
                            format_, ref=reference_genome)
    # 
    #   def __init__(self, input_bkpt, target_df, format, out_path, ref=None):
    pyfuse_obj.run_pipeline()
    print('fixing PyFuse output columns')
    test_out = fix_cols(pyfuse_obj.annotated_df, reference_genome)
    print('fixing TRUTH output columns')
    truth_output = fix_cols(truth_output, reference_genome)
    
    #python_out.to_csv('scriptout.csv')
    #TRUTH_OUTPUT.to_csv('TRUTH_OUTPUT.csv')
    #python_bkpt_df = pd.read_csv(test_bkpt, sep='\t')
    #print(TRUTH_OUTPUT)

    grp_cols = ["5'co-ordinate", "3'co-ordinate"]#, "5'_Exon_Annotation", "3'_Exon_Annotation"]
    test_gpby = test_out.groupby(grp_cols)
    truth_gpby = truth_output.groupby(grp_cols)
    #print(truth_gpby.describe())

    for c, (bkpt, python_gp) in enumerate(test_gpby):
        try:
            perl_gp = fix_df(truth_gpby.get_group(bkpt))
            #print(bkpt, perl_gp)
        except KeyError as e:
            perl_gp = pd.DataFrame(columns=python_gp.columns)
            #print('Keyerror', e, '\nBKPT:', bkpt, perl_gp)
        python_gp = fix_df(python_gp)

        #print(f'\n\nTest case-{c}: {bkpt}', python_gp["5'-3'Gene_Partners"].unique())
        #print(f'\nperl:\n{perl_gp}')
        #print(f'\npython:\n{python_gp}\n')
        case_id = f"fusion_{c:03d}_{bkpt[0]}_{bkpt[1]}"
        param_list.append((case_id, perl_gp, python_gp))


    '''for index, row in python_bkpt_df.iterrows():
        fus_id = row['Fusion_id']
        param_list.append((get_frame(python_out, fus_id), get_frame(script_out, fus_[perlid)))
        # assert_frame_equal(, , check_dtype=False)#, check_like=True)'''
    return param_list


def _fixture_root() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')


def _resolve_test_input_path() -> str:
    if os.environ.get('INPUT_BKPT'):
        return os.path.abspath(os.environ['INPUT_BKPT'])
    return os.path.join(_fixture_root(), 'input', DEFAULT_PYTEST_INPUT)


def _resolve_truth_output_path() -> str:
    if os.environ.get('TRUTH_OUTPUT'):
        return os.path.abspath(os.environ['TRUTH_OUTPUT'])
    return os.path.join(_fixture_root(), 'truth_output', DEFAULT_PYTEST_TRUTH_OUTPUT)


def _resolve_input_format() -> str:
    return os.environ.get('INPUT_FORMAT', DEFAULT_PYTEST_INPUT_FORMAT)


def _resolve_reference_genome():
    if os.environ.get('REFERENCE_GENOME'):
        return os.path.abspath(os.environ['REFERENCE_GENOME'])
    return None


def _resolve_output_path() -> str:
    if os.environ.get('OUTPUT_PATH'):
        base = os.path.abspath(os.environ['OUTPUT_PATH'])
    else:
        base = os.path.abspath(os.getcwd())

    out_path = os.path.join(base, f"pytest_results_{datetime.now().strftime('%d-%m-%Y_%I-%M-%S')}")
    os.makedirs(out_path, exist_ok=True)
    return out_path


def _configure_runtime_resources() -> None:
    res_manager = ResourceManager(cache_dir=os.environ.get('CACHE_DIR'))
    cfg_assets = res_manager.resolve_config_assets()
    config.update({k: str(v) for k, v in cfg_assets.items()})

    if os.environ.get('RESOURCE_PATH'):
        bundle = res_manager.resolve(user_dir=os.path.abspath(os.environ['RESOURCE_PATH']))
    else:
        genome = os.environ.get('GENOME', 'grch37').lower()
        default_bundle = f"default_{genome}"
        logger.info("Using default genome for tests: %s", genome)
        try:
            bundle = res_manager.resolve_from_cache(default_bundle=default_bundle)
        except FileNotFoundError as exc:
            pytest.skip(
                "No installed default resources found for integration tests. "
                "Set RESOURCE_PATH or run pyfuse resources install first. "
                f"Details: {exc}",
                allow_module_level=True,
            )

    config.update({k: str(v) for k, v in bundle.resource_files.items()})
    utils.res_path = str(bundle.root_dir)
    logger.info("Using resource bundle for tests: %s", bundle.root_dir)


def _build_param_cases():
    test_input_file = _resolve_test_input_path()
    truth_fus_out_file = _resolve_truth_output_path()
    format_ = _resolve_input_format()
    reference_genome = _resolve_reference_genome()

    if not os.path.exists(truth_fus_out_file):
        pytest.fail(f"Truth file {truth_fus_out_file} does not exist.", pytrace=False)

    if os.environ.get('FUSION2VCF') in {'True', 'TRUE'}:
        truth_output = parse_truth_output_file(truth_fus_out_file=truth_fus_out_file)
    else:
        truth_output = python_out(truth_fus_out_file=truth_fus_out_file)

    utils.out_path = _resolve_output_path()
    _configure_runtime_resources()
    return get_fusion_output(test_input_file, truth_output, format_, reference_genome)


_CASES = _build_param_cases()


@pytest.mark.integration
@pytest.mark.parametrize(
    "case_id, EXPECTED, PYFUSE_OUTPUT",
    _CASES,
    ids=[c[0] for c in _CASES],
)
def test_annotation_match(case_id, EXPECTED, PYFUSE_OUTPUT):
    show_pass_detail = os.environ.get("PYFUSE_SHOW_MATCH_DETAILS", "").lower() in {"1", "true", "yes", "y"}

    try:
        assert_frame_equal(EXPECTED, PYFUSE_OUTPUT, check_dtype=False, check_like=True)
        if show_pass_detail:
            msg = f"[PASS] {case_id}: exact dataframe match"
            print(msg)
            logger.info(msg)
    except AssertionError as primary_exc:
        expected_norm = _canonicalize_df(EXPECTED)
        output_norm = _canonicalize_df(PYFUSE_OUTPUT)
        try:
            assert_frame_equal(expected_norm, output_norm, check_dtype=False, check_like=True)
            if show_pass_detail:
                msg = f"[PASS] {case_id}: matched after canonicalization"
                print(msg)
                logger.info(msg)
        except AssertionError:
            detail = _build_diff_report(EXPECTED, PYFUSE_OUTPUT)
            pytest.fail(
                f"{primary_exc}\n\nDetailed mismatch for {case_id}:\n{detail}",
                pytrace=False,
            )
