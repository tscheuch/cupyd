import filecmp
from pathlib import Path


class DifferentFilesException(Exception):
    """Custom exception class to raise when different files are encountered"""


ROOT_DIRECTORY = Path(__file__).resolve().parent
LLANQUIHUE = ROOT_DIRECTORY / "llanquihue"
MODFLOW_WORKSPACE = LLANQUIHUE / "MODFLOW"
SWMM_WORKSPACE = LLANQUIHUE / "SWMM"
RESULTS_DIRECTORY = LLANQUIHUE / "results"
MODFLOW_RESULTS_WORKSPACE = RESULTS_DIRECTORY / "MODFLOW"
SWMM_RESULTS_WORKSPACE = RESULTS_DIRECTORY / "SWMM"

# Compare MODFLOW results
modflow_comparison = filecmp.dircmp(MODFLOW_WORKSPACE, MODFLOW_RESULTS_WORKSPACE)

if modflow_comparison.diff_files:
    raise DifferentFilesException(modflow_comparison.diff_files)
if modflow_comparison.left_only:
    raise DifferentFilesException(
        f"MODFLOW file missing on results: `{modflow_comparison.left_only}`"
    )
if modflow_comparison.right_only:
    raise DifferentFilesException(
        f"MODFLOW file missing on target: `{modflow_comparison.right_only}`"
    )


# Compare SWMM results
swmm_comparison = filecmp.dircmp(SWMM_WORKSPACE, SWMM_RESULTS_WORKSPACE)

if swmm_comparison.diff_files:
    diff_files = swmm_comparison.diff_files
    try:
        diff_files.remove("Llanquihue_base.rpt")
    except ValueError:
        print(
            "Make sure you run the simulation before testing the results. The `Llanquihue_base.rpt` file is identical and at least the analysis should have changed."
        )
        quit()
    with open(SWMM_WORKSPACE / "Llanquihue_base.rpt") as rpt_file:
        rpt_file_text = rpt_file.readlines()
    with open(SWMM_RESULTS_WORKSPACE / "Llanquihue_base.rpt") as result_rpt_file:
        result_rpt_file_text = result_rpt_file.readlines()
    # Skipping final analysis part on both files.
    if rpt_file_text[:-3] != result_rpt_file_text[:-3]:
        raise DifferentFilesException("The file `Llanquihue_base.rpt` are different.")
    if diff_files:
        raise DifferentFilesException(diff_files)
if swmm_comparison.left_only:
    raise DifferentFilesException(f"SWMM file missing on results: `{swmm_comparison.left_only}`")
if swmm_comparison.right_only:
    raise DifferentFilesException(f"SWMM file missing on target: `{swmm_comparison.right_only}`")
