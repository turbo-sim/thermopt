"""Reproduce the updated study, figures, validation and all report formats."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import argparse
from paths import LOGS

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-only", action="store_true", help="Rebuild the reports from existing saved tables and figures.")
    args = parser.parse_args()
    LOGS.mkdir(parents=True, exist_ok=True)
    stages = ("build_report.py",) if args.report_only else (
        "study.py", "refine_results.py", "run_sensitivity_efficiency_and_exploitation.py",
        "run_fluid_comparison.py", "validate_saved_results.py", "compare_results.py",
        "build_report.py")
    for script in stages:
        log = LOGS / (Path(script).stem + ".log")
        print(f"Running {script}; log: {log}", flush=True)
        with log.open("w", encoding="utf-8") as stream:
            subprocess.run([sys.executable, str(HERE / script)], stdout=stream, stderr=subprocess.STDOUT, check=True)
    python = Path(sys.executable)
    if importlib.util.find_spec("reportlab") is None:
        bundled = Path.home() / ".cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
        if not bundled.exists():
            raise RuntimeError("PDF generation requires reportlab and Pillow; install them or run render_report_pdf.py in a runtime that supplies them.")
        python = bundled
    with (LOGS / "render_report_pdf.log").open("w", encoding="utf-8") as stream:
        subprocess.run([str(python), str(HERE / "render_report_pdf.py")], stdout=stream, stderr=subprocess.STDOUT, check=True)
    print("Finished. Share report_partial_evaporation.html as a single file.", flush=True)


if __name__ == "__main__":
    main()
