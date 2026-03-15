#!/usr/bin/env python3
"""
run_pipeline.py
===============
Upload pipeline code to TSD, submit the Nextflow job, wait for completion,
and download the results.

Usage:
    python run_pipeline.py --project p212
    python run_pipeline.py --project p2521
    python run_pipeline.py --project all   # submits both in parallel

Requires TACL to be installed and registered for the project.
Must be run from a whitelisted IP.
"""

import argparse
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# tsd_api lives in the sibling tsd_hpc_api directory
sys.path.insert(0, str(Path(__file__).parent.parent / "tsd_hpc_api"))
import tsd_api

# =============================================================================
# Per-project configuration
# =============================================================================

PROJECTS = {
    "p212": {
        "username": "p212-andrewmo",
        "profile": "p212",
        "pipeline_dir": "/tsd/p212/home/p212-andrewmo/heritable_cancer",
        "import_dir": "/tsd/p212/data/durable/file-import/p212-member-group",
        "export_dir": "/tsd/p212/data/durable/file-export",
        "time_minutes": 120,
        "memory_mb": 16384,
    },
    "p2521": {
        "username": "p2521-andrewmo",
        "profile": "p2521",
        "pipeline_dir": "/tsd/p2521/home/p2521-andrewmo/heritable_cancer",
        "import_dir": "/tsd/p2521/data/durable/file-import/p2521-member-group",
        "export_dir": "/tsd/p2521/data/durable/file-export",
        "time_minutes": 120,
        "memory_mb": 16384,
    },
}

PIPELINE_ZIP = "pipeline.zip"
JOB_SCRIPT_TEMPLATE = Path(__file__).parent / "tsd_run.sh"

# =============================================================================
# Helpers
# =============================================================================

def build_zip() -> None:
    """Create pipeline.zip containing main.nf, nextflow.config, and src/."""
    Path(PIPELINE_ZIP).unlink(missing_ok=True)
    cmd = ["zip", "-r", PIPELINE_ZIP, "main.nf", "nextflow.config", "src/"]
    print(f"Building {PIPELINE_ZIP}...")
    subprocess.run(cmd, check=True)
    print(f"Created {PIPELINE_ZIP}")


def make_job_script(project: str, cfg: dict) -> str:
    """
    Prepend project-specific variable assignments to the generic tsd_run.sh
    and write to a temp file. Returns the temp file path.
    """
    header = (
        f"#!/bin/bash\n"
        f"PROJECT={project}\n"
        f"PROFILE={cfg['profile']}\n"
        f"PIPELINE_DIR={cfg['pipeline_dir']}\n"
        f"IMPORT_DIR={cfg['import_dir']}\n"
        f"EXPORT_DIR={cfg['export_dir']}\n"
    )
    # Strip the shebang from the template so we don't duplicate it
    template = JOB_SCRIPT_TEMPLATE.read_text()
    template_body = "\n".join(
        line for line in template.splitlines()
        if not line.startswith("#!/")
    )
    script_content = header + "\n" + template_body + "\n"

    tmp = tempfile.NamedTemporaryFile(
        mode="w",
        prefix=f"tsd_run_{project}_",
        suffix=".sh",
        delete=False,
        dir=".",
    )
    tmp.write(script_content)
    tmp.close()
    return tmp.name


def run_project(project: str) -> str:
    """Run the full upload-submit-wait-download pipeline for one project."""
    cfg = PROJECTS[project]
    job_script_local = make_job_script(project, cfg)
    job_script_name = Path(job_script_local).name
    result_file = f"public_{project}.zip"

    try:
        print(f"[{project}] Uploading pipeline files...")
        tsd_api.upload_file(project, PIPELINE_ZIP)
        tsd_api.upload_file(project, job_script_local)

        remote_script = f"{cfg['import_dir']}/{job_script_name}"
        print(f"[{project}] Submitting job: {remote_script}")
        job_id = tsd_api.submit_job(
            project,
            cfg["username"],
            remote_script,
            job_name="nextflow",
            nodes=1,
            tasks=1,
            memory_mb=cfg["memory_mb"],
            time_minutes=cfg["time_minutes"],
            working_dir=cfg["pipeline_dir"],
        )

        print(f"[{project}] Waiting for job {job_id}...")
        tsd_api.wait_for_job(
            project, cfg["username"], job_id, timeout_minutes=cfg["time_minutes"]
        )

        print(f"[{project}] Downloading {result_file}...")
        tsd_api.download_file(project, result_file)
        print(f"[{project}] Done. Results saved to {result_file}")
        return result_file

    finally:
        Path(job_script_local).unlink(missing_ok=True)


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Run heritable_cancer pipeline on TSD.")
    parser.add_argument(
        "--project",
        choices=["p212", "p2521", "all"],
        default="p212",
        help="TSD project to run on (default: p212). Use 'all' to submit both in parallel.",
    )
    args = parser.parse_args()

    build_zip()

    targets = list(PROJECTS.keys()) if args.project == "all" else [args.project]

    if len(targets) == 1:
        run_project(targets[0])
    else:
        print(f"Submitting to {targets} in parallel...")
        with ThreadPoolExecutor(max_workers=len(targets)) as pool:
            futures = {pool.submit(run_project, p): p for p in targets}
            for future in as_completed(futures):
                project = futures[future]
                try:
                    result = future.result()
                    print(f"[{project}] Finished: {result}")
                except Exception as exc:
                    print(f"[{project}] FAILED: {exc}")


if __name__ == "__main__":
    main()
