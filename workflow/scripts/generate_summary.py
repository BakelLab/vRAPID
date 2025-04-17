#!/usr/bin/env python3

import subprocess
import json
import os

def safe_run(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, text=True).strip()
    except subprocess.CalledProcessError:
        return f"[ERROR] Failed to run: {cmd}"

def main(output_file):
    report_lines = []

    # 1. Conda Environments
    #report_lines.append("=== Conda Environments ===")
    #conda_env_dir = ".snakemake/conda"
    #if os.path.isdir(conda_env_dir):
    #    for env_dir in os.listdir(conda_env_dir):
    #        env_path = os.path.join(conda_env_dir, env_dir)
    #        if os.path.isdir(env_path):
    #            report_lines.append(f"\n--- Environment: {env_dir} ---")
    #            conda_list_cmd = f"{env_path}/bin/conda list"
    #            report_lines.append(safe_run(conda_list_cmd))
    #else:
    #    report_lines.append("No Conda environments found.")

    # 2. Containers
    report_lines.append("\n=== Containers ===")
    container_map_path = ".snakemake/container-mapping.json"
    if os.path.exists(container_map_path):
        with open(container_map_path) as f:
            container_map = json.load(f)
            for container in sorted(set(container_map.values())):
                report_lines.append(f"- {container}")
    else:
        report_lines.append("No container mapping found. Containers may not have been used.")

    # 3. Tool Versions
    report_lines.append("\n=== Tool Versions ===")
    tools = {
        "cutadapt": "cutadapt --version",
        "fastqc": "fastqc --version",
        "kraken2": "krarken2 --version",
        "minimap2": "minimap2 --version",
        "multiqc": "multiqc --version",
        "pilon": "pilon --version",
        "prokka": "prokka --version",
        "qualimap": "qualimap --version | head -n 4 | tail -1",
        "samtools": "samtools --version | head -n 1",
        "shovill": "shovill --version",
    }

    for tool, cmd in tools.items():
        version = safe_run(cmd)
        report_lines.append(f"{tool}: {version}")

    # 4. Workflow Version
    version_file = "workflow/VERSION"
    if os.path.exists(version_file):
        with open(version_file) as f:
            report_lines.append(f"\nWorkflow version: {f.read().strip()}")

    # Save report
    with open(output_file, "w") as f:
        f.write("\n".join(report_lines))

if __name__ == "__main__":
    main(snakemake.output.txt)

