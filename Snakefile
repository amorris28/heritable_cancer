# Snakefile
import yaml

# Optional: default configfile
# configfile: "configs/local.yaml"

# Check that config was provided and has required keys
required_keys = ["phs_file", "bim_template", "build", "output", "docker_image", "docker_platform"]

missing_keys = [k for k in required_keys if k not in config]
if missing_keys:
    raise ValueError(
        f"No config file found or missing keys: {missing_keys}\n"
        "Please run with --configfile <path_to_config.yaml>"
    )
rule all:
    input:
        config["output"]

rule harmonize_phs:
    output:
        config["output"]
    shell:
        """
        docker run --platform {config[docker_platform]} --rm -it \
          -v "{workflow.basedir}:{workflow.basedir}" -w "{workflow.basedir}" \
          {config[docker_image]} python3 \
          src/01_preprocess/PHS_PC.py \
	  --phs_file {config[phs_file]} \
          --bim_template {config[bim_template]} \
          --build {config[build]} \
          --out {output}
        """

