# vRAPID Workflow Environments



## 🧪 Overview

The `envs` directory includes YAML files that define the Conda environments required for various stages of the `vRAPID`pipeline. These environments ensure consistent and reproducible execution of bioinformatics tools across different systems.

------

## 📁 Contents

- `envs/env.yml`: The primary environment file containing dependencies for general pipeline operations.
- `envs/vadr.yml`: An environment tailored for genome annotation tasks using VADR.

------

## 🛠️ Usage

To create a Conda environment from one of these YAML files, use the following command:

```
conda env create -f path/to/env.yml
```

Replace `path/to/env.yml` with the appropriate file path. After installation, activate the environment using:

```
conda activate environment_name
```

Ensure that all dependencies are correctly installed before running the pipeline.

------

## 📄 License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/BakelLab/vRAPID/blob/snakefile/LICENSE) file for more details.

------

For detailed information on the pipeline's components and usage, refer to the [vRAPID Wiki](https://github.com/BakelLab/vRAPID/wiki).