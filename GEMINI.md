# GEMINI.md

This file provides context for the Gemini AI agent to effectively assist with the comparative genomics analysis project in this repository.

## Directory Overview

This is a comprehensive bioinformatics research project focused on comparative genomics, with a specific emphasis on *Aspiorhynchus laticeps*. The project is structured as a series of Jupyter notebooks that guide the user through a complete analysis pipeline, from environment setup to advanced phylogenetic and synteny analysis.

The entire project is designed to be run within a Docker container, ensuring a reproducible and consistent environment with all necessary bioinformatics tools and libraries pre-configured.

## Key Files

*   **`README.md`**: The primary documentation for the project. It provides a detailed overview of the analysis pipeline, setup instructions, and a breakdown of the project structure.
*   **`docker-compose.yml`**: Defines the Docker container service, including resource limits and volume mounts.
*   **`.build/Dockerfile`**: Specifies the base Docker image for the analysis environment.
*   **`init.sh`**: The main script for initializing the project environment. It sets up the Docker container for either CPU or GPU-based analysis.
*   **`Notebooks/`**: This directory contains the core of the project—a series of Quarto (`.qmd`) notebooks that constitute the analysis pipeline. They are numbered in the order they should be executed.
*   **`env/`**: Contains Conda environment files (`.yml`) that define the project's dependencies. `env_biopython.yml` is the primary environment for Python-based analysis.
*   **`CLAUDE.md` / `CRUSH.md`**: These files contain specific instructions for AI agents working within this repository, outlining common commands and coding conventions.

## Usage

The project is intended to be used by running the Jupyter notebooks in the `Notebooks/` directory in sequential order.

### Building and Running

1.  **Initialize the Environment**:
    *   For CPU-only analysis:
        ```bash
        ./init.sh
        ```
    *   For GPU-accelerated analysis:
        ```bash
        ./init.sh GPU
        ```

2.  **Access Jupyter Lab**:
    *   Once the container is running, access the Jupyter Lab interface at `http://localhost:38888`. A token may be required, which will be printed in the console upon startup.

3.  **Execute Notebooks**:
    *   Follow the notebooks in numerical order, starting with `0-ultimate-env_construction.qmd` to ensure all environments and tools are correctly configured.

### Development Conventions

*   **Environment Management**: The project uses Conda for managing environments. The primary environment for Python analysis is `env_biopython`.
*   **Code Style**:
    *   **Python**: Adheres to PEP 8, with `snake_case` for variables and functions and `PascalCase` for classes. Type hints are encouraged.
    *   **Shell Scripts**: Should be robust, using `set -euo pipefail` where appropriate.
*   **AI Agent Interaction**: The `CLAUDE.md` and `CRUSH.md` files provide specific instructions for AI agents, including commands for environment activation and running key tools like OrthoFinder.
