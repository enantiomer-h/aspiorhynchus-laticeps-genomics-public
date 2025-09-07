# CRUSH.md

This file provides guidelines for agentic coding agents working in this repository.

## Commands

### Environment Activation
This project uses Conda for environment management. To activate the primary environment, run:
```bash
conda activate env_biopython
```
Other environments exist in the `env/` directory and can be activated similarly.

### Running OrthoFinder
OrthoFinder is a key tool in this project. To run it on a directory of proteomes:
```bash
conda run -n env_orthofinder orthofinder -f /path/to/proteomes/
```

### Running Jupyter Notebook
To start the Jupyter Notebook server:
```bash
./scripts/start-biopython-kernel.sh
```

## Code Style

### Python
- **Imports**: Group imports as follows: standard library, third-party libraries, and then local application imports. Sort them alphabetically.
- **Formatting**: Adhere to PEP 8. Use an autoformatter like Black if possible.
- **Types**: Use type hints for function signatures.
- **Naming Conventions**: Use `snake_case` for variables and functions. Use `PascalCase` for classes.
- **Error Handling**: Use `try...except` blocks for code that might raise exceptions, such as file I/O or network requests.

### Shell Scripts
- **Formatting**: Use a consistent style. Tools like `shellcheck` and `shfmt` are recommended.
- **Error Handling**: Use `set -euo pipefail` to make scripts more robust.
