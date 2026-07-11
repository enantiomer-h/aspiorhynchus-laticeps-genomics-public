"""
Path Configuration Loader for Comparative Genomics Notebooks

This module provides functions to load and resolve path configurations
from the central YAML configuration file (paths.yaml).

Usage:
    # Option 1: Import specific functions
    from load_config import load_paths, get_path, verify_paths

    # Load all paths
    paths = load_paths()

    # Access paths via dictionary
    orthogroups_file = paths['orthofinder']['orthogroups_tsv']

    # Option 2: Use convenience function with dot notation
    orthogroups_file = get_path('orthofinder.orthogroups_tsv')

Template Variables:
    The YAML configuration supports template variables in the format ${section.key}
    These are automatically resolved when the configuration is loaded.

    Example in paths.yaml:
        version: "Jan09"
        orthofinder:
          results_dir: "${base.outputs}/OrthoFinder/Results_${version}"

    Resolves to: /home/jovyan/Outputs/OrthoFinder/Results_Jan09

Author: Generated for ComparativeGenomics4AL_Publication
"""

import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

# Check for yaml availability
try:
    import yaml

    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False


# Module-level cache for loaded configuration
_config_cache: Optional[Dict[str, Any]] = None
_config_path: Optional[Path] = None


def _find_config_file() -> Path:
    """
    Find the configuration file in expected locations.

    Search order (paths-set2.yaml is preferred over paths.yaml):
    1. Environment variable GENOMICS_CONFIG_PATH if set
    2. ./config/paths-set2.yaml  (set3 - priority)
    3. ./config/paths.yaml       (set1 - fallback)
    4. ../config/paths-set2.yaml
    5. ../config/paths.yaml
    6. /home/jovyan/Notebooks/config/paths-set2.yaml  (Docker absolute)
    7. /home/jovyan/Notebooks/config/paths.yaml
    8. Relative to this module's location (paths-set2.yaml first, then paths.yaml)

    Returns:
        Path to the configuration file

    Raises:
        FileNotFoundError: If config file cannot be found
    """
    search_paths = []

    # Check environment variable first
    env_path = os.environ.get("GENOMICS_CONFIG_PATH")
    if env_path:
        search_paths.append(Path(env_path))

    # Standard search paths — paths-set2.yaml takes priority over paths.yaml
    search_paths.extend(
        [
            Path("./config/paths-set2.yaml"),
            Path("./config/paths.yaml"),
            Path("../config/paths-set2.yaml"),
            Path("../config/paths.yaml"),
            Path("/home/jovyan/Notebooks/config/paths-set2.yaml"),
            Path("/home/jovyan/Notebooks/config/paths.yaml"),
        ]
    )

    for path in search_paths:
        if path.exists():
            return path.resolve()

    # Try relative to this file's location (paths-set2.yaml first)
    module_dir = Path(__file__).parent
    for filename in ("paths-set2.yaml", "paths.yaml"):
        config_path = module_dir / filename
        if config_path.exists():
            return config_path.resolve()

    raise FileNotFoundError(
        f"Configuration file not found. Searched:\n"
        + "\n".join(f"  - {p}" for p in search_paths)
        + f"\n  - {module_dir / 'paths-set2.yaml'}"
        + f"\n  - {module_dir / 'paths.yaml'}\n"
        "Set GENOMICS_CONFIG_PATH environment variable or place paths-set2.yaml in config directory."
    )


def _get_nested_value(d: Dict, key_path: str) -> Optional[Any]:
    """
    Get a value from a nested dictionary using dot notation.

    Args:
        d: Dictionary to search
        key_path: Dot-separated path (e.g., 'orthofinder.results_dir')

    Returns:
        The value if found, None otherwise
    """
    parts = key_path.split(".")
    current = d

    try:
        for part in parts:
            if isinstance(current, dict):
                current = current[part]
            else:
                return None
        return current
    except (KeyError, TypeError):
        return None


def _resolve_template(
    value: str, config: Dict[str, Any], resolved: Dict[str, Any]
) -> str:
    """
    Resolve template variables in a string value.

    Handles ${section.key} and ${section.nested.key} syntax.
    Also handles simple ${variable} references.

    Args:
        value: String potentially containing template variables
        config: Original configuration dictionary
        resolved: Already resolved values for reference

    Returns:
        String with all template variables resolved
    """
    if not isinstance(value, str):
        return value

    # Pattern matches ${...} with nested dot notation
    pattern = r"\$\{([a-zA-Z_][a-zA-Z0-9_]*(?:\.[a-zA-Z_][a-zA-Z0-9_]*)*)\}"

    max_iterations = 10
    prev_value = None

    for _ in range(max_iterations):
        if prev_value == value:
            break
        prev_value = value

        def replace_var(match: re.Match) -> str:
            var_path = match.group(1)

            # First try resolved values, then original config
            for source in [resolved, config]:
                result = _get_nested_value(source, var_path)
                if result is not None and isinstance(result, str):
                    return result

            # If not found, return original (will be resolved in later pass)
            return match.group(0)

        value = re.sub(pattern, replace_var, value)

    return value


def _update_nested_dict(d: Dict, keys: List[str], value: Any) -> None:
    """Update a nested dictionary with a value at the specified key path."""
    for key in keys[:-1]:
        d = d.setdefault(key, {})
    d[keys[-1]] = value


def _resolve_dict(
    d: Dict[str, Any],
    config: Dict[str, Any],
    resolved: Dict[str, Any],
    path_prefix: str = "",
) -> Dict[str, Any]:
    """
    Recursively resolve all template variables in a dictionary.

    Args:
        d: Dictionary to resolve
        config: Original full configuration
        resolved: Already resolved values
        path_prefix: Current path for nested resolution

    Returns:
        Dictionary with all templates resolved
    """
    result = {}

    for key, value in d.items():
        current_path = f"{path_prefix}.{key}" if path_prefix else key

        if isinstance(value, dict):
            result[key] = _resolve_dict(value, config, resolved, current_path)
        elif isinstance(value, str):
            result[key] = _resolve_template(value, config, resolved)
        elif isinstance(value, list):
            result[key] = [
                _resolve_template(item, config, resolved)
                if isinstance(item, str)
                else item
                for item in value
            ]
        else:
            result[key] = value

        # Update resolved dict for subsequent references
        _update_nested_dict(resolved, current_path.split("."), result[key])

    return result


def load_paths(
    config_path: Optional[Union[str, Path]] = None, force_reload: bool = False
) -> Dict[str, Any]:
    """
    Load and resolve path configuration from YAML file.

    Args:
        config_path: Optional explicit path to config file.
                    If not provided, auto-detects location.
        force_reload: If True, reload even if cached.

    Returns:
        Dictionary with all paths resolved.

    Raises:
        FileNotFoundError: If config file not found
        ImportError: If PyYAML is not installed
        yaml.YAMLError: If config file is invalid YAML

    Example:
        >>> paths = load_paths()
        >>> print(paths['orthofinder']['results_dir'])
        /home/jovyan/Outputs/OrthoFinder/Results_Jan09
    """
    global _config_cache, _config_path

    if not YAML_AVAILABLE:
        raise ImportError(
            "PyYAML is required for configuration loading. "
            "Install with: pip install pyyaml"
        )

    # Use cached version if available
    if _config_cache is not None and not force_reload and config_path is None:
        return _config_cache

    # Find config file
    if config_path:
        path = Path(config_path)
        if not path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
    else:
        path = _find_config_file()

    # Load YAML
    with open(path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    if config is None:
        raise ValueError(f"Configuration file is empty: {path}")

    # Resolve all templates
    resolved: Dict[str, Any] = {}
    result = _resolve_dict(config, config, resolved)

    # Cache the result
    _config_cache = result
    _config_path = path

    return result


def get_path(key_path: str, default: Optional[str] = None) -> Optional[str]:
    """
    Get a specific path by dot-notation key.

    Args:
        key_path: Dot-separated path to the configuration value
                  e.g., "orthofinder.orthogroups_tsv"
        default: Default value if key not found

    Returns:
        The path string, or default if not found

    Example:
        >>> get_path('orthofinder.orthogroups_tsv')
        '/home/jovyan/Outputs/OrthoFinder/Results_Jan09/Orthogroups/Orthogroups.tsv'

        >>> get_path('nonexistent.path', default='/fallback/path')
        '/fallback/path'
    """
    config = load_paths()
    result = _get_nested_value(config, key_path)
    return result if result is not None else default


def get_cds_files() -> Dict[str, str]:
    """
    Get dictionary of all CDS file paths keyed by species name.

    Returns:
        Dictionary mapping species names to CDS file paths

    Example:
        >>> cds_files = get_cds_files()
        >>> print(cds_files['Aspiorhynchus_laticeps'])
        /home/jovyan/Outputs/Preprocessing/BRAKER_UNMASKED/Aspiorhynchus_laticeps/braker.codingseq
    """
    config = load_paths()
    cds_config = config.get("cds", {})

    # Filter out non-species keys (directory references used for path templating)
    exclude_keys = {"base_dir", "unmasked_dir", "masked_dir"}
    return {
        species: path
        for species, path in cds_config.items()
        if species not in exclude_keys and isinstance(path, str)
    }


def get_species_list() -> List[str]:
    """
    Get list of all species from configuration.

    Returns:
        List of species names

    Example:
        >>> species = get_species_list()
        >>> print(len(species))
        13
    """
    config = load_paths()
    return config.get("species", {}).get("all", [])


def get_focal_species() -> str:
    """
    Get the focal species name from configuration.

    Returns:
        Focal species name (default: 'Diptychus_maculatus')
    """
    config = load_paths()
    return config.get("species", {}).get("focal", "Diptychus_maculatus")


def get_focal_genus() -> str:
    """Get the focal species genus name (e.g., 'Diptychus')."""
    config = load_paths()
    genus = config.get("species", {}).get("genus")
    if genus:
        return genus
    return get_focal_species().split("_")[0]


def get_focal_species_epithet() -> str:
    """Get the focal species epithet (e.g., 'maculatus')."""
    config = load_paths()
    epithet = config.get("species", {}).get("species_epithet")
    if epithet:
        return epithet
    return get_focal_species().split("_")[1]


def get_focal_tax_id() -> Optional[str]:
    """Get the focal species NCBI taxonomy ID (e.g., '2743814')."""
    config = load_paths()
    return config.get("species", {}).get("tax_id")


def get_focal_abbreviation() -> str:
    """Get the focal species 3-letter abbreviation (e.g., 'Dip')."""
    config = load_paths()
    abbr = config.get("species", {}).get("abbreviation")
    if abbr:
        return abbr
    return get_focal_genus()[:3]


def get_focal_display_name(abbreviated: bool = False) -> str:
    """
    Get the focal species display name with spaces.

    Args:
        abbreviated: If True, abbreviate genus to first letter (e.g., 'D. maculatus')

    Returns:
        Display name string
    """
    genus = get_focal_genus()
    epithet = get_focal_species_epithet()
    if abbreviated:
        return f"{genus[0]}. {epithet}"
    return f"{genus} {epithet}"


def get_orgdb_package_name() -> str:
    """Get the OrgDb R package name (e.g., 'org.Dmaculatus.eg.db')."""
    config = load_paths()
    org_db_path = config.get("database", {}).get("org_db")
    if org_db_path:
        return os.path.basename(org_db_path)
    genus = get_focal_genus()
    epithet = get_focal_species_epithet()
    return f"org.{genus[0]}{epithet}.eg.db"


def get_focal_mcscanx_prefix() -> str:
    """Get the focal species 2-letter MCScanX prefix (e.g., 'dm')."""
    config = load_paths()
    prefix = config.get("species", {}).get("mcscanx_prefix")
    if prefix:
        return prefix
    genus = get_focal_genus()
    epithet = get_focal_species_epithet()
    return (genus[0] + epithet[0]).lower()


def get_species_prefix(species_name: str) -> str:
    """
    Get a 2-letter MCScanX prefix from any species name.

    Args:
        species_name: e.g., 'Oxygymnocypris_stewartii'

    Returns:
        Prefix, e.g., 'os'
    """
    parts = species_name.split("_")
    return (parts[0][0] + parts[1][0]).lower()


def get_mcscanx_comparisons() -> List[Dict[str, str]]:
    """
    Get MCScanX comparison configurations from YAML.

    Returns:
        List of dicts, each with keys: name, focal_species, comp_species,
        focal_prefix, comp_prefix
    """
    config = load_paths()
    comps_raw = config.get("mcscanx", {}).get("comparisons", [])

    focal_prefix = get_focal_mcscanx_prefix()
    result = []
    for comp in comps_raw:
        comp_species = comp["comparison"]
        comp_prefix = get_species_prefix(comp_species)
        result.append({
            "name": f"{focal_prefix}_{comp_prefix}",
            "focal_species": comp["focal"],
            "comp_species": comp_species,
            "focal_prefix": focal_prefix,
            "comp_prefix": comp_prefix,
        })
    return result


def get_version() -> str:
    """
    Get the current OrthoFinder results version from configuration.

    Returns:
        Version string (e.g., 'Jan08')
    """
    config = load_paths()
    return config.get("version", "unknown")


def ensure_output_dirs() -> None:
    """
    Create all output directories defined in configuration.

    This function creates directories that don't exist, making it
    safe to call at notebook startup.
    """
    dir_paths = [
        get_path("outputs.kaks.base_dir"),
        get_path("outputs.kaks.codon_alignments"),
        get_path("outputs.kaks.yn00_results"),
        get_path("outputs.kaks.codeml_results"),
        get_path("outputs.kaks.integration"),
        get_path("outputs.kaks.figures"),
        get_path("outputs.go_enrichment"),
        get_path("outputs.ko_enrichment"),
        get_path("outputs.figures"),
        get_path("outputs.tables"),
        get_path("cafe.base_dir"),
        get_path("cafe.single_lambda.dir"),
        get_path("cafe.diptychus_analysis"),
    ]

    for dir_path in dir_paths:
        if dir_path:
            os.makedirs(dir_path, exist_ok=True)


def verify_paths(categories: Optional[List[str]] = None) -> Dict[str, Dict[str, bool]]:
    """
    Verify that configured paths exist.

    Args:
        categories: Optional list of categories to check
                   (e.g., ['orthofinder', 'cafe']). If None, checks all.

    Returns:
        Dictionary with verification results (path -> exists)

    Example:
        >>> results = verify_paths(['orthofinder'])
        >>> for path_key, exists in results['orthofinder'].items():
        ...     status = 'OK' if exists else 'MISSING'
        ...     print(f"  [{status}] {path_key}")
    """
    config = load_paths()

    def check_paths(d: Dict, prefix: str = "") -> Dict[str, bool]:
        path_results = {}
        for key, value in d.items():
            full_key = f"{prefix}.{key}" if prefix else key
            if isinstance(value, str) and ("/" in value or "\\" in value):
                path_results[full_key] = os.path.exists(value)
            elif isinstance(value, dict):
                path_results.update(check_paths(value, full_key))
        return path_results

    results: Dict[str, Dict[str, bool]] = {}

    if categories:
        for cat in categories:
            if cat in config:
                results[cat] = check_paths(config[cat], cat)
    else:
        for cat in config:
            if isinstance(config[cat], dict):
                results[cat] = check_paths(config[cat], cat)

    return results


def print_config_summary() -> None:
    """Print a summary of the loaded configuration."""
    config = load_paths()

    print("=" * 60)
    print("PATH CONFIGURATION SUMMARY")
    print("=" * 60)
    print(f"Config file: {_config_path}")
    print(f"Version: {get_version()}")
    print(f"Project root: {get_path('base.project_root')}")
    print()
    print("Key paths:")
    print(f"  OrthoFinder results: {get_path('orthofinder.results_dir')}")
    print(f"  CAFE base: {get_path('cafe.base_dir')}")
    print(f"  Focal species: {get_focal_species()}")
    print(f"  Species count: {len(get_species_list())}")
    print("=" * 60)


def get_paths_dict() -> Dict[str, str]:
    """
    Get a flat dictionary of common paths for backward compatibility.

    This provides a PATHS-style dictionary similar to the pattern used
    in 2-ultimate-KaKs-analysis.qmd.

    Returns:
        Dictionary with common path keys
    """
    return {
        "orthofinder_base": get_path("orthofinder.results_dir"),
        "orthogroups_tsv": get_path("orthofinder.orthogroups_tsv"),
        "orthogroups_count": get_path("orthofinder.orthogroups_count"),
        "msa_dir": get_path("orthofinder.msa_dir"),
        "orthogroup_sequences": get_path("orthofinder.orthogroup_sequences"),
        "species_tree": get_path("orthofinder.species_tree"),
        "cafe_base": get_path("cafe.base_dir"),
        "cafe_input_tree": get_path("cafe.input_tree"),
        "cafe_family_results": get_path("cafe.single_lambda.family_results"),
        "cafe_change": get_path("cafe.single_lambda.change_tab"),
        "cds_base": get_path("cds.base_dir"),
        "database_base": get_path("database.genome_comparison"),
        "output_figures": get_path("outputs.figures"),
        "output_tables": get_path("outputs.tables"),
    }


# Run tests if executed directly
if __name__ == "__main__":
    print("Testing configuration loader...")
    print()

    try:
        print_config_summary()
        print()

        print("Testing get_path():")
        print(f"  orthofinder.results_dir = {get_path('orthofinder.results_dir')}")
        print(
            f"  orthofinder.orthogroups_tsv = {get_path('orthofinder.orthogroups_tsv')}"
        )
        print(
            f"  cafe.single_lambda.family_results = {get_path('cafe.single_lambda.family_results')}"
        )
        print()

        print("Testing get_species_list():")
        species = get_species_list()
        print(f"  Found {len(species)} species")
        print(f"  First: {species[0] if species else 'N/A'}")
        print(f"  Focal: {get_focal_species()}")
        print()

        print("Testing verify_paths():")
        results = verify_paths(["orthofinder", "cafe"])
        for category, paths in results.items():
            print(f"\n  {category}:")
            ok_count = sum(1 for v in paths.values() if v)
            missing_count = sum(1 for v in paths.values() if not v)
            print(f"    OK: {ok_count}, Missing: {missing_count}")
            for path_key, exists in list(paths.items())[:5]:
                status = "OK" if exists else "MISSING"
                print(f"    [{status}] {path_key}")
            if len(paths) > 5:
                print(f"    ... and {len(paths) - 5} more")

        print()
        print("All tests passed!")

    except Exception as e:
        print(f"Error: {e}")
        raise
