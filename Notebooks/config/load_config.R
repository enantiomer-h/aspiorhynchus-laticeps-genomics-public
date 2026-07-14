#' Path Configuration Loader for Comparative Genomics Notebooks
#'
#' This module provides functions to load and resolve path configurations
#' from the central YAML configuration file (paths.yaml).
#'
#' Usage:
#'   source("./config/load_config.R")
#'
#'   # Load all paths
#'   paths <- load_paths()
#'
#'   # Access paths via list
#'   orthogroups_file <- paths$orthofinder$orthogroups_tsv
#'
#'   # Or use convenience function with dot notation
#'   orthogroups_file <- get_path("orthofinder.orthogroups_tsv")
#'
#' Template Variables:
#'   The YAML configuration supports template variables in the format ${section.key}
#'   These are automatically resolved when the configuration is loaded.
#'

# Check for yaml package
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required. Install with: install.packages('yaml')")
}

library(yaml)

# Module-level cache for configuration
.config_env <- new.env(parent = emptyenv())
.config_env$cache <- NULL
.config_env$path <- NULL

#' Null coalescing operator
#'
#' Returns the first argument if not NULL, otherwise the second.
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Find the configuration file
#'
#' Searches for paths.yaml in standard locations.
#'
#' @return Path to the configuration file
#' @keywords internal
find_config_file <- function() {
  search_paths <- c(
    "./config/paths-set2.yaml",
    "./config/paths.yaml",
    "../config/paths-set2.yaml",
    "../config/paths.yaml",
    "/home/jovyan/Notebooks/config/paths-set2.yaml",
    "/home/jovyan/Notebooks/config/paths.yaml"
  )

  # Check environment variable first
  env_path <- Sys.getenv("GENOMICS_CONFIG_PATH", "")
  if (nzchar(env_path)) {
    search_paths <- c(env_path, search_paths)
  }

  for (path in search_paths) {
    if (file.exists(path)) {
      return(normalizePath(path))
    }
  }

  # Try relative to this script's location
  script_path <- tryCatch({
    # When sourced, get the script location
    if (sys.nframe() > 0) {
      frame_files <- sapply(sys.frames(), function(f) {
        tryCatch(get("ofile", envir = f), error = function(e) NULL)
      })
      frame_files <- frame_files[!sapply(frame_files, is.null)]
      if (length(frame_files) > 0) {
        dirname(frame_files[[1]])
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, error = function(e) NULL)

  if (!is.null(script_path)) {
    for (filename in c("paths-set2.yaml", "paths.yaml")) {
      config_path <- file.path(script_path, filename)
      if (file.exists(config_path)) {
        return(normalizePath(config_path))
      }
    }
  }

  # Also try current working directory (paths-set2.yaml first)
  for (filename in c("paths-set2.yaml", "paths.yaml")) {
    if (file.exists(filename)) {
      return(normalizePath(filename))
    }
  }

  stop(paste(
    "Configuration file not found. Searched:",
    paste(search_paths, collapse = ", "),
    "\nSet GENOMICS_CONFIG_PATH environment variable or place paths-set2.yaml in config directory."
  ))
}

#' Get a nested value from a list using dot notation
#'
#' @param lst List to search
#' @param key_path Dot-separated path (e.g., 'orthofinder.results_dir')
#' @return The value if found, NULL otherwise
#' @keywords internal
get_nested_value <- function(lst, key_path) {
  parts <- strsplit(key_path, "\\.")[[1]]
  current <- lst

  tryCatch({
    for (part in parts) {
      if (is.list(current) && part %in% names(current)) {
        current <- current[[part]]
      } else {
        return(NULL)
      }
    }
    return(current)
  }, error = function(e) {
    return(NULL)
  })
}

#' Set a nested value in a list using dot notation
#'
#' @param lst List to modify
#' @param key_path Dot-separated path
#' @param value Value to set
#' @return Modified list
#' @keywords internal
set_nested_value <- function(lst, key_path, value) {
  parts <- strsplit(key_path, "\\.")[[1]]

  if (length(parts) == 1) {
    lst[[parts[1]]] <- value
    return(lst)
  }

  if (is.null(lst[[parts[1]]]) || !is.list(lst[[parts[1]]])) {
    lst[[parts[1]]] <- list()
  }

  remaining_path <- paste(parts[-1], collapse = ".")
  lst[[parts[1]]] <- set_nested_value(lst[[parts[1]]], remaining_path, value)
  return(lst)
}

#' Resolve template variables in a string
#'
#' Handles ${section.key} and ${section.nested.key} syntax.
#'
#' @param value String potentially containing template variables
#' @param config Original configuration list
#' @param resolved Already resolved values for reference
#' @return String with all template variables resolved
#' @keywords internal
resolve_template <- function(value, config, resolved) {
  if (!is.character(value) || length(value) != 1) {
    return(value)
  }

  # Pattern matches ${...} with nested dot notation
  pattern <- "\\$\\{([a-zA-Z_][a-zA-Z0-9_]*(?:\\.[a-zA-Z_][a-zA-Z0-9_]*)*)\\}"

  max_iterations <- 10
  prev_value <- NULL

  for (i in seq_len(max_iterations)) {
    if (identical(prev_value, value)) {
      break
    }
    prev_value <- value

    # Find all matches
    matches <- gregexpr(pattern, value, perl = TRUE)[[1]]

    if (matches[1] == -1) {
      break
    }

    # Extract match positions
    match_starts <- as.vector(matches)
    match_lengths <- attr(matches, "match.length")

    # Process matches in reverse order to preserve positions
    for (j in rev(seq_along(match_starts))) {
      start <- match_starts[j]
      len <- match_lengths[j]
      matched_str <- substr(value, start, start + len - 1)

      # Extract variable path (remove ${ and })
      var_path <- substr(matched_str, 3, nchar(matched_str) - 1)

      # Try to resolve from resolved first, then config
      replacement <- NULL

      # Try resolved values
      resolved_val <- get_nested_value(resolved, var_path)
      if (!is.null(resolved_val) && is.character(resolved_val) && length(resolved_val) == 1) {
        replacement <- resolved_val
      }

      # Try original config if not found
      if (is.null(replacement)) {
        config_val <- get_nested_value(config, var_path)
        if (!is.null(config_val) && is.character(config_val) && length(config_val) == 1) {
          replacement <- config_val
        }
      }

      if (!is.null(replacement)) {
        value <- paste0(
          substr(value, 1, start - 1),
          replacement,
          substr(value, start + len, nchar(value))
        )
      }
    }
  }

  return(value)
}

#' Recursively resolve all template variables in a list
#'
#' @param lst List to resolve
#' @param config Original full configuration
#' @param resolved Already resolved values (environment for mutation)
#' @param path_prefix Current path for nested resolution
#' @return List with all templates resolved
#' @keywords internal
resolve_list <- function(lst, config, resolved, path_prefix = "") {
  result <- list()

  for (key in names(lst)) {
    value <- lst[[key]]
    current_path <- if (nzchar(path_prefix)) paste(path_prefix, key, sep = ".") else key

    if (is.list(value) && !is.null(names(value))) {
      # Named list (dictionary-like)
      result[[key]] <- resolve_list(value, config, resolved, current_path)
    } else if (is.character(value)) {
      if (length(value) == 1) {
        result[[key]] <- resolve_template(value, config, as.list(resolved))
      } else {
        result[[key]] <- sapply(value, function(v) {
          resolve_template(v, config, as.list(resolved))
        }, USE.NAMES = FALSE)
      }
    } else {
      result[[key]] <- value
    }

    # Update resolved for subsequent references
    assign(current_path, result[[key]], envir = resolved)
  }

  return(result)
}

#' Load and resolve path configuration from YAML file
#'
#' @param config_path Optional explicit path to config file
#' @param force_reload If TRUE, reload even if cached
#' @return List with all paths resolved
#' @export
#' @examples
#' paths <- load_paths()
#' orthogroups_file <- paths$orthofinder$orthogroups_tsv
load_paths <- function(config_path = NULL, force_reload = FALSE) {
  # Use cached version if available
  if (!is.null(.config_env$cache) && !force_reload && is.null(config_path)) {
    return(.config_env$cache)
  }

  # Find config file
  if (!is.null(config_path)) {
    if (!file.exists(config_path)) {
      stop(paste("Configuration file not found:", config_path))
    }
    path <- normalizePath(config_path)
  } else {
    path <- find_config_file()
  }

  # Load YAML
  config <- yaml::read_yaml(path)

  if (is.null(config)) {
    stop(paste("Configuration file is empty:", path))
  }

  # Resolve all templates
  resolved <- new.env(hash = TRUE, parent = emptyenv())
  result <- resolve_list(config, config, resolved)

  # Cache the result
  .config_env$cache <- result
  .config_env$path <- path

  return(result)
}

#' Get a specific path by dot-notation key
#'
#' @param key_path Dot-separated path to the configuration value
#'   e.g., "orthofinder.orthogroups_tsv"
#' @param default Default value if key not found
#' @return The path string, or default if not found
#' @export
#' @examples
#' get_path("orthofinder.orthogroups_tsv")
#' # Returns: "/home/jovyan/Outputs/OrthoFinder/Results_Jan09/Orthogroups/Orthogroups.tsv"
get_path <- function(key_path, default = NULL) {
  config <- load_paths()
  result <- get_nested_value(config, key_path)
  return(result %||% default)
}

#' Get dictionary of all CDS file paths keyed by species name
#'
#' @return Named list mapping species names to CDS file paths
#' @export
get_cds_files <- function() {
  config <- load_paths()
  cds_config <- config$cds

  if (is.null(cds_config)) {
    return(list())
  }

  # Filter out non-species keys (like base_dir)
  exclude_keys <- c("base_dir")
  species_files <- cds_config[!names(cds_config) %in% exclude_keys]

  # Filter to only string values (actual paths)
  species_files <- species_files[sapply(species_files, function(x) {
    is.character(x) && length(x) == 1
  })]

  return(species_files)
}

#' Get list of all species from configuration
#'
#' @return Character vector of species names
#' @export
get_species_list <- function() {
  config <- load_paths()
  species_list <- config$species$all

  if (is.null(species_list)) {
    return(character(0))
  }

  return(species_list)
}

#' Get the focal species name from configuration
#'
#' @return Focal species name
#' @export
get_focal_species <- function() {
  config <- load_paths()
  focal <- config$species$focal
  return(focal %||% "Diptychus_maculatus")
}

#' Get the focal species genus name
#'
#' @return Genus name (e.g., "Diptychus")
#' @export
get_focal_genus <- function() {
  config <- load_paths()
  genus <- config$species$genus
  if (!is.null(genus)) return(genus)
  focal <- get_focal_species()
  return(strsplit(focal, "_")[[1]][1])
}

#' Get the focal species epithet
#'
#' @return Species epithet (e.g., "maculatus")
#' @export
get_focal_species_epithet <- function() {
  config <- load_paths()
  epithet <- config$species$species_epithet
  if (!is.null(epithet)) return(epithet)
  focal <- get_focal_species()
  return(strsplit(focal, "_")[[1]][2])
}

#' Get the focal species NCBI taxonomy ID
#'
#' @return Taxonomy ID string (e.g., "2743814"), or NULL if not configured
#' @export
get_focal_tax_id <- function() {
  config <- load_paths()
  return(config$species$tax_id)
}

#' Get the focal species abbreviation (3-letter)
#'
#' @return Abbreviation (e.g., "Dip")
#' @export
get_focal_abbreviation <- function() {
  config <- load_paths()
  abbr <- config$species$abbreviation
  if (!is.null(abbr)) return(abbr)
  genus <- get_focal_genus()
  return(substr(genus, 1, 3))
}

#' Get the focal species display name
#'
#' @param abbreviated If TRUE, abbreviate genus to first letter (e.g., "D. maculatus")
#' @return Display name with spaces instead of underscores
#' @export
get_focal_display_name <- function(abbreviated = FALSE) {
  genus <- get_focal_genus()
  epithet <- get_focal_species_epithet()
  if (abbreviated) {
    return(paste0(substr(genus, 1, 1), ". ", epithet))
  }
  return(paste(genus, epithet))
}

#' Get the OrgDb package name for the focal species
#'
#' @return Package name (e.g., "org.Dmaculatus.eg.db")
#' @export
get_orgdb_package_name <- function() {
  config <- load_paths()
  org_db_path <- config$database$org_db
  if (!is.null(org_db_path)) {
    return(basename(org_db_path))
  }
  genus <- get_focal_genus()
  epithet <- get_focal_species_epithet()
  return(paste0("org.", substr(genus, 1, 1), epithet, ".eg.db"))
}

#' Get the focal species MCScanX 2-letter prefix
#'
#' @return Prefix (e.g., "dm")
#' @export
get_focal_mcscanx_prefix <- function() {
  config <- load_paths()
  prefix <- config$species$mcscanx_prefix
  if (!is.null(prefix)) return(prefix)
  genus <- get_focal_genus()
  epithet <- get_focal_species_epithet()
  return(tolower(paste0(substr(genus, 1, 1), substr(epithet, 1, 1))))
}

#' Get a 2-letter MCScanX prefix from an arbitrary species name
#'
#' @param species_name Species name (e.g., "Oxygymnocypris_stewartii")
#' @return Prefix (e.g., "os")
#' @export
get_species_prefix <- function(species_name) {
  parts <- strsplit(species_name, "_")[[1]]
  return(tolower(paste0(substr(parts[1], 1, 1), substr(parts[2], 1, 1))))
}

#' Get MCScanX comparison configurations from YAML
#'
#' Each element has: name, focal_species, comp_species, focal_prefix, comp_prefix
#'
#' @return List of comparison configurations
#' @export
get_mcscanx_comparisons <- function() {
  config <- load_paths()
  comps_raw <- config$mcscanx$comparisons
  if (is.null(comps_raw)) return(list())

  focal_prefix <- get_focal_mcscanx_prefix()

  lapply(comps_raw, function(comp) {
    comp_species <- comp$comparison
    comp_prefix <- get_species_prefix(comp_species)
    list(
      name = paste0(focal_prefix, "_", comp_prefix),
      focal_species = comp$focal,
      comp_species = comp_species,
      focal_prefix = focal_prefix,
      comp_prefix = comp_prefix
    )
  })
}

#' Get the current OrthoFinder results version
#'
#' @return Version string (e.g., 'Jan08')
#' @export
get_version <- function() {
  config <- load_paths()
  return(config$version %||% "unknown")
}

#' Create all output directories defined in configuration
#'
#' @export
ensure_output_dirs <- function() {
  dir_paths <- c(
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
    get_path("cafe.diptychus_analysis")
  )

  for (dir_path in dir_paths) {
    if (!is.null(dir_path) && nzchar(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
}

#' Verify that configured paths exist
#'
#' @param categories Optional character vector of categories to check
#'   (e.g., c("orthofinder", "cafe")). If NULL, checks all.
#' @return Named list with verification results (TRUE/FALSE for each path)
#' @export
verify_paths <- function(categories = NULL) {
  config <- load_paths()

  check_paths <- function(lst, prefix = "") {
    path_results <- list()

    for (key in names(lst)) {
      value <- lst[[key]]
      full_key <- if (nzchar(prefix)) paste(prefix, key, sep = ".") else key

      if (is.character(value) && length(value) == 1 && grepl("/", value)) {
        path_results[[full_key]] <- file.exists(value)
      } else if (is.list(value) && !is.null(names(value))) {
        path_results <- c(path_results, check_paths(value, full_key))
      }
    }

    return(path_results)
  }

  if (!is.null(categories)) {
    results <- list()
    for (cat in categories) {
      if (cat %in% names(config)) {
        results[[cat]] <- check_paths(config[[cat]], cat)
      }
    }
    return(results)
  } else {
    results <- list()
    for (cat in names(config)) {
      if (is.list(config[[cat]]) && !is.null(names(config[[cat]]))) {
        results[[cat]] <- check_paths(config[[cat]], cat)
      }
    }
    return(results)
  }
}

#' Print a summary of the loaded configuration
#'
#' @export
print_config_summary <- function() {
  config <- load_paths()

  cat(strrep("=", 60), "\n")
  cat("PATH CONFIGURATION SUMMARY\n")
  cat(strrep("=", 60), "\n")
  cat(sprintf("Config file: %s\n", .config_env$path %||% "N/A"))
  cat(sprintf("Version: %s\n", get_version()))
  cat(sprintf("Project root: %s\n", get_path("base.project_root")))
  cat("\n")
  cat("Key paths:\n")
  cat(sprintf("  OrthoFinder results: %s\n", get_path("orthofinder.results_dir")))
  cat(sprintf("  CAFE base: %s\n", get_path("cafe.base_dir")))
  cat(sprintf("  Focal species: %s\n", get_focal_species()))
  cat(sprintf("  Species count: %d\n", length(get_species_list())))
  cat(strrep("=", 60), "\n")
}

#' Get a flat list of common paths for backward compatibility
#'
#' @return Named list with common path keys
#' @export
get_paths_list <- function() {
  list(
    orthofinder_base = get_path("orthofinder.results_dir"),
    orthogroups_tsv = get_path("orthofinder.orthogroups_tsv"),
    orthogroups_count = get_path("orthofinder.orthogroups_count"),
    msa_dir = get_path("orthofinder.msa_dir"),
    orthogroup_sequences = get_path("orthofinder.orthogroup_sequences"),
    species_tree = get_path("orthofinder.species_tree"),
    cafe_base = get_path("cafe.base_dir"),
    cafe_input_tree = get_path("cafe.input_tree"),
    cafe_family_results = get_path("cafe.single_lambda.family_results"),
    cafe_change = get_path("cafe.single_lambda.change_tab"),
    cds_base = get_path("cds.base_dir"),
    database_base = get_path("database.genome_comparison"),
    output_figures = get_path("outputs.figures"),
    output_tables = get_path("outputs.tables")
  )
}

# Try to load configuration when script is sourced (silent failure)
tryCatch({
  .default_config <- load_paths()
  message(sprintf("Configuration loaded: OrthoFinder Results_%s", get_version()))
}, error = function(e) {
  # Silent failure - user can call load_paths() explicitly with path
})
