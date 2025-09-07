# Comparative Genomics Analysis Project

This repository contains a comprehensive comparative genomics analysis pipeline for studying evolutionary relationships across multiple fish species, with particular focus on *Aspiorhynchus laticeps*. The project implements a complete workflow from environment setup through advanced phylogenetic and synteny analysis.

## Project Overview

The analysis pipeline identifies and characterizes orthogroups (groups of homologous genes) across 13 fish species, providing insights into:
- **Evolutionary relationships** through phylogenetic tree construction
- **Gene family dynamics** through orthogroup analysis  
- **Functional enrichment** through GO and KEGG pathway analysis
- **Ancestral state reconstruction** using PastML
- **Syntenic relationships** through MCScanX analysis

### Key Results
- Comprehensive functional enrichment analysis for species-specific gene families
- Phylogenetic reconstruction with ancestral state analysis
- Synteny analysis revealing chromosomal evolution patterns

## Environment Setup

This project uses a containerized environment with Docker to ensure reproducibility and consistent tool availability.

### Prerequisites
- Docker and Docker Compose
- At least 12GB RAM available for containers
- Optional: NVIDIA GPU support for accelerated analysis

### Quick Start

1. **Clone the repository** (if applicable):
   ```bash
   # Navigate to your project directory
   cd /path/to/ComparativeGenomics4AL_Publication
   ```

2. **Initialize the environment**:
   ```bash
   # For CPU-only analysis
   ./init.sh
   
   # For GPU-accelerated analysis (requires NVIDIA Docker)
   ./init.sh GPU
   ```

3. **Access Jupyter Lab**:
   - The container will start automatically
   - Access the notebook server at `http://localhost:38888` (check console for token)
   - Use port `30000` if specified in your configuration

### Manual Environment Setup

The Docker environment includes multiple specialized conda environments:

| Environment | Purpose | Key Tools |
|-------------|---------|-----------|
| `env_biopython` | Python analysis | BioPython, pandas, matplotlib |
| `env_orthofinder` | Orthology analysis | OrthoFinder, DIAMOND |
| `bioconductor-clusterprofiler_env` | R enrichment analysis | clusterProfiler, GO.db |
| `eggnog_env` | Functional annotation | EggNOG-mapper |
| `pastml_env` | Ancestral reconstruction | PastML |
| `gff2bed_env` | Format conversion | BEDOPS utilities |

### Available Jupyter Kernels
- **Python (Biopython)**: `env_biopython` - Main analysis kernel
- **R clusterprofiler**: `ir_clusterprofiler_4100_env` - GO/KEGG enrichment analysis
- **R clusterprofiler (legacy)**: `ir_clusterprofiler_481_env` - Alternative enrichment analysis kernel
- **Base R**: `ir` - General R analysis
- **Julia**: `julia-1.7` - High-performance computing
- **Base Python**: `python3` - Standard Python kernel

#### Setting Default Kernel
You can change the default kernel by editing `.jupyter/jupyter_notebook_config.py`.

```py
## The name of the default kernel to start
#  Default: 'python3'
c.MultiKernelManager.default_kernel_name = 'python3'
# c.MultiKernelManager.default_kernel_name = 'env_biopython'
# c.MultiKernelManager.default_kernel_name = 'ir_clusterprofiler_4100_env'
```

## Analysis Pipeline

The project follows a structured 7-step analysis pipeline:

### Step 0: Environment Construction
**Notebook**: `0-ultimate-env_construction.qmd`
**Title**: "Ultimate Environment Construction Guide - Complete setup for comparative genomics analysis"
- Comprehensive environment setup for all bioinformatics tools and dependencies
- Installation and configuration of multiple specialized conda environments
- Jupyter kernel setup and verification
- Tool accessibility testing and troubleshooting guides

### Step 1: OrthoFinder Preparation
**Notebooks**: 
- `1-ultimate-OrthoFinder_preparation-step1.qmd` - "Protein Sequence Extraction"
- `1-ultimate-OrthoFinder_preparation-step2.qmd` - "File Preprocessing and Header Standardization"

**Process**:
- **Step 1**: Extract protein sequences (FAA files) from genomic sequences (FNA files) using GFF3 annotations
- **Step 2**: Standardize sequence identifiers between FNA and GFF3 files to resolve header mismatches
- Address common issues: case sensitivity differences, naming convention variations, complex header formats
- Prepare properly formatted and validated input files for OrthoFinder analysis

### Step 2: OrthoFinder Analysis  
**Notebook**: `2-ultimate-AL_OrthoFinder.qmd`
**Title**: "Ultimate OrthoFinder Analysis - Comparative genomics analysis with comprehensive troubleshooting guide"

**Process**:
- Run comparative genomics analysis on fish species (14 species initially, 13 species final after troubleshooting)
- Perform all-vs-all sequence comparisons using DIAMOND
- Identify orthogroups and construct species phylogeny with comprehensive error handling
- Generate orthology matrices with detailed troubleshooting documentation

**Key Outputs**:
- Species tree with gene duplication events
- Orthogroup gene count matrices  
- Statistical summaries of orthology relationships
- **Results**: 52,388 orthogroups identified from 644,648 genes across 13 species

### Step 3: Visualization & Statistical Analysis
**Notebook**: `3-ultimate-Visualization4OrthoFinder.qmd`
**Title**: "Ultimate OrthoFinder Visualization Analysis - Comprehensive visualization suite for comparative genomics results"

**Visualizations**:
- **Phylogenetic trees** with branch support values and publication-quality formatting
- **Heatmaps** showing gene family size variations with advanced filtering options
- **Upset plots** for orthogroup overlaps with statistical analysis
- **Statistical plots** for orthogroup distributions and duplication events
- **Species overlap analysis** with normalized gene counts and comparative ratios
- **Publication-ready outputs** with 300 DPI resolution and professional typography

### Step 4: Functional Enrichment Analysis
**Notebooks**:
- `4-ultimate-EnrichmentAnalysis-preparing.qmd` - "Annotation and Database Setup"
- `4-ultimate-EnrichmentAnalysis4SpecificOGs.qmd` - "Enrichment Analysis for Specific Orthogroups"

**Preparation Process**:
- Generate functional annotations using eggnog-mapper
- Create custom organism databases (OrgDb) for species-specific analysis
- Extract GO and KEGG pathway information with validation

**Analysis Categories**:
- **Increased gene families**: Where *Aspiorhynchus laticeps* has maximum gene count
- **Decreased gene families**: Where *Aspiorhynchus laticeps* has minimum gene count  
- **Species-specific families**: Unique to *Aspiorhynchus laticeps* (non-zero exclusive)

**Methods**:
- Gene Ontology (GO) enrichment analysis (BP, CC, MF) with clusterProfiler
- KEGG Orthology (KO) pathway enrichment with custom databases
- Statistical significance testing with multiple correction methods
- Publication-quality barplots and dotplots for each category

### Step 5: Ancestral State Reconstruction
**Notebook**: `5-ultimate-PastML.qmd`
**Title**: "Ultimate PastML Analysis for Phylogenetic Trees"

**Process**:
- Enhanced phylogenetic tree visualization from OrthoFinder results
- Comprehensive gene duplication analysis with statistical plotting
- PastML-based ancestral state reconstruction for orthogroups
- Batch processing workflows for large-scale analysis with automatic file naming
- Integration with OrthoFinder species trees and enhanced error handling
- Publication-ready tree visualizations with automatic dual-format output

### Step 6: Synteny Analysis
**Notebook**: `6-ultimate-MCScanX.qmd`
**Title**: "Ultimate MCScanX Synteny Analysis - Comparative Genomics"

**Process**:
- Comprehensive synteny analysis using MCScanX for syntenic block identification
- Pairwise BLAST searches between genomes with optimized parameters
- Multiple visualization approaches: dot plots, dual synteny plots, circle plots, bar plots
- Chromosomal collinearity assessment and rearrangement detection
- Gene duplication mechanism analysis and evolutionary relationship mapping
- HTML reports for detailed synteny examination and interactive visualization

## Project Structure

```
ComparativeGenomics4AL_Publication/
├── README.md                    # This file
├── CLAUDE.md                   # Claude Code project instructions
├── init.sh                     # Environment initialization script
├── docker-compose.yml          # Docker container configuration
├── .build/
│   └── Dockerfile             # Container build instructions
├── .jupyter/                   # Jupyter configuration files
├── Notebooks/                  # Analysis notebooks and outputs
│   ├── 0-ultimate-env_construction.qmd
│   ├── 1-ultimate-OrthoFinder_preparation-step1.qmd
│   ├── 1-ultimate-OrthoFinder_preparation-step2.qmd
│   ├── 2-ultimate-AL_OrthoFinder.qmd
│   ├── 3-ultimate-Visualization4OrthoFinder.qmd
│   ├── 4-ultimate-EnrichmentAnalysis-preparing.qmd
│   ├── 4-ultimate-EnrichmentAnalysis4SpecificOGs.qmd
│   ├── 5-ultimate-PastML.qmd
│   ├── 6-ultimate-MCScanX.qmd
│   ├── Figures/               # Generated visualization outputs
│   │   ├── GO_Enrichment/     # GO enrichment plots
│   │   ├── KO_Enrichment/     # KEGG pathway plots
│   │   ├── OG_bars/          # Orthogroup distribution plots
│   │   └── *.pdf/*.png       # Tree and statistical plots
│   └── Tables/                # Analysis result tables
│       ├── GO_enrichment_*.csv
│       ├── KO_enrichment_*.csv
│       └── *.tsv             # Orthogroup and ratio tables
└── DB/                        # Input genomic data (not shown)
    └── GENOME_Comparison/FAA/ # Protein sequence files
```

## Usage Instructions

### Running the Complete Pipeline

1. **Start the environment**:
   ```bash
   ./init.sh  # or ./init.sh GPU for GPU support
   ```

2. **Access Jupyter Lab**:
   - Navigate to `http://localhost:38888` (check console for token)
   - Open the notebook files in numerical order (0-6)

3. **Execute notebooks sequentially**:
   - Follow the numbered sequence for proper dependency handling
   - Each notebook contains detailed documentation and error handling
   - Results are automatically saved to `Notebooks/Figures/` and `Notebooks/Tables/`

### Individual Analysis Steps

#### Environment Setup Only
```bash
# Access the container
docker exec -it gpu-jupyter bash
# Activate specific environment
conda activate env_biopython
```

#### OrthoFinder Analysis Only
```bash
conda run -n env_orthofinder orthofinder -f ./DB/GENOME_Comparison/FAA/
```

#### Custom Analysis
- Use the `env_biopython` kernel for Python-based analysis
- Use the `ir_clusterprofiler_4100_env` kernel for R-based enrichment analysis
- Refer to existing notebooks for code examples and best practices

## Key Dependencies

### Core Bioinformatics Tools
- **OrthoFinder** (v2.5+): Orthology inference and comparative genomics
- **DIAMOND** (v2.0+): Fast sequence alignment for large datasets  
- **PastML**: Maximum likelihood ancestral state reconstruction
- **MCScanX**: Synteny analysis and visualization
- **EggNOG-mapper**: Functional annotation

### Python Libraries  
- **BioPython**: Sequence analysis and file parsing
- **pandas/numpy**: Data manipulation and numerical analysis
- **matplotlib/seaborn**: Statistical visualization
- **gffutils**: GFF3 file processing

### R Libraries
- **clusterProfiler**: GO and KEGG enrichment analysis
- **ggplot2**: Advanced statistical plotting
- **dplyr**: Data manipulation and filtering

## Output Files

### Primary Results
- `Notebooks/Figures/SpeciesTree_rooted.pdf` - Main phylogenetic tree
- `Notebooks/Tables/orthologs-ratio.tsv` - Orthogroup quantification matrix
- `Notebooks/Tables/*enrichment*.csv` - GO and KEGG enrichment results

### Visualization Outputs
- Species overlap heatmaps and normalized distributions
- Individual orthogroup bar plots (300+ files in `OG_bars/`)
- Upset plots showing orthogroup intersections
- Phylogenetic trees with gene duplication annotations

## Figure Generation Reference

This section provides a comprehensive reference for all plotting functions and their output figures across the analysis notebooks.

### 3-ultimate-Visualization4OrthoFinder.qmd

#### Enhanced Phylogenetic Tree Functions
- **`plot_phylogenetic_tree_with_lengths()`** (lines 49-234)
  - **Output**: Publication-quality species phylogenetic trees with branch lengths
  - **Files**: `./Notebooks/Figures/SpeciesTree_*.png|pdf`
  - **Features**: 300 DPI resolution, Times New Roman typography, dual PNG/PDF output, automatic scaling
  - **Used for**: Rooted species tree, node-labeled tree, gene duplication tree

#### Advanced R ggtree Visualization Functions
- **`assign_clade_colors()`** (lines 320-388)
  - **Purpose**: Automatic color assignment for phylogenetic clades using clustering algorithms
  - **Methods**: Height-based, depth-based, and support-based clustering

- **`plot_phylogenetic_tree_ggtree()`** (lines 391-652)
  - **Output**: Advanced colorful phylogenetic trees with clustering
  - **Files**: `*_ggtree.png|pdf|svg`
  - **Features**: Multiple layouts (rectangular, circular, fan, radial), automatic clade coloring, publication themes

- **`plot_circular_tree_ggtree()`** (lines 655-664)
  - **Output**: Circular layout phylogenetic trees
  - **Features**: Specialized circular visualization with enhanced tip labeling

- **`plot_fan_tree_ggtree()`** (lines 667-675)
  - **Output**: Fan-shaped phylogenetic trees with gradient coloring
  - **Features**: Gradient branch coloring based on branch lengths

- **`create_tree_comparison_plots()`** (lines 678-705)
  - **Output**: Multiple visualization schemes for comparative analysis
  - **Features**: Auto_Clade, Viridis_Gradient, Plasma_Clade, Branch_Length schemes

#### Orthogroup Analysis Functions  
- **`plot_orthogroup_bar_chart()`** (lines 878-908)
  - **Output**: Gene count bar charts for individual orthogroups
  - **Files**: `./Notebooks/Figures/OG_bars/OG*_OG_genecount_barplot.png|pdf`
  - **Batch Usage**: Generates 300+ individual orthogroup plots

#### Gene Tree Functions
- **`select_representatives()`** (lines 945-968)
  - **Purpose**: Representative clade selection based on tree depth

- **`prune_and_plot_tree()`** (lines 970-1016)
  - **Output**: Pruned gene trees for specific orthogroups
  - **Files**: `./Notebooks/Figures/OG*_gene_tree.png|pdf`
  - **Used for**: Gene tree visualization with depth-based pruning

#### Comprehensive Heatmap Analysis Functions
- **`load_and_process_orthogroup_data()`** (lines 1044-1068)
  - **Purpose**: Centralized data loading with error handling
  - **Features**: Handles TSV files, optional Total column removal

- **`create_publication_heatmap()`** (lines 1069-1544)
  - **Output**: Advanced publication-quality heatmaps with comprehensive styling
  - **Files**: Various heatmap outputs with customizable paths
  - **Features**: 300 DPI resolution, journal-specific themes, statistical overlays, multiple transforms
  - **Advanced Options**: Data transformations, clustering, annotations, custom colormaps

- **`filter_orthogroups_by_variance()`** (lines 1546-1568)
  - **Output**: `heatmap_top_variance.png|pdf`
  - **Purpose**: Variance-based orthogroup filtering

- **`filter_orthogroups_by_species_extreme()`** (lines 1569-1608)
  - **Output**: `heatmap_Aspiorhynchus_max*.png|pdf`, `heatmap_Aspiorhynchus_min*.png|pdf`
  - **Purpose**: Species-specific maximum/minimum gene count analysis

- **`filter_species_specific_orthogroups()`** (lines 1609-1642)
  - **Output**: `heatmap_unique_nonzero_Aspiorhynchus*.png|pdf`
  - **Purpose**: Species-unique orthogroup identification

#### Species Overlap Analysis Functions
- **`create_species_overlap_heatmap()`** (lines 1747-2207)
  - **Output**: `Orthogroups_SpeciesOverlaps.png|pdf`, `Orthogroups_SpeciesOverlaps_Normalized.png|pdf`
  - **Features**: Advanced symmetric matrix visualization, clustering, statistical overlays
  - **Advanced Options**: Multiple transformations, journal-specific styling, diagonal highlighting

- **`_create_overlap_statistics_text()`** (lines 2208-2232)
  - **Purpose**: Helper function for statistical text generation

- **`analyze_species_overlap_statistics()`** (lines 2234-2261)
  - **Purpose**: Statistical summary generation with mean, median, range calculations

#### Gene and Orthogroup Ratio Analysis Functions
- **`create_gene_orthogroup_ratio_plot()`** (lines 2310-2457)
  - **Output**: `orthologs-ratio.png|pdf`, `orthologs-ratio-grouped.png|pdf`
  - **Features**: Multi-type plotting (stacked, grouped, line), professional color palette

- **`calculate_ratio_statistics()`** (lines 2458-2477)
  - **Purpose**: Statistical summaries with coefficient of variation

- **`create_comparative_ratio_plot()`** (lines 2478-2602)
  - **Output**: `orthologs-ratio_comparative.png|pdf`
  - **Features**: Scatter plots and ratio analysis for gene vs. orthogroup relationships

#### Statistical Analysis Functions
- **`plot_duplications_histogram()`** (lines 2666-2698)
  - **Output**: Histograms of gene duplication events
  - **Files**: `./Notebooks/Figures/duplications_histogram.png|pdf`

- **`plot_duplications_box_and_violin()`** (lines 2699-2731)
  - **Output**: Box and violin plots for duplication statistics
  - **Files**: `./Notebooks/Figures/duplications_box_violin.png|pdf`

#### UpSet Plot Analysis Functions  
- **`create_publication_upset_plot()`** (lines 2755-2863)
  - **Output**: `upset_plot.png|pdf`, `upset_plot_filtered.png|pdf`
  - **Features**: Advanced filtering (min_subset_size, max_degree), publication-ready styling

- **`analyze_upset_statistics()`** (lines 2864-2901)
  - **Purpose**: Comprehensive intersection statistics with degree distribution analysis

- **`create_upset_summary_plot()`** (lines 2902-2968)
  - **Output**: `upset_plot_summary.png|pdf`
  - **Features**: Species coverage and degree distribution visualizations

### 4-ultimate-EnrichmentAnalysis4SpecificOGs.qmd

#### GO Enrichment Functions
- **`enrich_orthogroups()`** (R function)
  - **Output**: GO enrichment barplots and dotplots for each analysis category
  - **Files**: 
    - `./Notebooks/Figures/GO_Enrichment/*/GO_enrichment_barplot.png`
    - `./Notebooks/Figures/GO_Enrichment/*/GO_enrichment_dotplot.png`
  - **Categories**: Increased, decreased, and species-specific gene families
  - **Ontologies**: BP (Biological Process), CC (Cellular Component), MF (Molecular Function)

#### KEGG Pathway Functions  
- **`enrich_orthogroups_with_enricher()`** (R function)
  - **Output**: KEGG pathway enrichment visualizations
  - **Files**: 
    - `./Notebooks/Figures/KO_Enrichment/*/KO_enrichment_barplot.png`
    - `./Notebooks/Figures/KO_Enrichment/*/KO_enrichment_dotplot.png`
  - **Categories**: Increased, decreased, and species-specific gene families

### 4-ultimate-EnrichmentAnalysis-preparing.qmd

#### Test Enrichment Plots
- **Test GO Analysis** (R code cells)
  - **Output**: Validation plots for GO enrichment setup
  - **Files**: 
    - `./Outputs/GO_Enrichment/GO_enrichment_test_barplot.png|pdf`
    - `./Outputs/GO_Enrichment/GO_enrichment_test_dotplot.png|pdf`

### 5-ultimate-PastML.qmd

#### Enhanced Phylogenetic Functions
- **Enhanced versions** of visualization functions from notebook 3
  - **Output**: Enhanced species trees with auto-generated file paths
  - **Files**: `./Notebooks/Figures/*.png|pdf` (auto-generated names)
  - **Features**: Automatic output path generation, dual format saving
  - **Enhanced gene trees**: Automatic orthogroup ID extraction for file naming
  - **Enhanced duplication plots**: Custom paths and enhanced styling

### Figure Output Directory Structure

```
Notebooks/Figures/
├── GO_Enrichment/          # GO enrichment analysis plots
│   ├── Increased_*_Plots/  # Increased gene family GO plots  
│   ├── decreased_*_Plots/  # Decreased gene family GO plots
│   └── unique_nonzero_*_Plots/ # Species-specific GO plots
├── KO_Enrichment/          # KEGG pathway enrichment plots
│   ├── increased_KO_Plots/ # Increased gene family KO plots
│   ├── decreased_KO_Plots/ # Decreased gene family KO plots  
│   └── unique_nonzero_KO_Plots/ # Species-specific KO plots
├── OG_bars/               # Individual orthogroup bar charts (300+ files)
│   └── OG*_OG_genecount_barplot.png|pdf
├── PastMLTrees/           # PastML ancestral reconstruction results  
│   └── OG*_map.html       # Interactive HTML visualizations
├── SpeciesTree_*.png|pdf  # Publication-quality species phylogenetic trees
├── SpeciesTree_*_ggtree.png|pdf # Advanced colorful ggtree visualizations
├── SpeciesTree_comparison_*.png|pdf # Multi-scheme tree comparisons
├── SpeciesTree_circular_*.png|pdf # Circular tree layouts
├── SpeciesTree_fan_*.png|pdf # Fan-shaped tree layouts
├── *gene_tree*.png|pdf    # Gene trees for specific orthogroups
├── heatmap_*.png|pdf      # Advanced comparative heatmaps with transforms
├── heatmap_*_print.png    # High-resolution print versions (600 DPI)
├── duplications_*.png|pdf # Gene duplication analysis plots
├── upset_plot*.png|pdf    # Advanced set intersection plots
├── upset_plot_summary.png|pdf # UpSet analysis summary plots
├── orthologs-ratio*.png|pdf # Gene/orthogroup ratio analysis
├── Orthogroups_SpeciesOverlaps*.png|pdf # Species overlap matrices
└── Orthogroups_SpeciesOverlaps*_print.png # High-res overlap matrices
```

### Function Usage Examples

```python
# Generate publication-quality species tree with enhanced styling
plot_phylogenetic_tree_with_lengths(
    newick_file_path="./tree.txt", 
    save_path="./Notebooks/Figures/my_tree",
    publication_ready=True,
    dpi=300
)

# Create advanced colorful ggtree visualization (R)
plot_phylogenetic_tree_ggtree(
    newick_file_path="./tree.txt",
    layout="rectangular",
    color_scheme="clade_auto",
    color_palette="Set1",
    save_path="./Notebooks/Figures/my_tree_ggtree"
)

# Create publication-quality heatmap with advanced options
create_publication_heatmap(
    data=orthogroup_data,
    title="Gene Counts Heatmap",
    save_path="./Notebooks/Figures/my_heatmap",
    publication_ready=True,
    journal_style="nature",
    data_transform="log10"
)

# Generate UpSet plot for orthogroup intersections
create_publication_upset_plot(
    file_path="./one-hot-orthogroups.tsv",
    save_path="./Notebooks/Figures/upset_analysis",
    min_subset_size=50,
    publication_ready=True
)
```

## Troubleshooting

### Common Issues

1. **Docker container fails to start**:
   - Ensure Docker has sufficient memory allocation (≥12GB)
   - Check port availability (38888, 30000)

2. **OrthoFinder analysis fails**:
   - Verify protein sequence file quality in `DB/GENOME_Comparison/FAA/`
   - Check for corrupted or empty FAA files

3. **Jupyter kernels not available**:
   - Restart the container: `docker-compose down && ./init.sh`
   - Manually install kernels following notebook 0 instructions

4. **Memory issues during analysis**:
   - Reduce batch sizes in analysis notebooks
   - Consider using GPU acceleration with `./init.sh GPU`

### Performance Optimization

- **Use GPU acceleration** when available for DIAMOND searches
- **Process datasets in batches** for memory-intensive operations  
- **Monitor container resources** using `docker stats`

## Citation and Acknowledgments

This pipeline integrates multiple established bioinformatics tools:
- OrthoFinder for orthology inference
- PastML for ancestral state reconstruction  
- MCScanX for synteny analysis
- clusterProfiler for functional enrichment

Please cite the relevant tools and methods when using this pipeline in your research.

## Contact and Support

For technical issues or questions about the analysis pipeline, please refer to:
- Individual notebook documentation for specific steps
- Tool-specific documentation for parameter optimization
- Container logs for debugging environment issues

---

**Last Updated**: 2025-08-28  
**Analysis Version**: Ultimate Pipeline v1.0  
**Container Version**: huang/comparative-genomics:20250828
