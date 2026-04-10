# Code Repo for "Understanding polymer zip length through the lens of mechanistic modeling"

## Overview

This repository contains computational models for understanding zip length of depolymerization through a mechanistic modeling lens. The work compares two computational approaches to analytical expressions:

1. **Kinetic Monte Carlo (KMC)**
2. **Method of Moments (MoM)**
3. **Boyd (1959)**

## Publication

**Title:** Understanding polymer zip length through the lens of mechanistic modeling

**DOI:** [To be completed]

## Contents

### `kmc_src/` — Kinetic Monte Carlo Simulations

- **`spedup_cpp_files/`** — KMC source code for simple mechanisms
  - `main.cpp`, `reaction_functions.cpp/h`, `save_variables.cpp/h` — Core simulation engine in C++17
  - `Makefile` — Compilation instructions (requires g++)
  - `params.inp` — Temperature-dependent kinetic parameters (A, Ea for I/D/T/SR reactions)
  - `T_260/` — Simulations at multiple DPn values (100, 150, ..., 4000) for scaling analysis

- **`polyME_cpp_files_conv/`** — KMC source code for modeling polyME depolymerization
  - `run_conv_4_2_26/` — Simulation outputs at three conditions (T=250°C/300°C, ρ=5/25 g/L)
  - `no_sr_run_conv_4_2_26/` — Control runs without side reactions


### `mom_src/` — Method of Moments Simulations

- **`T260/run400/`** — MoM ODE solutions for comparison with KMC at DPn=400

### `paper_analysis/` — Figure Generation & Analysis

- **`plot_figures.ipynb`** — Jupyter notebook generating all publication figures (Main Figs 2–6 + SI Figs S3–S8)
- **`plot_helpers.py`** — Data loading and processing functions
  - `read_data()` — Loads KMC CSV outputs (time_based_data, zip_length, chain_weights, params)
  - `process_data()` — Calculates derived quantities (radical concentrations, Mn, zip length corrected values)
- **`presentation.mplstyle`** — Matplotlib style file for publication-quality figures
- **`fig_pdfs/`** — Output directory for generated figures

### `parameter_exploration/` — Parameter Space Analysis

- **`radical_KMC_runs/`** — Output directory for Scheme 2 radicals
- **CSV files** results across DPn and dispersity ranges
- **`rs_methods_all_dpn_og_v*.csv`** — Comprehensive parameter search results for Scheme 1
- **`eci_summary_with_plots.csv`** — Comprehensive parameter search results for Scheme 2

## Installation & Requirements

### For Analysis Only (Jupyter Notebook)

**Python:** 3.8 or higher

**Dependencies:** pandas numpy matplotlib scipy scikit-learn jupyter

### For KMC Simulations (C++ Compilation)

**Compiler:** g++ with C++17 support  

**Compile:**
```bash
cd cpp_file_src
make
```

### Hardware
- Analysis notebook: ~5 minutes on modern laptop
- KMC simulation: 4+ GB RAM, 1 core recommended

## Usage

### Reproducing Figures

**All publication figures are generated from a single notebook:**

Run cells sequentially to generate:
- **Figure 2** — Model comparison (KMC vs MoM) at DPn=400
- **Figure 3** — Scaling with initial chain length (DPn = 100, 400, 2000)
- **Figure 4** — Parameter space exploration (KMC + ECI data)
- **Figure 5** — ECI special cases across DPn values
- **Figure 6** — Best-fit kinetics at three conditions (T, ρ variations)
- **Figures S3–S8** — Supplementary analyses

### Running KMC Simulations from Scratch

1. **Modify parameters** in `kmc_src/polyME_cpp_files_conv/params.inp`
2. **Compile:**
   ```bash
   cd kmc_src/polyME_cpp_files_conv
   make
   ```
3. **Run** (example):
   ```bash
   mkdir -p run_custom && cd run_custom
   cp ../params.inp ../program .
   ./program > simulation.log
   ```
   Generates: `time_based_data.csv`, `zip_length.csv`, `chain_weights_*.csv`

4. **Analyze** with `plot_helpers.read_data()` in your script


## Data

### Included in Repository

- **KMC Outputs** — CSV time-series for key conditions (T=250/300°C, ρ=5/25 g/L)
  - `time_based_data.csv` — Time evolution: Mn, conversion, radical counts
  - `zip_length.csv` — Detected backbiting events: zip length, occurrence time
  - `chain_weights_*pct_conversion.csv` — Molecular weight snapshots at conversion milestones
  
- **Parameter Exploration** — Summary CSV files with computed ratio metrics across parameter space

### Large Data / Excluded Files

Extended simulation data (all DPn, all dispersities) available upon request. Contact authors for:
- Full `spedup_cpp_files/T_260/` tree outputs
- ECI parameter sweep raw data
- Method of Moments integration results

## License

MIT License — See LICENSE file for details

## Authors & Contact

**Author(s):** Shivani Kozarekar, Dachey Lin (PI: Linda Broadbelt)
**Affiliation:** Northwestern University, Department of Chemical and Biological Engineering
**Contact:** shivanikozarekar2026@u.northwestern.edu

## Citation

If you use this code or data, please cite:

```bibtex
@article{[author_year],
  title={Understanding polymer zip length through the lens of mechanistic modeling},
  author={[Authors]},
  journal={[Journal]},
  year={[Year]},
  doi={[DOI]}
}
```