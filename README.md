[![DOI](https://img.shields.io/badge/10.26434/chemrxiv.15001415/v2-blue)](https://doi.org/10.26434/chemrxiv.15001415/v2) ![SLURM Compatible](https://img.shields.io/badge/HPC%20scheduler-SLURM-green)

# SCOPE: A Chemically-Aware Workflow-Automation Software for Molecules and Molecular Crystals

  This repository contains the SCOPE core package, its optional add-ons, benchmarks, and development material. SCOPE is a Python package designed to orchestrate computational chemistry workflows: it prepares, submits, and analyses quantum chemistry computations of individual molecules or molecule-based crystals, and organizes the results to simplify analysis.
  
  SCOPE has four core pillars.
  1) Chemical Species:        Dedicated to encode the chemistry of the systems of interest. It can be expanded through add-ons
  2) Computational Workflow:  Dedicated to define dynamic computational workflows.
  3) States:                  Dedicated to the analysis of results
  4) Environment:             Dedicated to file management and job execution in HPC clusters.

---

# Documentation

  - The associated preprint is available [here](https://doi.org/10.26434/chemrxiv.15001415/v2)
  - Repository and source [here](https://github.com/QTC-IQAC/Scope)
  - Tutorials are available [here](https://github.com/QTC-IQAC/Scope_Tutorials)

## Developer Docs

  For contributors and editor-based AI tools, the repository also includes:
  - [Architecture Notes](docs/architecture.md)
  - [Concepts And Terminology](docs/concepts.md)
  - [Coding Style](docs/coding_style.md)
  - [Bug-fixing Procedure](docs/bug_fixing.md)
  - [Benchmarks](benchmarks/README.md)
  - [Agent Guidance](AGENTS.md)

---

# Features

  - Configure computational environments for your projects, storing paths, software, and queues  
  - Run quantum chemistry workflows for molecules or periodic structures, using Quantum ESPRESSO (v7.0) or Gaussian16
  - Parse and analyse results of computations, and connect the data with SCOPE's molecule- and cell-class objects
  - Integration with [cell2mol](https://github.com/lcmd-epfl/cell2mol)  
  - Optional `sco` and `azo` add-ons extend SCOPE with additional capabilities
  - CLI tools

---

# Installation

Python 3.12 is a strict requirement for SCOPE and its add-ons.

  ```bash
  # 1-create and activate conda environment and install pip
  conda create --name scope python=3.12 
  conda activate scope
  conda install pip

  # 2-install external prerequisite
  pip install "cell2mol @ git+https://github.com/lcmd-epfl/cell2mol.git@55649ba5f444895846a5d049adcb056c01bb0421"
  ``` 
  
  ### Option 1 (preferred): from pip
  ```bash
  # 1-install core
  pip install scope-qc

  # 2-install add-ons (optional)
  conda install openbabel -c conda-forge             # only needed for the azo add-on
  pip install scope-azo
  pip install scope-sco
  ```


  ### Option 2 (alternative): from repository
  ```bash
  # 0-Download repo
  git clone https://github.com/QTC-IQAC/Scope.git
  cd Scope

  # 1-install core
  pip install -e core

  # 2-install add-ons (optional) 
  conda install openbabel -c conda-forge             # only needed for the azo add-on
  pip install -e azo
  pip install -e sco
  ```

  ### Optional Quantum ESPRESSO pseudopotential Libraries:

  - `scope-qc` ships only the `Vanderbilt_USPP` library
  - The larger `Efficiency` and `Precision` libraries are distributed separately through the GitHub releases page: <https://github.com/QTC-IQAC/Scope/releases>
  - if you want to use them, download the corresponding release asset and extract the `Efficiency/` or `Precision/` folder into `scope/software/quantum_espresso/PP_Libraries/` inside your installed `scope-qc` package

  ## Dependencies

  The core package depends on:

  - `numpy`
  - `networkx < 3.3`
  - `scipy`
  - `rdkit`
  - `ipykernel`
  - `plotly`
  - `nbformat`
  - `jupyter`
  - `matplotlib`
  - `platformdirs`
  - `PyYAML >= 6.0.1`

  External prerequisite:

  - `cell2mol`, installed separately from <https://github.com/lcmd-epfl/cell2mol.git>

  Python requirement:

  - Python 3.12

---

  # Usage

  Ideally, CLI commands are used to configure SCOPE, create systems, and run tasks. 
  Systems are stored in binary files, and are conceived to be inspected interactively in notebooks.   

  ## Command Line Interface: 

  SCOPE provides a single top-level command, "scope", with multiple subcommands.

  ```bash
  scope config         # Configure the SCOPE environment
  scope create_many    # Create many systems from multiple xyz data
  scope create_single  # Create one system from a single xyz source
  scope run            # Run a SCOPE task for a given system
  scope set_path       # Set the current directory as the system main path
  ```

  All subcommands have a dedicated --help with the intended use. For instance:
  ```bash
  scope config -h 
  ```

  ## Quantum ESPRESSO

  For Quantum ESPRESSO runs, `pp_library = vanderbilt` works out of the box.
  If you set `pp_library = efficiency` or `pp_library = precision`, make sure the corresponding library has first been downloaded from the GitHub releases page and extracted into `scope/software/quantum_espresso/PP_Libraries/`.

  ## Interactive
  Comprehensive tutorials are available in their own GitHub repository [here](https://github.com/QTC-IQAC/Scope_Tutorials.git):

---

  # License

  SCOPE is source-available software distributed under the PolyForm Noncommercial License 1.0.0.
  
  The license permits use, modification, and redistribution for non-commercial purposes, including use by educational institutions and public research organizations. Commercial use is not permitted under this license and requires separate authorization from the copyright holder, the Spanish National Research Council (CSIC)
  
  Users developing new functionality are encouraged to contribute to the official SCOPE repository or to distribute functionality through independent SCOPE add-ons. Modified distributions of the core should be clearly distinguished from the official SCOPE distribution.
  
  See [LICENSE](LICENSE) and [BRANDING](BRANDING.md) for details.
  
---

  # Acknowledgements
- The Spanish Ministerio de Ciencia, Innovación y Universidades for funding (Project PID2022-138265NA-I00)
- The EuroHPC Development Access Call (Project: EHPC-DEV-2024D11-031)
- The Centre de Supercomputació de Catalunya (CSUC) for Computational Resources

<p align="center">
  <img src="mciu_logo.png" alt="Logo" width="400">
</p>
