[![DOI](https://img.shields.io/badge/10.26434/chemrxiv.15001415/v1-blue)](https://doi.org/10.26434/chemrxiv.15001415/v1) ![SLURM Compatible](https://img.shields.io/badge/HPC%20scheduler-SLURM-green)

# SCOPE: A Chemically-Aware Workflow-Automation Software for Molecules and Molecular Crystals

  This repository contains the tutorials associated with SCOPE, a Python package designed to orchestrate computational chemistry workflows. SCOPE prepares, submits, and analyses quantum chemistry computations of individual molecules or molecule-based crystals, and organizes the results to simplify analysis.  
  
  SCOPE has four core pillars.
  1) Chemical Species:        Dedicated to encode the chemistry of the systems of interest. It can be expanded through add-ons
  2) Computational Workflow:  Dedicated to define dynamic computational workflows.
  3) States:                  Dedicated to the analysis of results
  4) Environment:             Dedicated to file management and job execution in HPC clusters.

---

# Documentation

  - The associated manuscript is under preparation. Preprint available [here](https://doi.org/10.26434/chemrxiv.15001415/v1)
  - Repository and source [here](https://github.com/QTC-IQAC/Scope)
  - Tutorials are available [here](https://github.com/QTC-IQAC/Scope_Tutorials)

## Developer Docs

  For contributors and editor-based AI tools, the repository also includes:
  - [Architecture Notes](docs/architecture.md)
  - [Concepts And Terminology](docs/concepts.md)
  - [Agent Guidance](AGENTS.md)

---

# Features

  - Configure computational environments for your projects, storing paths, software, and queues  
  - Run quantum chemistry workflows for molecules or periodic structures, using Quantum Espresso or Gaussian16
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
  pip install git+https://github.com/lcmd-epfl/cell2mol.git
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

  ### Optional Quantum Espresso pseudopotential Libraries:
  - `scope-qc` ships only the `Vanderbilt_USPP` library
  - The larger `Efficiency` and `Precision` libraries are distributed separately through the GitHub releases page: <https://github.com/QTC-IQAC/Scope/releases>
  - if you want to use them, download the corresponding release asset and extract the `Efficiency/` or `Precision/` folder into `scope/software/quantum_espresso/PP_Libraries/` inside your installed `scope-qc` package

  ## Dependencies

  - python 3.12
  - numpy
  - networkx < 3.3
  - rdkit
  - platformdirs
  - [cell2mol](https://github.com/lcmd-epfl/cell2mol.git)

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

  For Quantum Espresso runs, `pp_library = vanderbilt` works out of the box.
  If you set `pp_library = efficiency` or `pp_library = precision`, make sure the corresponding library has first been downloaded from the GitHub releases page and extracted into `scope/software/quantum_espresso/PP_Libraries/`.

  ## Interactive
  Comprehensive Tutorials are available in their own github repository [here](https://github.com/QTC-IQAC/Scope_Tutorials.git):

---

  # License

This project is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International license (`CC BY-NC-ND 4.0`). See the [LICENSE](LICENSE) file for the full text.

---

  # Acknowledgements
- The Spanish Ministerio de Ciencia, Innovación y Universidades for funding (Project PID2022-138265NA-I00)
- The EuroHPC Development Access Call (Project: EHPC-DEV-2024D11-031)
- The Centre de Supercomputació de Catalunya (CSUC) for Computational Resources

<p align="center">
  <img src="mciu_logo.png" alt="Logo" width="400">
</p>
