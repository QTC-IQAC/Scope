[![DOI](https://img.shields.io/badge/DOI-10.1234%2Fabcd1234-blue)]() ![SLURM Compatible](https://img.shields.io/badge/HPC%20scheduler-SLURM-green)

# SCOPE: A Chemically-Aware Workflow-Automation Software for Molecules and Molecular Crystals

  This repository contains the source code of SCOPE, a Python package designed to handle computational chemistry workflows. SCOPE prepares, submits, and analyses Gaussian and Quantum Espresso computations of individual molecules or molecule-based crystals.
  
  SCOPE has three core functionalities.
  1) Chemical Species:        Dedicated to the definition of chemistry-related classes (e.g. Molecules, Ligands, Atoms, Bonds, Cells) that can be navigated and interacted with.
  2) Computational Workflow:  Dedicated to the orchestration of dynamic computational workflows.
  3) States:                  Dedicated to the analysis of results
  4) Environment:             Dedicated to file management and job execution in HPC clusters.

---

# Documentation

  - The associated manuscript is under preparation. Partial documentation is available at: ???
  - Tutorials are available [here](https://github.com/QTC-IQAC/Scope_Tutorials)

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

  ```bash
  
  # create and activate conda environment and install pip
  conda create --name scope python=3.12 
  conda activate scope
  conda install pip
  
  # clone repository and enter
  git clone https://github.com/QTC-IQAC/Scope.git
  cd Scope
  
  # install core with pip 
  pip install -e core

  # install add-ons with pip if you want them
  conda install openbabel -c conda-forge             ## only needed for the azo add-on
  pip install -e azo
  pip install -e sco
  ```

  Installation when deployed to pip

  ```bash

  # create and activate conda environment and install pip
  conda create --name scope python=3.12 
  conda activate scope
  conda install pip

  # install core with pip
  pip install scope

  # install add-ons with pip if you want them
  conda install openbabel -c conda-forge             ## only needed for the azo add-on
  pip install scope-azo
  pip install scope-sco
  ```

  ## Dependencies

  - python 3.12
  - numpy
  - networkx < 3.3
  - rdkit
  - platformdirs
  - [cell2mol](https://github.com/lcmd-epfl/cell2mol.git)

---

  # Usage

  Ideally, CLI commands are used to configure and run SCOPE. The results are stored in binary files containing System-class objects
  These objects are conceived to be inspected interactively in notebooks.   

  ## Command Line Interface: 

  SCOPE provides a single top-level command, "scope", with multiple subcommands.
   - config         Configure the SCOPE environment
   - create_many    Create many systems from xyz data
   - create_single  Create one system from a xyz source
   - run            Run a SCOPE task for a given system
   - set_path       Set the current directory as the system main path

  All subcommands have a dedicated --help with the intended use. For instance:
  ```bash
  scope config -h 
  ```

  ## Interactive
  Comprehensive Tutorials are available in their own github repository [here](https://github.com/QTC-IQAC/Scope_Tutorials.git):

---

  # License

MIT

---

  # Acknowledgements
- The Spanish Ministerio de Ciencia, Innovación y Universidades for funding (Project PID2022-138265NA-I00)
- The EuroHPC Development Access Call (Project: EHPC-DEV-2024D11-031)
- The Centre de Supercomputació de Catalunya (CSUC) for Computational Resources

<p align="center">
  <img src="mciu_logo.png" alt="Logo" width="400">
</p>
