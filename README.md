[![DOI](https://img.shields.io/badge/DOI-10.1234%2Fabcd1234-blue)]() ![SLURM Compatible](https://img.shields.io/badge/HPC%20scheduler-SLURM-green)

# SCOPE

  This repository contains the source code of SCOPE, a Python package designed to handle computational chemistry workflows.  
  It prepares, submits, and analyses Gaussian and Quantum Espresso computations of individual molecules or unit cells.  
  
  SCOPE has three main modules.
  1) Chemical Species:        Enables the creation of chemistry-related classes (e.g. Molecules, Ligands, Atoms, Bonds, Cells) that can be navigated and interacted with.
  2) Computational Workflow:  Enables setting any type of computational task that chemical species must go through.
  3) Environment:             Enables the management of files within and between computers, and the submission of jobs to computational clusters. 

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
  - CLI tools

---

# Installation

  ```bash
  
  # create and activate conda environment and install pip
  conda create --name scope 
  conda activate scope
  conda install pip
  
  # clone repository and enter
  git clone https://github.com/QTC-IQAC/Scope.git
  cd Scope
  
  # install with pip 
  pip install -e .
  ```

  ## Dependencies

  - numpy
  - networkx
  - rdkit
  - platformdirs
  - [cell2mol](https://github.com/lcmd-epfl/cell2mol.git)

---

  # Usage

Comprehensive Tutorials are available in their own github repository [here](https://github.com/QTC-IQAC/Scope_Tutorials.git):

---

  # License

MIT

---

  # Acknowledgements
- Manel Serrano and Raul Santiago (IQAC-CSIC) for their coding contribution, and help at setting the repository
- The Spanish Ministerio de Ciencia, Innovación y Universidades for funding (Project PID2022-138265NA-I00)
- The EuroHPC Development Access Call (Project: EHPC-DEV-2024D11-031)
- The Centre de Supercomputació de Catalunya (CSUC) for Computational Resources

<p align="center">
  <img src="mciu_logo.png" alt="Logo" width="400">
</p>

