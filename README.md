# About SCOPE

This repository contains the SCOPE package developed by Sergi Vela and collaborators, at the IQAC(CSIC)
SCOPE is a Python package designed to handle computational chemistry workflows. It prepares, submits, and analyses Gaussian and Quantum Espresso computations of individual molecules or unit cells. 

SCOPE has three main modules.
1) Chemical Species:        Enables the creation of chemistry-related classes (e.g. Molecules, Ligands, Atoms, Bonds, Cells) that can be navigated and interacted with.
2) Computational Workflow:  Enables setting any type of computational task that chemical species must go through.
3) Environment:             Enables the management of files within and between computers, and the submission of jobs to computational clusters. 

# Documentation

The associated manuscript is under preparation. Partial documentation is available at: ???

---

# Features

- Configure computational environments for your projects, storing paths, software, and queues  
- Run quantum chemistry workflows for one or multiple molecules or periodic structures, using Quantum Espresso or Gaussian16
- Parse and analyse results of computations, and connect the data with SCOPE's molecule- and cell-class objects
- CLI tools for common workflows  
- Integration with [cell2mol](https://github.com/lcmd-epfl/cell2mol)  

---

# Installation

  ## Safest way
```bash
# Create Environment
conda create -n scope python=3.10   ## Cosymlib does not accept higher versions 
conda activate scope

# Install dependencies
pip install --upgrade pip setuptools wheel
pip install "numpy<=1.22.4"
pip install scipy rdkit platformdirs ipykernel plotly nbformat

# Download Code and Install
git clone --branch dev https://github.com/QTC-IQAC/Scope.git
cd Scope
pip install -e .
```

  ## Fastest way
```bash
conda create --name scope python=3.10  ## Cosymlib does not accept higher versions
conda activate scope
conda install pip
git clone --branch dev https://github.com/QTC-IQAC/Scope.git
cd Scope
pip install -e .
```

  ## Dependencies

- numpy
- scikit learn
- rdkit
- platformdirs
- cell2mol
- ipykernel ?

---

  # Usage

Tutorials are available at the github repository (https://github.com/QTC-IQAC/Scope.git):
Alternatively, you can download and install those tutorials together with the main package, by doing:

```bash
pip install -e .[tutorials]
```

rather than the bare "pip install -e ." that is suggested above: 

  # License

MIT

  # Acknowledgements
- Manel Serrano and Raul Santiago (IQAC-CSIC), for their coding contribution, and help at setting the repository
- The Spanish Ministerio de Ciencia, Innovación y Universidades for funding (Project PID2022-138265NA-I00)
