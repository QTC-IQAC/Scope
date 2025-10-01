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

- ✅ Configure environments with paths, software, and queues  
- ✅ Run single or multiple spin-crossover jobs  
- ✅ CLI tools for common workflows  
- ✅ Integration with [cell2mol](https://github.com/lcmd-epfl/cell2mol)  

---

# Installation

```bash
conda create --name scope python=3.10  ## Cosymlib does not accept higher versions
conda activate scope
conda install pip

git clone https://github.com/QTC-IQAC/Scope.git
cd scope
pip install -e .
```
  ## Cosymlib does not accept higher versions
## Dependencies

- numpy
- scikit learn
- rdkit
- platformdirs
- cell2mol
- ipykernel ?

---

# Usage

The following tutorials are available at the github repository (https://github.com/QTC-IQAC/Scope.git):

- Configuration of an Environment
- 

# License

MIT

# Acknowledgements
- Raul Santiago (IQAC-CSIC), for his coding contribution and help at setting the repository
- The Spanish Ministerio de Ciencia, Innovación y Universidades for funding (Project PID2022-138265NA-I00)