# scope-qc

`scope-qc` is the core SCOPE package: A Chemically-Aware Workflow-Automation Software for Molecules and Molecular Crystals

It provides the main object model and workflow engine used across the SCOPE project.
Optional add-ons such as `scope-azo` and `scope-sco` build on top of this core package.

## Documentation

- Repository and source code: <https://github.com/QTC-IQAC/Scope>
- Preprint: <https://doi.org/10.26434/chemrxiv.15001415/v1>
- Tutorials: <https://github.com/QTC-IQAC/Scope_Tutorials>

## Installation

create and activate conda environment and install pip:
```bash
conda create --name scope python=3.12 
conda activate scope
conda install pip

# `cell2mol` is an external dependency and must currently be installed separately from its source repository.
pip install git+https://github.com/lcmd-epfl/cell2mol.git
```

### Option 1: from PyPI 
```bash
pip install scope-qc    # Core Modules
```

Quantum Espresso pseudopotentials:
- `scope-qc` ships only the `Vanderbilt_USPP` library
- `Efficiency` and `Precision` are distributed separately through the GitHub releases page: <https://github.com/QTC-IQAC/Scope/releases>
- To use either of them, download the corresponding release asset and extract the `Efficiency/` or `Precision/` folder into `scope/software/quantum_espresso/PP_Libraries/` inside your installed `scope-qc` package

### Option 2: from repository 
```bash
git clone https://github.com/QTC-IQAC/Scope.git
cd Scope
pip install -e core
```

### Option 3: from TestPyPi: 
```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ scope-qc
```

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

External prerequisite:

- `cell2mol`, installed separately from <https://github.com/lcmd-epfl/cell2mol.git>

Python requirement:

- Python 3.12

## Command Line Interface

The core package provides the `scope` command with subcommands for configuring SCOPE environments, creating systems, and running tasks.

For help:

```bash
scope -h
```

## Usage

SCOPE is typically used through the command line to:
- configure environments: `scope config`
- create systems:         `scope create_single` and `scope create_many`
- execute workflows:      `scope run`

For Quantum Espresso jobs, `pp_library = vanderbilt` is available by default.
If you select `efficiency` or `precision`, install that library first from the GitHub releases page above.

`System` objects are saved in binary files that can then be inspected later in notebooks or other interactive sessions.

## License

See the repository-level [LICENSE](../LICENSE) file for licensing information.
