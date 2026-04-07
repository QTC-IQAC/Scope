# scope

`scope` is the core SCOPE package: a chemically aware workflow-automation toolkit for molecules and molecular crystals.

It provides the main object model and workflow engine used across the SCOPE project:

- chemistry-aware classes for molecules, ligands, atoms, bonds, and cells
- workflow orchestration for preparing, submitting, and tracking computations
- parsing and registration of Gaussian and Quantum Espresso results
- environment, path, and queue management for HPC execution
- a CLI for configuring projects, creating systems, and running tasks

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
```

`cell2mol` is an external dependency and must currently be installed separately from its source repository.
```bash
pip install git+https://github.com/lcmd-epfl/cell2mol.git
```

### Install from PyPI:

```bash
pip install scope-qc    # Core Modules
```

### For a local editable installation from this repository:

```bash
pip install -e core
```

### Install from TestPyPi:
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

Shell completion is available through `argcomplete`. After installation, enable it in your shell with:

```bash
activate-global-python-argcomplete
```

## Usage

SCOPE is typically used through the command line to configure environments, create systems, and execute computational workflows. Serialized `System` objects can then be inspected later in notebooks or other interactive sessions.

## License

See the repository-level [LICENSE](../LICENSE) file for licensing information.
