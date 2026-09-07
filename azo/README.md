# scope-azo

`scope-azo` is the azo add-on for SCOPE. It extends the core `scope` package with tools for building, running, and analysing azo-based systems.

This add-on is intended for users who already work with SCOPE's core chemistry and workflow abstractions and need azo-specific functionality on top of them.

## Documentation

- Repository and source code: <https://github.com/QTC-IQAC/Scope>
- Preprint: <https://doi.org/10.26434/chemrxiv.15001415/v2>
- Tutorials: <https://github.com/QTC-IQAC/Scope_Tutorials>

## Installation

Installing `scope-azo` automatically installs `scope-qc`, but the external `cell2mol` prerequisite must be installed separately.

```bash
# create and activate conda environment and install pip:
conda create --name scope python=3.12
conda activate scope
conda install pip

pip install "cell2mol @ git+https://github.com/lcmd-epfl/cell2mol.git@55649ba5f444895846a5d049adcb056c01bb0421"
conda install openbabel -c conda-forge
```

### Option 1 (preferred): from PyPI
```bash
pip install scope-azo
```

### Option 2 (alternative): from repository
```bash
git clone https://github.com/QTC-IQAC/Scope.git
cd Scope
pip install -e core  # optional, otherwise the core package (scope-qc) will be installed from pip 
pip install -e azo
```

## Command Line Interface

The package provides the `scope_azo` command. For help:

```bash
scope_azo -h
```

## Usage

`scope-azo` complements the main `scope` package and is meant to be used together with the core SCOPE workflow and data model.

## License

See the repository-level [LICENSE](../LICENSE) file for licensing information.
