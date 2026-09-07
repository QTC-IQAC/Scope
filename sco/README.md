# scope-sco

`scope-sco` is the spin-crossover add-on for SCOPE. It extends the core `scope` package with utilities for building, running, and analysing spin-crossover workflows.

This add-on is intended for users who already work with SCOPE's core chemistry and workflow abstractions and need spin-crossover-specific functionality on top of them.

## Documentation

- Repository and source code: <https://github.com/QTC-IQAC/Scope>
- Preprint: <https://doi.org/10.26434/chemrxiv.15001415/v2>
- Tutorials: <https://github.com/QTC-IQAC/Scope_Tutorials>

## Installation

Installing `scope-sco` automatically installs `scope-qc`, but the external `cell2mol` prerequisite must be installed separately.

```bash
# create and activate conda environment and install pip
conda create --name scope python=3.12
conda activate scope
conda install pip

pip install "cell2mol @ git+https://github.com/lcmd-epfl/cell2mol.git@55649ba5f444895846a5d049adcb056c01bb0421"
```

### Option 1 (preferred): from PyPI
```bash
pip install scope-sco
```

### Option 2 (alternative): from repository
```bash
git clone https://github.com/QTC-IQAC/Scope.git
cd Scope
pip install -e core  # optional, otherwise the core package (scope-qc) will be installed from pip 
pip install -e sco
```

## Command Line Interface

The package provides the `scope_sco` command. For help:

```bash
scope_sco -h
```

## Usage

`scope-sco` complements the main `scope` package and is meant to be used together with the core SCOPE workflow and data model.

## License

See the repository-level [LICENSE](../LICENSE) file for licensing information.
