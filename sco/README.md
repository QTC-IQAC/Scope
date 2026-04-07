# scope-sco

`scope-sco` is the spin-crossover add-on for SCOPE. It extends the core `scope` package with utilities for building, running, and analysing spin-crossover workflows.

This add-on is intended for users who already work with SCOPE's core chemistry and workflow abstractions and need spin-crossover-specific functionality on top of them.

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

### Option 1: from PyPI
```bash
pip install scope-sco
```

### Option 2: from repository
```bash
git clone https://github.com/QTC-IQAC/Scope.git
cd Scope
pip install -e sco
```

### Option 3: from TestPyPi
```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ scope-sco
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
