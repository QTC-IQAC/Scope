# scope-azo

`scope-azo` is the azo add-on for SCOPE. It extends the core `scope` package with tools for building, running, and analysing azo systems.

This add-on is intended for users who already work with SCOPE's core chemistry and workflow abstractions and need azo-specific functionality on top of them.

## Documentation

- Repository and source code: <https://github.com/QTC-IQAC/Scope>
- Preprint: <https://doi.org/10.26434/chemrxiv.15001415/v1>
- Tutorials: <https://github.com/QTC-IQAC/Scope_Tutorials>

## Installation

Python 3.12 is required.

Install the core package first, then install the add-on from PyPI:

```bash
pip install scope-qc
conda install openbabel -c conda-forge
pip install scope-azo
```

For a local editable installation from this repository:

```bash
conda create --name scope python=3.12
conda activate scope
pip install -e core
conda install openbabel -c conda-forge
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
