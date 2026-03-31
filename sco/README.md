# scope-sco

`scope-sco` is the spin-crossover add-on for SCOPE. It extends the core `scope` package with utilities for building, running, and analysing spin-crossover workflows.

This add-on is intended for users who already work with SCOPE's core chemistry and workflow abstractions and need spin-crossover-specific functionality on top of them.

## Documentation

- Repository and source code: <https://github.com/QTC-IQAC/Scope>
- Preprint: <https://doi.org/10.26434/chemrxiv.15001415/v1>
- Tutorials: <https://github.com/QTC-IQAC/Scope_Tutorials>

## Installation

Python 3.12 is required.

Install the core package first, then install the add-on from PyPI:

```bash
pip install scope-qc
pip install scope-sco
```

For a local editable installation from this repository:

```bash
conda create --name scope python=3.12
conda activate scope
pip install -e core
pip install -e sco
```

## Command Line Interface

The package provides the `scope_sco` command. For help:

```bash
scope_sco -h
```

Shell completion is available through `argcomplete`. After installation, enable it in your shell with:

```bash
activate-global-python-argcomplete
```

## Usage

`scope-sco` complements the main `scope` package and is meant to be used together with the core SCOPE workflow and data model.

## License

See the repository-level [LICENSE](../LICENSE) file for licensing information.
