# scope-azo

`scope-azo` is the azo add-on for SCOPE. It extends the core `scope` package with tools for building and working with azo systems.

## Installation

Install the core package first, then install the add-on:

```bash
pip install scope
conda install openbabel -c conda-forge
pip install scope-azo
```

For a local editable installation from this repository:

```bash
pip install -e core
conda install openbabel -c conda-forge
pip install -e azo
```

## Command Line Interface

The package provides the `scope_azo` command. For help:

```bash
scope_azo -h
```

## License

See the repository-level [LICENSE](../LICENSE) file for licensing information.
