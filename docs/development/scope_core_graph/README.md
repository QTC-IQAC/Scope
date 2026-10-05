# SCOPE Core Call Graph

This directory contains an interactive call map of the Python modules under `core/src/scope`. It is intended as an architectural and debugging aid for SCOPE developers; it is not part of SCOPE.

## Contents

- `generate.ipynb`: parses the core source and regenerates the graph.
- `scope_core_graph.html`: interactive developer-facing graph.
- `scope_core_graph.json`: definitions, resolved call relations, call sites, and unresolved calls used by the HTML.

The generated data records the Git commit, generation time, and whether the working tree contained local changes.

## Using The Graph

Open `scope_core_graph.html` in a browser.

- **Module overview** shows calls between core modules. Double-click a module to inspect its contents.
- **Functions in selected module** shows the functions and methods defined in one module.
- **Entire core** displays every definition and can be visually dense.
- **Global function search** focuses the graph on one qualified function name, independently of the selected module.
- **Depth 1** shows the selected function and its immediate callers and callees.
- **Depth 2** additionally shows the immediate connections of those neighboring functions. This is the default because it generally provides useful workflow context without expanding to the entire package.
- **Depth 3** provides broader context but may become cluttered.

Clicking a function displays its source location, resolved callers and callees, and calls that could not be resolved safely.

## Regenerating The Snapshot

From the SCOPE repository root, run:

```bash
conda run -n scope_dev jupyter nbconvert \
  --to notebook \
  --execute \
  --inplace \
  docs/development/scope_core_graph/generate.ipynb
```

The notebook locates the repository from its current working directory and writes the HTML and JSON files back to this directory.

The graph is a developer snapshot and does not need to be regenerated after every commit. Refresh it for releases or after substantial changes to package architecture, object relationships, or call paths.

## Requirements And Limitations

Generating the graph requires Python, Jupyter, and IPython, but adds no dependency to the SCOPE package. The parser itself uses only the Python standard library. The HTML currently loads `vis-network` from a public CDN, so internet access is required when opening the graph.

The map is based on static analysis. Direct functions, imported functions, class methods, `self`, `cls`, and many `super` calls can be resolved, but Python features such as dynamic dispatch, runtime imports, and callbacks cannot always be safely assigned to one implementation. Such calls are reported as unresolved, but that does not necessarily indicate an error.
