# SCOPE Benchmarks

This folder contains three reproducible benchmarks covering molecular import, ligand filtering and workflow execution, and molecular overlap. The datasets are not distributed with SCOPE; the notebooks use the supplied relative-path selections to identify the same inputs after the datasets have been downloaded locally.

## Available benchmarks

| Benchmark | Purpose | Main Files |
| --- | --- | --- |
| [1 – Import](1-Import/1-Import.ipynb) | Imports and validates cell2mol `Cell` objects, FORMED XYZ structures, and GEOM-QM9 RDKit molecules. | Three `selected_*.txt` files and three `*_import_results.csv` reports. |
| [2 – Filtering and Execution](2-Filtering_and_Execution/2-Filtering_and_Execution.ipynb) | Finds unique ligands in cell2mol Cells, collects them in one SCOPE `System`, and prepares and inspects Quantum ESPRESSO calculations. | `selected_cell_paths.txt`, `unique.npy`, `task1.scope`, `task2.scope`, and `unique_finished.npy`. |
| [3 – Molecule Overlap](3-Molecule_Overlap/3-Molecule_Comparison.ipynb) | Compares the SCOPE overlap routine with RDKit constructions with and without bond information, recording timing and structural metrics. | `selected_cell_paths.txt`, `molecule_comparison_results.csv`, and optional XYZ files in `overlap_comparison_tests/`. |

Shared utilities that are not part of the SCOPE package are kept in [benchmark_functions.py](benchmark_functions.py).

## 1. Install SCOPE

Python 3.12 is required. Create an environment and follow the main [SCOPE installation instructions](../README.md#installation). A minimal installation for the notebooks is:

```bash
conda create --name scope_benchmarks python=3.12
conda activate scope_benchmarks
conda install pip
pip install "cell2mol @ git+https://github.com/lcmd-epfl/cell2mol.git@55649ba5f444895846a5d049adcb056c01bb0421"
pip install scope-qc
```

The `scope-azo` and `scope-sco` add-ons are not required. To benchmark a local source checkout instead of the published package, replace the final command with `pip install -e core` from the repository root.

Benchmark 2 additionally requires access to Quantum ESPRESSO 7.0 and a configured SCOPE execution environment only if the calculations themselves are to be reproduced. Its supplied `unique_finished.npy` can be inspected locally without rerunning Quantum ESPRESSO.

## 2. Download and arrange the datasets

The benchmarks use these external datasets:

- [cell2mol CSD dataset](https://archive.materialscloud.org/records/cdt4h-yvv48), used by all three benchmarks.
- [FORMED dataset](https://archive.materialscloud.org/records/20dc5-mvx03), used by Benchmark 1.
- [GEOM QM9 dataset](https://github.com/learningmatter-mit/geom), used by Benchmark 1. Extract its Python-specific `rdkit_folder.tar.gz` archive so that the RDKit pickle files and QM9 summary are available.

Arrange the extracted data under one common parent folder:

```text
<datasets_folder>/
├── cell2mol/
│   ├── 1-Iron/
│   ├── 2-Manganese/
│   ├── 3-Ruthenium/
│   ├── 4-Rhenium/
│   ├── 5-Chromium/
│   ├── 6-Cobalt/
│   ├── 7-Nickel/
│   └── 8-Copper/
├── formed/
│   └── *.xyz
└── rdkit_folder/
    ├── summary_qm9.json
    └── qm9/
        └── *.pickle
```

The GEOM summary may originally be named `qm9_summary.json`; Benchmark 1 accepts either `summary_qm9.json` or `qm9_summary.json`.

## 3. Configure and run a benchmark

Open the desired notebook with the SCOPE environment and execute it with its own folder as the working directory. In the Part 0 configuration cell, change only `datasets_folder` to the common parent folder shown above—not to an individual dataset folder.

For example:

```bash
conda activate scope_benchmarks
cd benchmarks/1-Import
jupyter notebook 1-Import.ipynb
```

Use the equivalent notebook directory for Benchmarks 2 and 3. Each notebook derives `benchmark_folder` from its current working directory and imports the shared helper module from the parent `benchmarks/` folder.

The notebooks store relative dataset paths, so the supplied selections remain valid when `<datasets_folder>` changes. Before reproducing the retained Benchmark 1 results, set all three selection controls to `False`:

```python
overwrite_cell_selection = False
overwrite_xyz_selection  = False
overwrite_mol_selection  = False
```

With those settings, Benchmark 1 reuses both an existing selection and its corresponding CSV report. Set an overwrite option to `True` only when intentionally generating a new random selection and replacing that part's results. Benchmark 2 similarly reuses its existing Cell selection when `overwrite_cell_selection = False`, although its ligand comparison and `unique.npy` creation are performed again.

Benchmark 3 always uses its supplied 25-Cell selection. Running its main comparison recreates `molecule_comparison_results.csv`; the selected inputs remain unchanged. Elapsed times depend on the machine and current load, whereas the RMSD and atom-displacement results are the structural comparison metrics intended for reproduction.

If a dataset required by Benchmark 1 or 2 is absent, its corresponding section is skipped. Benchmark 3 requires the cell2mol dataset and stops with a clear error if it cannot be found.

## Benchmark 2 execution files

After Benchmark 2 creates `unique.npy`, copy it together with `task1.scope` and `task2.scope` to a configured SCOPE project on the execution machine. The notebook lists the complete submission procedure. Run the two tasks with Quantum ESPRESSO 7.0, download the resulting `System`, and save it beside the notebook as `unique_finished.npy`. The final cells load that file and report finished and unfinished jobs.

For reproducible timing reports, record the processor or compute node, SCOPE version, RDKit version, and Quantum ESPRESSO version alongside the benchmark results.
