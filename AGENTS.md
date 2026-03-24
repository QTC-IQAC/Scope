# AGENTS.md

## Project Summary
SCOPE is a Python codebase for computational chemistry workflows on molecules and molecular crystals.
It combines four responsibilities:
- chemistry-aware data structures such as molecules, ligands, cells, atoms, and bonds
- workflow orchestration for generating, submitting, and tracking computations
- parsing and registration of Gaussian and Quantum Espresso results
- environment and queue management for HPC execution

Users often interact with SCOPE through the CLI, while serialized `System` objects are a main persistent artifact for later inspection in notebooks.

The repository is split into a core package and two optional add-ons:
- `core/`: main SCOPE package
- `azo/`:  azobenzene-focused extensions
- `sco/`:  spin-crossover-focused extensions

## Core Object Model
- `System` is the top-level project container.
- `System.sources` stores chemistry sources, typically molecules or cells.
- Sources may come from native SCOPE objects, `.xyz` data, RDKit molecules, or cell2mol objects.
- Sources usually get an `"initial"` state during setup.
- `State._source` links a state to its parent source.
- `Environment` and `Queue` store runtime, filesystem, and scheduler information.
- `Data`, `Collection`, and `VNM` store parsed or derived results for later analysis.
- Workflow execution follows this hierarchy: `System -> Branch -> Workflow -> Job -> Computation`

## Key Invariants
- Keep chemistry containers internally consistent.
  `labels`, coordinates, atoms, and adjacency/connectivity data should describe the same structure.
- When attaching a source to a `System`, keep parent links and source registration in sync.
- Preserve the expected relationship between a `State` and its source.
- Many objects use lazy or incrementally populated attributes. Do not replace these patterns casually with eager recomputation.
- Path propagation across `System`, branches, workflows, jobs, and computations is part of the workflow model.
  Avoid changing it unless the task is explicitly about path handling.

## Domain Notes
- Fragmentation, reconstruction, and reference-molecule logic are scientifically sensitive.
- Transition-metal and ligand handling are central in species and cell operations.
- The `azo` and `sco` packages add domain-specific workflows and assumptions on top of the core model.
- Gaussian and Quantum Espresso parsers feed results back into SCOPE objects.
  Changes in parser output should preserve the structure expected by downstream classes.
- Scope input files and `run_task` drive workflow execution step by step.
  A run may submit pending computations, skip still-running work, or register finished jobs before advancing.

## Coding Guidance
- Prefer small, local, backward-compatible edits unless a broader refactor is requested.
- Follow the existing naming style and object model.
- Keep docstrings concise.
  Preferred sections are `Parameters:`, `Returns:`, `Attributes:`, and `Methods:` when useful.
- Avoid rewriting stable scientific logic just to make it look cleaner.
- If you add a new public concept, place it near the relevant package layer instead of creating parallel abstractions.

## Validation
- For touched Python files, a lightweight check is:   `python -m py_compile <files>`
- If you change workflow execution logic, inspect:
  - state creation and source linkage
  - branch/workflow/job/computation hierarchy
  - path generation and propagation
  - environment and queue interactions
- If you change chemistry classes, sanity-check that structure metadata stays aligned.

## Documentation
- The main user-facing overview lives in `README.md`.
- Keep architecture and concept docs aligned with the code when introducing new abstractions.
