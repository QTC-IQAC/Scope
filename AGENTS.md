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
- `Environment` stores runtime, filesystem, and scheduler information.
- Results are stored as `Data` and `Collection` classes
- Some quantum-chemistry (QC) objects exists, such as `VNM` or `ExcitedState`. They store parsed or derived results for later analysis.
- Workflow execution follows this hierarchy: `System -> Branch -> Workflow -> Job -> Computation`

## Objects are interconnected
- In principle, all objects associated with a `System` are interconnected.
- Chemistry-aware objects can be navigated both top-down and bottom-up 
- Workflow-related objects can be navigated both top-down and bottom-up 
- Sources have `States`. `State._source` links a state to its parent source.
- Sources are connected to the computational workflow through `Workflow.source`
- Computations that modified `States` are stored in `State.computations`
- State, Branch, Workflow and System have a dictionary in `self.results` where results can be stored.

## Key Invariants
- Keep chemistry containers (e.g. Molecule, Ligand, Group, Atom) internally consistent.
  `labels`, coordinates, atoms, and adjacency/connectivity data should describe the same structure.
- When attaching a source to a `System`, keep parent links and source registration in sync.
- Preserve the expected relationship between a `State` and its source.
- Path propagation across `System`, branches, workflows, jobs, and computations is part of the workflow model.
  Avoid changing it unless the task is explicitly about path handling.

## Domain Notes
- Fragmentation, reconstruction, and reference-molecule logic are scientifically sensitive.
- Transition-metals and ligands are central in species and cell operations.
- Gaussian and Quantum Espresso parsers feed results back into SCOPE objects.
  First, parsers add information to the G16_output and QE_output classes
  Second, the output classes feed the State class.
- Scope input files and `run_task` drive workflow execution step by step.
  A run may submit pending computations, skip still-running work, or register finished jobs before advancing.
- The `azo` and `sco` packages add domain-specific workflows and assumptions on top of the core model.

## Coding Guidance
- Prefer small, local, backward-compatible edits unless a broader refactor is requested.
- Follow the existing naming style and object model.
- Before creating a new helper function, search the repository for existing functions with the same or closely related behavior. Pay particular attention to shared modules such as `other.py`, `connectivity.py`, and `operations/`.
- Reuse or extend an existing function when its responsibility and semantics match the required behavior. Do not create private wrappers that duplicate existing SCOPE functions or standard-library/dependency utilities.
- When a new helper is necessary, ensure that it represents a distinct operation rather than merely renaming or forwarding another function. After implementing it, audit the new private functions against the rest of the repository.
- Prefer single-line function calls, function signatures, and assignments, even when the resulting line is long. Do not split an argument list across several lines solely to satisfy a conventional line-length limit.
- Keep logically complete operations visible on one physical line whenever practical. Use multiline formatting only for structures whose contents are genuinely easier to understand vertically, such as large dictionaries, lists, or matrices.
- When successive assignments have the same structure or purpose, align their equal signs and corresponding values into visual columns. Preserve this alignment when modifying an existing block, but do not force column alignment across unrelated statements.
- Keep docstrings concise.
  Preferred sections are `Parameters:`, `Returns:`, `Attributes:`, and `Methods:` when useful.
- Avoid rewriting stable scientific logic just to make it look cleaner.
- If you add a new public concept, place it near the relevant package layer instead of creating parallel abstractions.
- `debug` must be the last optional attribute to be passed to a function, and it must default to 0

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
