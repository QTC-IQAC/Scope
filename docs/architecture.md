# Architecture

## Overview
SCOPE is organized around a chemistry model, a workflow model, an analysis model, and a runtime model.
The core package handles workflows for generic molecules and molecule-based crystals, while `azo` and `sco` add domain-specific behavior.
Users commonly build and update `System` objects through the CLI and inspect them later in notebooks.

## Package Layout
- `core/src/scope/`:    main package, CLI, chemistry classes, workflows, output parsers, and utilities
- `azo/src/scope_azo/`: azobenzene-oriented extensions and CLI helpers
- `sco/src/scope_sco/`: spin-crossover-oriented extensions and CLI helpers

## Chemistry Layer
The chemistry layer represents structures and their relationships.
The two main classes are:
- `Molecules` (subclass of `Specie`) 
- `Cell`

Substructures are automatically recognized, and include:
- `Ligand`
- `Group`
- `Atom`
- `Metal`

Normally, the user works with the main classes (Molecule and Cell), while substructures are used internally by SCOPE to execute tasks. 
Those classes hold structural data, connectivity, and chemical metadata used by later workflow and analysis stages.
They must be sourced to a System in order to execute a computational workflow. There are System-class functions to facilitate this task.
Species can be created from `.xyz` data, RDKit molecules, or cell2mol molecules.
Cells can be created from `.xyz` data and cell parameters, or cell2mol Cell objects.

## Workflow Layer
The workflow layer organizes computations from project level down to individual runs.
Its main hierarchy is:

`System -> Branch -> Workflow -> Job -> Computation`

- `System` is the top-level container for sources, branches, and file management.
- Branches connect different workflows of a system. Branches are associated with a folder were all computations will be located.
- Workflows group related jobs and are associated with a specific source of the system.
- Jobs group one or more computations.
- Computations represent concrete software runs and their inputs, outputs, and status.

Scope input files describe which actions should be taken on a system.
The CLI command `scope run` and the underlying `run_task` logic traverse this hierarchy, decide what can be submitted or registered, and advance the workflow incrementally.

## Runtime Layer
The runtime layer describes where and how computations are executed.
Key classes include:
- `Environment`
- `Queue`
- `Node`

The `Environment` class manages paths, software configuration, and scheduler-specific settings for HPC usage.
`Queue` and `Node` are mainly used for queue selection (if more than one is available) when submitting jobs

## Parsing And Registration
Software-specific modules under `core/src/scope/software/` generate inputs and parse outputs for Gaussian and Quantum Espresso.
Parsed results are then registered back into SCOPE objects, typically through states and workflow objects.
Data is typically organized in two classes: `Data` and `Collection` with dedicated attributed and methods.
Some QC objects have dedicated classes, such as `VNM` (Vibrational Normal Modes) and `ExcitedState`. 

## Data Flow
The usual execution path is:
1. configure the `Environment` associated with a project (`scope config`).
2. create a `System`.
3. add sources (`Molecule` or `Cell`).
4. define the desired workflow through input files.
5. execute workflow (`scope run`). Run as many times as necessary
6. results are stored in binary `System` files. Download if necessary and analyse

Persistent outputs include those binary `System` files, and the computations. Each are located under the main project path, in dedicated folders.

## Add-Ons
- `azo` extends SCOPE with azobenzene-specific structure handling and workflow utilities
- `sco` extends SCOPE with spin-crossover-specific structures, analysis, and workflow helpers

These add-ons build on the same core system, state, and workflow abstractions instead of replacing them.