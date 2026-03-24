# Architecture

## Overview
SCOPE is organized around a chemistry model, a workflow model, an analysis model, and a runtime model.
The core package handles generic molecule and crystal workflows, while `azo` and `sco` add domain-specific behavior.
Users commonly build and update serialized `System` objects through the CLI and inspect them later in notebooks.

## Package Layout
- `core/src/scope/`: main package, CLI, chemistry classes, workflows, parsers, and utilities
- `azo/src/scope_azo/`: azobenzene-oriented extensions and CLI helpers
- `sco/src/scope_sco/`: spin-crossover-oriented extensions and CLI helpers

## Chemistry Layer
The chemistry layer represents structures and their relationships.
Common classes include:
- `Atom`
- `Specie` and derived structure classes
- `Cell`
- `State`

These classes hold structural data, connectivity, and chemical metadata used by later workflow and analysis stages.
Sources can be created from native SCOPE objects, `.xyz` data, RDKit molecules, or cell2mol structures.

## Workflow Layer
The workflow layer organizes computations from project level down to individual runs.
Its main hierarchy is:

`System -> Branch -> Workflow -> Job -> Computation`

- `System` is the top-level container for sources, branches, and environment configuration.
- Branches separate alternative workflow paths for a system.
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

They manage paths, software configuration, and scheduler-specific settings for HPC usage.

## Parsing And Registration
Software-specific modules under `core/src/scope/software/` generate inputs and parse outputs for Gaussian and Quantum Espresso.
Parsed results are then registered back into SCOPE objects, typically through states, workflow objects, and data containers.
The analysis layer also includes `Data`, `Collection`, and `VNM` objects for storing parsed or derived results.

## Data Flow
The usual execution path is:
1. create or load a source structure
2. attach it to a `System`
3. generate an initial `State`
4. create branches, workflows, jobs, and computations
5. write inputs and run software
6. parse outputs
7. store results back in SCOPE objects for inspection and further steps

Persistent outputs commonly include binary `System` files plus the corresponding computation directory structure.

## Add-Ons
- `azo` extends SCOPE with azobenzene-specific structure handling and workflow utilities
- `sco` extends SCOPE with spin-crossover-specific structures, analysis, and workflow helpers

These add-ons build on the same core system, state, and workflow abstractions instead of replacing them.
