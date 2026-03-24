# Concepts

## Source
A source is the primary chemistry object attached to a `System`.
In practice, a source is often a molecule-like object or a crystal/cell object from which states and workflows are derived.
Sources may be created from native SCOPE objects, `.xyz` files, RDKit molecules, or cell2mol objects.

## State
A `State` stores a concrete structural and data snapshot associated with a source.
States are used to represent starting points, optimized geometries, excited configurations, or other workflow results.
The `"initial"` state is commonly used as the starting reference for later computations.

## System
A `System` is the top-level project container.
It groups sources, states, workflows, environment settings, and generated data for one study.
Users commonly save and reload systems as binary files while workflows progress.

## Workflow Hierarchy
SCOPE organizes execution through four nested concepts:
- `Branch`: an alternative line of work within a system
- `Workflow`: a grouped scientific procedure inside a branch
- `Job`: a set of related execution steps inside a workflow
- `Computation`: a concrete software run

This hierarchy helps track both structure evolution and calculation history.
During execution, SCOPE may submit pending computations, wait for running ones, or register finished jobs before moving to later steps.

## Environment And Queue
An `Environment` stores filesystem and software configuration.
A `Queue` stores scheduler-related options for HPC execution.
Together they describe where jobs run and how they are submitted.

## Molecules, Cells, And Fragmentation
SCOPE handles both isolated molecules and molecule-based crystals.
Cell-related logic may fragment a crystal into molecular units, identify reference molecules, and reconstruct chemically meaningful substructures.
This logic is important for downstream analysis and should be treated carefully.

## Registration
Registration is the process of taking parsed calculation results and attaching them to SCOPE objects.
This usually means updating states, computations, or data containers so later steps can reuse the results.

## Data, Collection, And VNM
These classes organize results after or alongside registration.
They store parsed quantities, derived values, and grouped analysis data.

## Scope Input Files
Scope input files describe tasks that `scope run` should carry out on a system.
They connect user intent to the workflow hierarchy and drive automated execution step by step.

## Azo Concepts
The `azo` add-on introduces azobenzene-specific workflows and analysis.
Common ideas include:
- trans and cis forms
- transition-state or isomerization workflows
- photostationary-state style analysis
- wavelength-dependent or cross-section-driven metrics

## SCO Concepts
The `sco` add-on introduces spin-crossover-specific workflows and analysis.
Common ideas include:
- high-spin and low-spin states
- reference molecules and reference cells
- structure comparisons between spin manifolds
- energetic or structural descriptors relevant to spin switching
