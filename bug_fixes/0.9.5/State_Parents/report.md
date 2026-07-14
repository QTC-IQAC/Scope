# State-derived molecule parents

## Summary

In SCOPE 0.9.5, molecules created from a `State` were incorrectly attached directly to the state's source. As a result, their atoms could inherit source-level parents and coordinates instead of remaining linked to the geometry of the state from which they were generated.

This was visible for the AJAXEP optimized structures: a metal atom had the correct optimized-state coordinates, while atoms returned by `get_coord_sphere()` could correspond to the source molecule rather than the optimized molecule.

## Ownership model

The applied correction uses the following relationship:

```text
Source -> State -> Molecule -> Atom
```

`State._source` remains the route from a state-derived molecule back to the original source. State-derived molecules and their atoms therefore no longer need a direct source parent.

## Applied changes

### `classes_state.py`

- `State.get_molecules()` now assigns the `State` as parent of every generated molecule.
- Molecular fractional coordinates are derived from the state's geometry and cell data, rather than from the source cell.
- `State.set_geometry()` and `State.set_cell()` rebuild state-derived molecules so their geometry and parent relationships remain synchronized.
- `State.set_geometry_from_molecules()` uses indices stored in the State parent relationship when reconstructing atom order.
- Molecules produced during `State.reconstruct()` are reattached to the State. Their atoms and complex subdivisions are then refreshed so the corrected parent relationship propagates through the molecular hierarchy.

### `classes_specie.py`

- Species can obtain fractional coordinates from a State parent that contains a cell vector.
- Atom creation recognizes the State as a valid provider of periodic-cell information.

### `reconstruct.py`

- Reconstruction obtains original cell indices from the State parent first.
- Direct Cell parents remain supported as a fallback for objects created outside a State.

## Compatibility considerations

- Access to the original source remains available through the State parent and `State._source`.
- Direct Cell parenting is retained as a fallback for independently created species.
- State geometry, fractional coordinates, molecule membership, and atom ordering are refreshed together after geometry or cell changes.

## Validation

The change was first applied to an isolated copy of SCOPE and then tested with:

- `state_parents.ipynb`, including the AJAXEP high-spin and low-spin optimized structures
- Tutorial 1 from the SCOPE tutorials repository
- a synthetic periodic-cell geometry round trip
- Python compilation checks for all modified source files

The regression notebook confirmed that coordination-sphere atoms belong to the optimized State molecule and that their coordinates match that molecule. Periodic State parents, fractional coordinates, atom ordering, and geometry round trips also passed.

Tutorial 1 completed all code cells. Two RDKit rich-display operations reported that `MolDraw2DCairo` was unavailable in the test environment; these display-only messages are unrelated to the State-parent changes.
