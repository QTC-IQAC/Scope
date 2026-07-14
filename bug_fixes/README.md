# Bug-fix notebooks

This directory contains reproducible debugging cases and regression notebooks used during SCOPE development.

Cases are grouped by the SCOPE version in which the bug was identified. For example, `0.9.5/` contains regressions discovered while reviewing SCOPE 0.9.5. A regression may remain useful after that version as a test against reintroduction of the bug.

It is intentionally located outside the installable packages. The core package is discovered only under `core/src`, while the add-ons are discovered under their respective `src` directories. Files in `bug_fixes/` are therefore not included when SCOPE is installed through pip.

Notebooks in this directory may use local or optional datasets. Portable synthetic fixtures should be included whenever possible so that the underlying regression can still be tested without those datasets.
