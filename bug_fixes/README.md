# Bug-fix notebooks

This directory contains reproducible debugging cases and notebooks used during SCOPE development.

Cases are grouped by the SCOPE version in which the bug was identified. For example, `0.9.5/` contains bugs discovered while reviewing SCOPE 0.9.5. 

It is intentionally located outside the installable packages. The core package is discovered only under `core/src`, while the add-ons are discovered under their respective `src` directories. Files in `bug_fixes/` are therefore not included when SCOPE is installed through pip.