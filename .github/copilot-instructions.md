# Copilot Instructions

- Follow the existing SCOPE object model and naming conventions.
- Prefer minimal, targeted edits over broad refactors.
- Keep docstrings concise.
  Use a short summary plus `Parameters:` and `Returns:` only when they add value.
- Preserve backward compatibility unless a breaking change is explicitly requested.
- Keep the workflow hierarchy intact:
  `System -> Branch -> Workflow -> Job -> Computation`
- Remember that serialized `System` objects are a common user-facing artifact.
- Preserve the relationship between systems, sources, states and workflows. 
- Similarly, preserve the relationship between environments, queues and nodes.
- Keep support for common source imports such as `.xyz`, RDKit, and cell2mol.
- Avoid changing lazy attribute patterns or parser result structures unless necessary.
- In chemistry objects, keep labels, coordinates, atom containers, and connectivity data consistent.
- Keep `Data` and `Collection` compatible with the parser and workflow layers.
- In `azo` and `sco`, preserve domain-specific assumptions unless the task is about updating them.
- For quick validation of touched Python files, use:
  `python -m py_compile <files>`
- Take advantage of physical constants and unit conversions implemented in constants.py. Add new ones if necessary.
