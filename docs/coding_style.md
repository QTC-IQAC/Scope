# SCOPE Coding Style

This document records the coding conventions used in SCOPE tutorials and development notebooks. It is descriptive rather than a strict formatter configuration. Scientific clarity and consistency with the object model take priority over abstraction.

## Notebook Structure

- Introduce the scientific or technical purpose in Markdown before executing code.
- Organize notebooks into clearly named parts and discuss one concept at a time.
- Keep code cells short and sequential so the state of the notebook is easy to follow.
- Explain the expected result near the code that produces it.
- Prefer visible intermediate values and printed SCOPE objects over hidden processing.

## Imports And Paths

- Use direct imports near the beginning of the notebook, such as `import os`, `import numpy as np`, and `import scope`.
- Import a specific SCOPE class or function only when direct access improves readability.
- Define data paths explicitly using a short descriptive variable such as `data_folder`, `system_path`, or `ajaxep_path`.
- Assume SCOPE has been installed in the active environment. Avoid repository-discovery variables and `sys.path` manipulation in user-facing notebooks.

## Working With SCOPE Objects

- Use the methods already provided by SCOPE instead of wrapping them in notebook-specific helper functions.
- Navigate objects explicitly with methods such as `find_source()`, `find_state()`, `find_branch()`, `find_workflow()`, and `get_parent()`.
- Preserve the usual `(found, object)` return pattern:

  ```python
  found, source = sys.find_source('ref_hs_mol')
  found, state = source.find_state('b3lyp_opt')
  ```

- Use concise domain names where the context is clear: `sys`, `cell`, `mol`, `lig`, `met`, `state`, `job`, and `comp`.
- Print objects directly to use their `__repr__` summaries and make hierarchy inspection visible.

## Results And Debugging

- Prefer direct comparisons, printed values, and small loops over generic result dictionaries or reporting frameworks.
- Keep a debugging calculation close to the SCOPE objects it examines.
- Use assertions for regression checks, but keep them in a dedicated cell when the notebook must also illustrate the unfixed behavior.
- When a bug depends on private data, include a small portable example whenever possible and keep the real-data reproduction as a separate section.
- Do not save or mutate persistent `System` files unless that is the explicit purpose of the notebook.

## Functions And Abstraction

- Do not create a helper function merely to rename or wrap one existing SCOPE method.
- Define functions when they represent a genuine reusable algorithm, a substantial repeated operation, or an independent scientific calculation.
- Prefer explicit code when it makes object relationships or scientific assumptions easier to inspect.

## Core Code

- Keep changes small, local, and backward-compatible when possible.
- Maintain consistency between labels, coordinates, atoms, connectivity, and parent indices.
- Preserve bidirectional navigation through the chemistry and workflow hierarchies.
- Use concise docstrings with `Parameters:`, `Returns:`, `Attributes:`, and `Methods:` sections when useful.
- Avoid broad refactors of stable scientific logic unless the task explicitly requires them.
