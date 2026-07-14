# Bug-fixing procedure

This document describes how reproducible bug investigations and their fixes are recorded in SCOPE. The aim is to preserve the original problem, explain the reasoning behind the correction, and provide a regression case that can be rerun in later versions.

## Directory structure

Bug cases are stored outside the installable packages under `bug_fixes/`. They are first grouped by the SCOPE version in which the bug was identified and then by a descriptive case name:

```text
bug_fixes/
└── <version>/
    └── <Bug_Name>/
        ├── bug_details.ipynb
        └── report.md
```

For example:

```text
bug_fixes/
└── 0.9.5/
    └── State_Parents/
        ├── bug_details.ipynb
        └── report.md
```

Use the version under review when the problem is discovered, even if the fix is applied before that version is released. Use a short, descriptive folder name with words separated by underscores.

Because `bug_fixes/` is outside `core/src` and the add-on `src` directories, these development files are not installed with SCOPE through pip.

## The bug-details notebook

`bug_details.ipynb` is the executable record of the problem. It should:

- describe the expected and observed behavior before discussing the implementation;
- reproduce the bug with the smallest useful example;
- include a portable synthetic example whenever private or external data is also required;
- inspect SCOPE objects directly using the existing API and normal object hierarchy;
- show the scientifically relevant values, relationships, or coordinates involved;
- contain assertions that verify the corrected behavior and guard against regression;
- state clearly when an optional dataset is unavailable and skip only the checks that require it;
- avoid modifying persistent user data.

The notebook should follow the conventions in `docs/coding_style.md`. In particular, use explicit scientific variables and existing SCOPE methods instead of repository-discovery machinery, generic result dictionaries, or small wrapper functions.

During diagnosis, assertions that describe the intended behavior may initially fail or be temporarily disabled. Once the correction is applied, they should be enabled by default so that executing the notebook performs the regression checks.

## The report

`report.md` documents the completed investigation. It should be concise but sufficient for a future developer to understand why the code changed. Include:

1. a summary of the bug and its visible effect;
2. the relevant object model, data flow, or scientific invariant;
3. the diagnosed cause;
4. the changes applied, grouped by source file or component;
5. compatibility considerations and retained fallback behavior;
6. the validation performed and its outcome;
7. any unrelated warnings, unavailable optional data, or environmental limitations observed during testing.

The report should describe behavior rather than merely list changed lines. If the fix changes ownership or navigation between SCOPE objects, show the intended relationship explicitly.

## Fixing workflow

Use the following sequence for a bug investigation:

1. **Reproduce the problem.** Create the versioned case folder and establish the failure in `bug_details.ipynb`.
2. **Identify the invariant.** Determine which object relationship, structural metadata, workflow behavior, or scientific result is incorrect.
3. **Trace the cause.** Inspect the smallest relevant part of the code and record how the observed result is produced.
4. **Design a local correction.** Prefer a small change consistent with SCOPE's existing object model. Consider parent navigation, atom ordering, coordinates, connectivity, and periodic information when chemistry classes are involved.
5. **Test in isolation when appropriate.** For changes with broad consequences, first apply the proposed patch to a temporary copy of SCOPE and run the regression notebook and relevant tutorials against it.
6. **Apply the correction.** Transfer the validated change to the main repository without altering unrelated work.
7. **Enable regression assertions.** Ensure `bug_details.ipynb` fails if the corrected behavior is reintroduced.
8. **Write `report.md`.** Record the cause, implementation, compatibility decisions, and test results.
9. **Run final validation.** Repeat the checks against the repository version, not only the temporary copy.

## Validation expectations

Validation should be proportional to the affected code. At minimum:

- execute `bug_details.ipynb` from beginning to end;
- compile each modified Python file with `python -m py_compile`;
- run the tutorial or workflow most closely related to the changed behavior;
- inspect relevant parent links and indices when object ownership changes;
- confirm that labels, coordinates, atoms, adjacency data, and fractional coordinates remain aligned when chemistry containers change;
- distinguish failures caused by the fix from optional dependency or display-only warnings.

When a bug uses real project data, retain both levels of validation when possible: a portable synthetic regression for future use and a real-data check demonstrating that the original problem is resolved.

## Scope of bug-fix cases

The contents of `bug_fixes/` are development records and regression material, not part of the public SCOPE API. Reusable automated checks may later be promoted into the formal test suite, but the original notebook and report should remain as the historical explanation of the bug.
