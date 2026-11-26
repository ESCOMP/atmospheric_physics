# Fortran code style

Distilled, agent-relevant rules from the [CAM coding standards](https://escomp.github.io/CAM-SIMA-docs/development/cam-coding-standards/), which apply to Fortran in both atmospheric_physics and CAM-SIMA. Legacy code may not comply; follow these standards for new and modified lines, but do not reformat code you are not otherwise changing.

Deterministic rules are checked on changed lines by Fortitude (`fortitude.toml`) in CI; run `python3 .github/scripts/fortitude_lint.py` to check them locally.

The below rules still require human (agent) attention.

## MUST

- Functions must not have side effects and should carry the `pure` keyword; if a function cannot be `pure`, say why in its preamble.
- Initialize local pointers (by default with `nullify`) before any non-initialization statement.
- Namelist variables (except logicals) are initialized to invalid sentinels: integer `-HUGE(1)`, real `-HUGE(1._kind_phys)`, character `'UNSET'`.
- Optional arguments are passed by keyword: `call sub(x, opt_y=y)`, never positionally.
- Bring `use` statements in at the smallest scope: inside subroutines, unless the symbol is needed by module-level declarations.

## SHOULD

- Avoid preprocessor directives (`#if`, `#ifdef`); prefer runtime logic.
- Break long formulas into readable pieces with temporary variables.
- Avoid pointers as dummy arguments.
- Module name matches the filename (minus `.F90`).
- One dummy argument declaration per line (related items may be grouped); declaration order matches the argument list.
- Avoid unnecessary statements (e.g., a bare `return` at the end of a subroutine).

## Layout

- Indentation follows scope (3 spaces recommended; each module at least self-consistent).
- Continuation lines indented 5 spaces or aligned with similar lines in the statement.
- Spaces around operators and `::`, after commas and after `if`/`else`/`end`/`do`/`while`; `only:` with no space before the colon.

## Comments

- Explain the purpose or non-obvious rationale of the following code; do not restate the code logic.
- When modifying code, check that adjacent comments remain correct.
- Do not keep commented-out code "for later".
