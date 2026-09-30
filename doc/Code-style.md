# Fortran code style

Distilled, agent-relevant rules from the [CAM coding standards](https://escomp.github.io/CAM-SIMA-docs/development/cam-coding-standards/), which apply to Fortran in both atmospheric_physics and CAM-SIMA. Legacy code may not comply; follow these standards for new and modified lines, but do not reformat code you are not otherwise changing.

## If you remember nothing else

1. Every `real` variable and literal has an explicit kind: `1.5_kind_phys`, never `1.5` or `1.5d0`. Literals include the decimal point.
2. Never initialize a local variable on its declaration line — it acquires implicit `SAVE` and is not thread-safe. Initialize at the top of the executable section. (Module-level variables already have implicit `SAVE` per the Fortran standard; an explicit `save` there is redundant and being removed in refactors.)
3. No naked `use`: always `use mod, only: symbol`.
4. Modules have `implicit none` in the preamble and are default `private`; public interfaces are declared explicitly.
5. Optional arguments are passed by keyword: `call sub(x, opt_y=y)`, never positionally.
6. Bring `use` statements in at the smallest scope: inside subroutines, unless the symbol is needed by module-level declarations.

## MUST

- Spaces, not tabs; no trailing whitespace.
- `intent` on all dummy arguments except pointers.
- Fortran 90+ declaration syntax: `::` in all variable, type, and procedure declarations; `character(len=*)`-style character declarations.
- `if` statements spanning more than one line use `if ... then`; no continued single-line `if`.
- No semicolons combining statements on one line.
- Functions must not have side effects and should carry the `pure` keyword; if a function cannot be `pure`, say why in its preamble.
- Initialize local pointers (by default with `nullify`) before any non-initialization statement.
- Namelist variables (except logicals) are initialized to invalid sentinels: integer `-HUGE(1)`, real `NaN`, character `'UNSET'`.

## SHOULD

- Avoid preprocessor directives (`#if`, `#ifdef`); prefer runtime logic.
- Break long formulas into readable pieces with temporary variables.
- Use symbolic comparison operators (`==`, `/=`, `<`) rather than `.eq.`-style.
- Avoid pointers as dummy arguments.
- Module name matches the filename (minus `.F90`).
- One dummy argument declaration per line (related items may be grouped); declaration order matches the argument list.
- Avoid unnecessary statements (e.g., a bare `return` at the end of a subroutine).

## Layout

- Indentation follows scope (3 spaces recommended; each module at least self-consistent).
- Lines under 133 characters.
- Continuation lines indented 5 spaces or aligned with similar lines in the statement.
- Spaces around operators and `::`, after commas and after `if`/`else`/`end`/`do`/`while`; `only:` with no space before the colon.

## Comments

- Explain the purpose or non-obvious rationale of the following code; do not restate the code logic.
- When modifying code, check that adjacent comments remain correct.
- Do not keep commented-out code "for later".
