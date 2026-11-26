# Code review rules

In addition to existing rules outlined in AGENTS.md and Code-style.md, the following includes rules that should be verified when you are an AI agent performing a code review task on pull requests within this repository.

Deterministic lint rules are checked on changed lines by the `Fortran lint` CI job (Fortitude, configured in `fortitude.toml`).
Do not spend review effort re-verifying them; read the CI result, or run `python3 .github/scripts/fortitude_lint.py` locally.
Other rules which are not checked or partially checked are marked herein and still need review.

1. Use-statements should be at the lowest necessary scope

2. Subroutines and functions should be labeled `pure` when they satisfy the "pure" requirements defined in the Fortran standard.

3. Use Modern Fortran syntax when applicable (e.g. `>` instead of `.gt.`). *(lint for relational operators, `(/ /)` array literals, `double precision`, and obsolescent features)*

4. Make sure namelist XML files have the standard comment block at the beginning.

5. Make sure all `allocate` calls return `stat` and `errmsg`,
   and if there is a non-zero `stat` in the physics scheme it returns with correctly set `ccpp_error_code` and `ccpp_error_message` variables.
   *(lint only checks that a given `stat` comes with `errmsg` and is checked; an `allocate` with no `stat` at all is not flagged.)*

6. Make sure all citations have a DOI, or if not, as close to a full citation as possible.

7. Check that all full-sentence comments are grammatically correct.

8. Convert any subroutines that have only a single output variable and no side effects to functions.

9. Make sure all character-type subroutine arguments are `len=*`, unless they are declared as an allocatable or pointer in the subroutine itself.

10. Make sure all `goto` and `go to` statements have been removed. *(lint covers computed `goto` and labelled `do` loops only; plain `goto` is not flagged.)*

11. Request the removal of any lines of code that are commented out unless that code is explicitly part of a larger comment block
    describing currently working code.  If the comments aren't clear then simply flag it for human review.

12. Check all variables in modified `*.meta` files and flag any variables
    where the `standard_name` attribute is the same, but the `units` attribute is different.
    While it may not be incorrect, it could indicate the
    possibility of a science bug, and thus should be flagged so that
    a human can look it over.

13. If the review is for a module in a scheme directory, flag for manual review any `use` statements referring to modules outside the scheme directory that are not enumerated in the AGENTS.md exceptions.
