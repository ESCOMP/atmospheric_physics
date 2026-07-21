# Code review rules

In addition to existing rules outlined in AGENTS.md and Code-style.md, the following includes rules that should be verified when you are an AI agent performing a code review task on pull requests within this repository:

1. No `save` statement at the module level, since module variables are implicitly `save`d.
   Subroutine level implicit `save` is usually unintentional and should be flagged.

2. Use-statements should be at the lowest necessary scope

3. Subroutines and functions should be labeled `pure` when they satisfy the “pure” requirements defined in the Fortran standard.

4. Use `end do` over `enddo`, and `end if` over `endif`.

5. Use Modern Fortran syntax when applicable (e.g. `>` instead of `.gt.`)

6. Make sure namelist XML files have the standard comment block at the beginning.

7. Make sure all `allocate` calls return `stat` and `errmsg`, and if there is a non-zero `stat` in the physics scheme it returns with correctly set `ccpp_error_code` and `ccpp_error_message` variables.

8. Make sure all citations have a DOI, or if not, as close to a full citation as possible.

9. Check that all full-sentence comments are grammatically correct.

10. Convert any subroutines that have only a single output variable to functions.

11. Make sure all character-type subroutine arguments are `len=*`, unless they are declared as an allocatable or pointer in the subroutine itself.

12. Make sure all `goto` and `go to` statements have been removed.
