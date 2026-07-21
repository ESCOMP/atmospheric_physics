# atmospheric_physics AI agent guidance

This repository (`ESCOMP/atmospheric_physics`) is the unified science codebase of CCPP-compliant atmospheric physics parameterizations, shared between **CAM** (submoduled at `src/atmos_phys`) and **CAM-SIMA** (submoduled at `src/physics/ncar_ccpp`). The default branch is `main` (the `development` branch is deprecated).

Schemes here must remain **portable**: usable by any CCPP-enabled host model, free of host-model code. Host-side concerns (registry, build system, snapshot testing) are documented in CAM-SIMA's `AGENTS.md`; porting considerations from CAM are in CAM-SIMA's `PORTING.md`. Full conversion and usage documentation: https://escomp.github.io/CAM-SIMA-docs/.

## If you remember nothing else

1. **No host-model code in schemes.** Never `use` a CAM or CAM-SIMA module (see [Portability](#portability) for the allowed list and the `sima_diagnostics` exception).
2. **State answer impact in every PR.** The default expectation for ports, refactors, and cleanups is bit-for-bit: do not reorder floating-point arithmetic, change parenthesization, or "simplify" expressions as a side effect of other work.
3. **Never invent standard names silently.** Reuse blessed names (see [Standard names](#standard-names)); explicitly flag any new name in the PR description.
4. **Errors return, never abort.** Set `errmsg`/`errflg` and `return`; there is no `endrun` here.
5. **Keep `.meta` and Fortran in sync.** Every argument change must be mirrored in the companion `.meta` file.

## AI disclosure

atmospheric_physics **requires** disclosure of AI-assisted contributions:

- In PR descriptions: state which model(s) were used and how.
- In commits made by agents or incorporating AI-assisted code: use `Assisted-by: model:version` (e.g., `Assisted-by: claude-opus:4.8`)
- The human contributor remains responsible for all committed code at all times.

## Portability

- Schemes may `use` only `ccpp_kinds`, `ccpp_constituent_prop_mod`, `ccpp_const_utils`, intra-package dependencies, and legacy utilities under `to_be_ccppized/` (excluding `ppgrid`).
- Real kind is `kind_phys` from `ccpp_kinds`, never `r8`/`shr_kind_mod`.
- Anything host-side a scheme needs (physical constants, fields, the log unit) arrives as a subroutine argument threaded by standard name, not via a `use` of a host module.
- Dummy arguments are assumed-shape (`(:,:)`), never sized with host parameters like `pcols`.
- Module-level state: minimal; prefer arguments. Module-level vertically-dimensioned arrays cannot be statically sized (the vertical dimension is a runtime value in CAM-SIMA): make them allocatable and allocate in `_init`.
- Logging: never `write` to a hardcoded or host-provided unit directly; thread the `log_output_unit` standard name as an integer argument if a scheme must print.
- **Exception:** `schemes/sima_diagnostics/` is a non-portable directory centralizing CAM-SIMA-specific history output. Code outside `sima_diagnostics/` that calls `cam_history` routines is a bug.
- `to_be_ccppized/` holds shared legacy utilities pending full CCPPization; do not add new code there without discussion.

## Key CCPP scheme conventions

- Subroutine names are `<scheme_name>_<phase>`; generally, the module name is the scheme name, although one file may contain multiple schemes (e.g., `physics_tendency_updaters.F90`).
- Valid phases: `register`, `init`, `timestep_init`, `run`, `timestep_final`, `final`. All phases are optional.
- Every phase subroutine has `errmsg` and `errflg` (or `errcode`) as `intent(out)`: standard names `ccpp_error_message` / `ccpp_error_code`. These are required in every scheme.
- Two required Doxygen lines precede each subroutine. Current form uses the `arg_table_` prefix on the html:
  ```
  !> \section arg_table_<scheme>_<phase> Argument Table
  !! \htmlinclude arg_table_<scheme>_<phase>.html
  ```
  Older schemes omit the `arg_table_` prefix on the html file but you should use the new convention.
- Every CCPP scheme subroutine argument must have a corresponding entry in the companion `.meta` file with a valid standard name. A mismatch between the Fortran arguments and the `.meta` entries is the most common source of capgen errors.
- Constituent registration (declaring advected or non-advected constituents) happens in the `_register` phase, not `_init`. The `_register` phase runs before `_init` and before constituent indices are available. Index lookup at runtime uses `ccpp_const_get_idx` with the constituent properties object passed as an argument.

## Standard names

The standard name is the interface: the CCPP framework connects producers and consumers by exact standard-name match.

- **Reuse before coining.** Search recently converted schemes in this repo and CAM-SIMA's `src/data/registry.xml` for current, blessed usage. The official dictionary ([ESCOMP/CCPPStandardNames](https://github.com/ESCOMP/CCPPStandardNames)) is authoritative but lags current usage.
- If a new name is unavoidable, follow the naming patterns of existing names and **flag it explicitly in the PR description** for SE review. Units must be consistent with existing usage.
- The horizontal dimension standard name depends on phase: `horizontal_dimension` in non-`run` phases (`init`, `timestep_init`, ...), `horizontal_loop_extent` in the `run` phase.
- Vertical dimensions: `vertical_layer_dimension` (layers) and `vertical_interface_dimension` (interfaces).

## Code style

Read `Code-style.md` before writing or modifying Fortran.

## Code reviews

Read `Code-review.md` for specific guidelines on code reviews when performing a code review task on pull requests within this repository.

## Testing

Build and run the pFUnit unit tests (also run by CI with coverage on every PR):

```
cmake -DCMAKE_PREFIX_PATH=<pfunit>/build/installed -DATMOSPHERIC_PHYSICS_ENABLE_TESTS=ON -S ./test/unit-test -B ./build
cd build && make && ctest -V
```

- pFUnit tests live at `test/unit-test/tests/<path>/test_<module>.pf`, mirroring the scheme location.
- Test (i.e., single physics scheme) SDFs are located in `test/test_suites`; production (i.e., CAM4, CAM5, CAM7) SDFs are located in `suites`.
- Full-model validation (snapshot testing against CAM) runs from CAM-SIMA and requires NCAR cluster access; see CAM-SIMA's `AGENTS.md`. It cannot be run from this repo alone.

## PRs and commits

- PRs target `main`. Branch from a personal fork. Never push branches under `ESCOMP`.
- Feature branches are squash-merged. Do not add version tags to commit subjects; tags are applied at release time.
- Follow `.github/pull_request_template.md`, which requires listing changed files, test results, and an explicit statement of whether the PR changes answers.
- `doc/ChangeLog` and `doc/ChangeLog_template` are deprecated (history lives in git); do not add entries.

## Example code to scope for conventions

Learn the current conventions from recently converted schemes, which are preferred over older schemes that predate convention updates:

* Minimal CCPP scheme: `schemes/cloud_fraction/convective_cloud_cover.{F90,meta}`
* Scheme with namelist XML: `schemes/vertical_diffusion/vertical_diffusion_sponge_layer.{F90,meta}`, `schemes/vertical_diffusion/vertical_diffusion_sponge_layer_namelist.xml`
* Diagnostics scheme: `schemes/sima_diagnostics/scheme_diagnostics_template.F90`
