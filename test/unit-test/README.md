# Unit Tests

To add or update unit tests, follow the instructions from the [development guide](https://escomp.github.io/CAM-SIMA-docs/atmospheric_physics/development_workflow/#5-unit-testing).  Also, make sure [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit) is built and installed following the [build directions](https://github.com/Goddard-Fortran-Ecosystem/pFUnit?tab=readme-ov-file#building-and-installing-pfunit) (see the associated workflow file for a detailed example).

To run the tests, from the root directory of your clone, run:

```bash
$ cmake \
    -DCMAKE_PREFIX_PATH=<path_to_pfunit>/build/installed \
    -DATMOSPHERIC_PHYSICS_ENABLE_CODE_COVERAGE=ON \
    -B./build \
    -S./test/unit-test
$ cd build
$ make
$ ctest -V --output-on-failure
```

Where `<path_to_pfunit>` is the path to your pfunit repository.  The install path of pFUnit may be different depending on how you've built your local verison.
