# IMPORTANT NOTE - This repository is currently under development and items contained within it should be used with caution.

# NCAR CCPP-enabled Atmospheric Physics

Physics packages which are CCPP_enabled (along with supporting routines) are stored in this repository.

# Requirements for ccpp physics schemes

* Physics schemes contain two files, a .F90 file which contains the routines and a .meta file which contains the CCPP metadata.

* Standard names for metadata can be found [here](https://github.com/gold2718/ccpp-framework/wiki/Metadata-standard-names)

* Physics modules may contain one or more of the following routines (where XXX is the name of the module containing the physics)
  * XXX_init - contains the startup/restart initialization required for the package.  This is unthreaded and may contain I/O of undistibuted data (e.g. parameter tables)
  * XXX_timestep_init - contains the initialization required for each timestep.  This is unthreaded and may contain I/O of undistibuted data (e.g. parameter tables)
  * XXX_run - the code which is run every timestep.  Threading is controlled via the host model and the number of available threads will be passed into the routine.  No I/O is permitted in this routine
  * XXX_timestep_final - contains any post timestep handling.  This is unthreaded
  * XXX_final -  contains any code logic which needs to be performed at the end of the run prior to the program exiting.  This is unthreaded.

# Steps to adding a physics package to the repository
  * Currently under discussion

# CCPP resources

* [ccpp framework repository](https://github.com/NCAR/ccpp-framework)

# Testing

## MUSICA Chemistry

To build and test chemistry schemes, a dockerfile is provided.

From the root directory, build the docker file

```
docker build -t atmo -f test/docker/Dockerfile.musica .
```

then run it

```
docker run --rm -it atmo
```