# NCAR CCPP-enabled Atmospheric Physics

Physics packages which are CCPP_enabled (along with supporting routines) are stored in this repository.

# Requirements for ccpp physics packages

* physics packages contain two files, a .F90 file which contains the routines and a .meta file which contains the CCPP metadata.

* Physics modules may contain one or more of the following routines (where XXX is the name of the module containing the physics)
** XXX_initl - contains the startup/restart initialization required for the package.  This is unthreaded and may contain I/O
** XXX_timestep_init - contains the initialization required for each timestep.  This is unthreaded and may contain I/O
** XXX_run - the code which is run every timestep.  Threading is controlled via the host model and the number of available threads will be passed into the routine.  No I/O is permitted in this routine
** XXX_timestep_final - contains any post timestep handling.  This is unthreaded
** XXX_final -  contains any code logic which needs to be performed at the end of the run prior to the program exiting.

# Steps to adding a physics package to the repository

# CCPP resources

* [ccpp framework repository] (https://github.com/NCAR/ccpp-framework)

