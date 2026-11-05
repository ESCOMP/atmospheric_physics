# Notes on running the TEMPO microphysics scheme

The TEMPO microphysics scheme requires the following four lookup tables, stored as Fortran binary files:

1. CCN activation table (`ccn_activate.bin`)
2. Rain-graupel collection table (`qr_acr_qg_data_tempo_v3`)
3. Rain-snow collection table (`qr_acr_qs_data_tempo_v3`)
4. Freeze water collection table (`freeze_water_data_tempo_v3`)

You can download these tables from the [release assets](https://github.com/NCAR/TEMPO/releases/tag/tempo_v3.1.0) of the TEMPO [GitHub repository](https://github.com/NCAR/TEMPO).

Alternatively, you can build the standalone version of the TEMPO microphysics scheme from source, then run the resulting `run_build_tables` executable to generate the tables yourself.

Once you have the tables, either set their paths with the corresponding namelist options in `user_nl_cam`, or place them in the CAM-SIMA `run` directory. See [`mp_tempo_namelist.xml`](mp_tempo_namelist.xml) for details.

Fortran binary files are sensitive to byte order, and the release assets are in big-endian format. On supported machines defined in `ccs_config`, CAM-SIMA is built with compiler flags that enforce big-endian I/O for Fortran binary files. However, if you use custom `ccs_config` for your own machine, make sure to enable the same behavior. Otherwise, the tables may fail to load, or worse, load silently with garbage values.
