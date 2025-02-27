MUSICA tests for CAM-SIMA Physics
=================================

To build and run the MUSICA tests for CAM-SIMA in a Docker container, from the
top-level folder run:

```
docker build -t atmos-phys . -f test/docker/Docker.musica
docker run -it atmos-phys bash
make test
```