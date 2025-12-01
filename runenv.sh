# docker build -t fortran-mpi:latest .
docker container rm -f fortran-mpi
docker run -it --name fortran-mpi -v $PWD/code:/app fortran-mpi:latest