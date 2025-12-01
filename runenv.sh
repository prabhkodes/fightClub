docker build -t atmosphere-model:latest .
docker run -it --name atmosphere-model -v $PWD/code/serial:/app atmosphere-model:latest