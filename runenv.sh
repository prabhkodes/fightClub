docker build -t atmosphere-model:latest .
docker rm atmosphere-model
docker run --rm -it --name atmosphere-model -v $PWD/code/serial:/app atmosphere-model:latest