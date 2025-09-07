#!/bin/zsh

mkdir -p Figures
# Set runtime based on USE_GPU
if [ "$1" = "GPU" ]; then
  docker compose up --build -d
else
  sed -i 's/        reservations: #/#         reservations:/g' ./docker-compose.yml
  sed -i 's/          devices: #/#           devices:/g' ./docker-compose.yml
  sed -i 's/            - driver: nvidia #/#             - driver: nvidia/g' ./docker-compose.yml
  sed -i 's/              capabilities: \[gpu\] #/#               capabilities: \[gpu\]/g' ./docker-compose.yml
  sed -i 's/              device_ids: \["0"\] #/#               device_ids: \["0"\]/g' ./docker-compose.yml
  docker compose up --build -d
fi


