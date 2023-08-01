#! /bin/bash

mkdir -p build
cd build
cmake ~/canoe -DNETCDF=ON -DRFM=ON -DDISORT=ON -DHYDROSTATIC=ON -DNVAPOR=1 -DPLANET=Earth
