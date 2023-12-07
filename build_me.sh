#! /bin/bash

mkdir -p build
cd build
cmake ~/canoe -DNETCDF=ON -DRFM=ON -DDISORT=ON -DHYDROSTATIC=ON -DNVAPOR=1 -DPLANET=Earth
#cmake ~/canoe_1 -DNETCDF=ON -DRFM=ON -DDISORT=ON -DHYDROSTATIC=ON -DNVAPOR=1 -DPYTHON_BINDINGS=ON
