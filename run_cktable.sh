#!/bin/bash
/home/linfel/LavaPlanet/build/bin/cktable.py --kcoeff kcoeff.lava_planet-B7.nc --atm /home/linfel/canoe/data/lava_SiO_atm_isothermal.txt --cia "" --input kcoeff.inp-B7 --output lava_planet-B7  > log.cktable-lava_planet-B7 &

/home/linfel/LavaPlanet/build/bin/cktable.py --kcoeff kcoeff.lava_planet-B8.nc --atm /home/linfel/canoe/data/lava_SiO_atm_isothermal.txt --cia "" --input kcoeff.inp-B8 --output lava_planet-B8  > log.cktable-lava_planet-B8 &

