#!/bin/sh

# script expects that createMeshes.sh was executed first

for setup in $(ls simpleFoam/); do
    # copy files to run folder
    cp -r simpleFoam/$setup/hydro_steady* run/$setup/
    cp simpleFoam/$setup/runHydroSteady.sh run/$setup/
    # execute run script for each case
    ./run/$setup/runHydroSteady.sh
done;
