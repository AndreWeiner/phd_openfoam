#!/bin/bash

for case in $(ls -d hydro_steady_*[0-9]); do
    ./$case/Allrun;
done
