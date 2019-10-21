#!/bin/sh

for case in $(ls -d run/*/hydro_steady_*); do
   ./$case/Allclean
done;
