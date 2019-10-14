#!/bin/sh

for case in $(ls -d run/*/refinement_[0-9]); do
   ./$case/Allclean
done;
