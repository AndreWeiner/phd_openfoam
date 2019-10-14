#!/bin/sh

cp -r snappyHexMesh/* run/

for case in $(ls -d run/*/refinement_[0-9]); do
  ./$case/Allrun
done;
