#!/bin/sh

mkdir -p run
cp -r snappyHexMesh/* run/

for case in $(ls -d run/*/refinement_[0-9]); do
  ./$case/Allrun
done;
