#!/bin/bash

# rh
grep '# subjectname' rh.aparc.stats/*.stats |
sed -r 's/.*([0-9]{7}).*/\1/' > rh.aparc.SUB_ID.txt

grep '# Measure Cortex, NumVert' rh.aparc.stats/*.stats |
sed -r 's/.*Vertices, ([0-9]+), .*/\1/' > rh.aparc.NumVert.txt

grep '# Measure Cortex, WhiteSurfArea' rh.aparc.stats/*.stats |
sed -r 's/.*Area, ([0-9]*\.?[0-9]*), .*/\1/' > rh.aparc.WhiteSurfArea.txt

grep '# Measure Cortex, MeanThickness' rh.aparc.stats/*.stats |
sed -r 's/.*Thickness, ([0-9]*\.?[0-9]*), .*/\1/' > rh.aparc.MeanThickness.txt

# lh
grep '# subjectname' lh.aparc.stats/*.stats |
sed -r 's/.*([0-9]{7}).*/\1/' > lh.aparc.SUB_ID.txt

grep '# Measure Cortex, NumVert' lh.aparc.stats/*.stats |
sed -r 's/.*Vertices, ([0-9]+), .*/\1/' > lh.aparc.NumVert.txt

grep '# Measure Cortex, WhiteSurfArea' lh.aparc.stats/*.stats |
sed -r 's/.*Area, ([0-9]*\.?[0-9]*), .*/\1/' > lh.aparc.WhiteSurfArea.txt

grep '# Measure Cortex, MeanThickness' lh.aparc.stats/*.stats |
sed -r 's/.*Thickness, ([0-9]*\.?[0-9]*), .*/\1/' > lh.aparc.MeanThickness.txt
