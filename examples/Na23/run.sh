#!/bin/bash

# Set the path appropriately
bigstick="../src/bigstick.x"

# Ensure the inputs files (int, sps, opme, etc) are ready before running
$bigstick < create_wfn.in
$bigstick < create_operator.in
$bigstick < create_strength.in


rm -rf *.bigstick fort.* ham.dat

echo
echo "DONE"
echo "See nucleus_s.res for the results"
