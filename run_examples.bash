#!/bin/bash -e

source /Applications/NGSolve/ngs-venv/bin/activate

python3 condition.py --steps 100
python3 convergence.py --example 1 --maxh 0.3 --nref 1 --vtk 1
python3 convergence.py --example 1 --maxh 0.3 --nref 6
python3 convergence.py --example 1 --maxh 0.3 --nref 6 --hogeo 0
python3 convergence.py --example 1 --maxh 0.3 --nref 6 --lm 2

python3 convergence.py --example 2 --maxh 0.3 --nref 6 --threshold 0.1
python3 convergence.py --example 2 --maxh 0.3 --nref 6 --threshold 0.1 --lm 2
python3 convergence.py --example 2 --maxh 0.3 --nref 6 --hogeo 0

deactivate