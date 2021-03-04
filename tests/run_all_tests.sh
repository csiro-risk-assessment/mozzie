#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

python ${DIR}/TestGrid.py
python ${DIR}/TestCellDynamicsStatic15_9_3_2.py
python ${DIR}/TestCellDynamicsLogistic1_1.py
python ${DIR}/TestCellDynamicsBeeton2_2.py
python ${DIR}/TestCellDynamicsMosquito23.py
python ${DIR}/TestPopulationsAndParameters.py
python ${DIR}/TestWind.py
python ${DIR}/TestDiffusion_1.py
python ${DIR}/TestDiffusion_2.py
python ${DIR}/TestAdvection_1.py
python ${DIR}/TestAdvection_2.py
python ${DIR}/TestEvolveCells.py
python ${DIR}/TestCellDynamicsMosquito26.py


