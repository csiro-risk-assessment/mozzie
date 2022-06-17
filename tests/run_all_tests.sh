#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

python3 ${DIR}/TestGrid.py
python3 ${DIR}/TestCellDynamicsStatic15_9_3_2.py
python3 ${DIR}/TestCellDynamicsLogistic1_1.py
python3 ${DIR}/TestCellDynamicsBeeton2_2.py
python3 ${DIR}/TestCellDynamicsMosquito23.py
python3 ${DIR}/TestPopulationsAndParameters.py
python3 ${DIR}/TestWind.py
python3 ${DIR}/TestDiffusion_1.py
python3 ${DIR}/TestDiffusion_2.py
python3 ${DIR}/TestAdvection_1.py
python3 ${DIR}/TestAdvection_2.py
python3 ${DIR}/TestAdvection_3.py
python3 ${DIR}/TestEvolveCells.py
python3 ${DIR}/TestCellDynamicsMosquito26.py
python3 ${DIR}/TestCellDynamicsDelayBase.py
python3 ${DIR}/TestCellDynamics26DelayBase.py
python3 ${DIR}/TestCellDynamicsMosquitoLogistic26Delay.py
python3 ${DIR}/TestCellDynamicsMosquitoBH26Delay.py





