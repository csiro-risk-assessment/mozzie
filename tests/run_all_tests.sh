#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

python3 ${DIR}/testGrid.py
python3 ${DIR}/testSpatialDynamics.py
python3 ${DIR}/testCellDynamicsStatic15_9_3_2.py
python3 ${DIR}/testCellDynamicsLogistic1_1.py
python3 ${DIR}/testCellDynamicsBeeton2_2.py
python3 ${DIR}/testCellDynamicsMosquito23.py
python3 ${DIR}/testPopulationsAndParameters.py
python3 ${DIR}/testWind.py
python3 ${DIR}/testDiffusion_1.py
python3 ${DIR}/testDiffusion_2.py
python3 ${DIR}/testAdvection_1.py
python3 ${DIR}/testAdvection_2.py
python3 ${DIR}/testAdvection_3.py
python3 ${DIR}/testEvolveCells.py
python3 ${DIR}/testCellDynamicsMosquito26.py
python3 ${DIR}/testCellDynamicsDelayBase.py
python3 ${DIR}/testCellDynamics26DelayBase.py
python3 ${DIR}/testCellDynamicsMosquitoLogistic26Delay.py
python3 ${DIR}/testCellDynamicsMosquitoBH26Delay.py
python3 ${DIR}/testCalcQm.py






