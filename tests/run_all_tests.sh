#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

python ${DIR}/TestGrid.py
python ${DIR}/TestWind.py
python ${DIR}/TestDiffusion_1.py
python ${DIR}/TestDiffusion_2.py

