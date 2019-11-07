#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

python ${DIR}/TestGrid.py -v
python ${DIR}/TestWind.py -v
python ${DIR}/TestDiffusion_1.py -v
python ${DIR}/TestDiffusion_2.py -v

