#!/bin/bash
# License           : GPL3
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 10/02/2020
# Last Modified Date: 10/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

echo "Creating a Python $PYTHON_VERSION environment"
conda create -n hosts python=$PYTHON_VERSION || exit 1
conda activate hosts
