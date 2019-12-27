#!/bin/bash

pst_version="0.0.1"

clear
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "... pstools package@${DIR}- start version ${pst_version} ..."

export PYTHONPATH="${DIR}/src/"
