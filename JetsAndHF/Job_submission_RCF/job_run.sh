#!/bin/bash

echo ${1}
echo ${2}
echo ${3}

cd ${1}

/eic/u/rongrong/bin/eic-shell  << EOF
root -b -q 'analysis.C("'"${2}"'","'"${3}"'")'
exit
EOF
