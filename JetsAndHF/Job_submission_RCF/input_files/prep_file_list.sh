#!/bin/bash

# change the following two lines to get the sample you need
datadir=/volatile/eic/EPIC/RECO/25.10.0/epic_craterlake/DIS/NC/10x100/minQ2=1
prod=DIS_10x100_minQ2_1_25.10.0

filelistall=file.all.list
cp blank $filelistall
files=`xrdfs root://dtn-eic.jlab.org  ls $datadir`
for file in $files; do
    echo root://dtn-eic.jlab.org/$file >> $filelistall
done

split -l 20 --numeric-suffixes --suffix-length=3 $filelistall --additional-suffix=.list subList_

if [ ! -d $prod ]; then
    mkdir -pv $prod
fi
rm $prod/*
mv $filelistall $prod
mv subList* $prod
