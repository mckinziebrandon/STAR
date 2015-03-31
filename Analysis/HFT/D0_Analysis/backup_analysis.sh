#!/bin/bash

cd ./backup

name_sep="_"
name_prefix="backup_"
name_hour=$(date +%H)
name_days=$(date +%j)
name_year=$(date +%Y)
name_full=$name_prefix$name_hour$name_sep$name_days$name_sep$name_year
echo Backup files to directory $name_full

mkdir $name_full
cd $name_full
cp ../../StStrangenessAna_macro.cc .
cp ../../StRoot/StStrangenessAna/*.* .
cd ../../

