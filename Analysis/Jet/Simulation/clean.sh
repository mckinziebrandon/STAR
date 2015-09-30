#!/bin/bash

echo "clean up log and script directories"
cd log/
rm *.log *.err
cd ../
cd script/
rm *.csh
cd ../
echo "log and script directories cleaned up"
