#!/bin/bash

for i in *.pbs
do 
    qsub $i
done