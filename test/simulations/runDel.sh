#!/bin/bash

for i in `seq 1 10`;
do
    make -f MakefileDEL all
    make -f MakefileDEL rest
    find DEL*out |  xargs -I {} perl processWHAMout.pl {} > DEL-res.txt
done