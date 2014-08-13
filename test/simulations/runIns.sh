#!/bin/bash

for i in `seq 1 10`;
do
    make -f MakefileINS all
    make -f MakefileINS rest
    find INS*out |  xargs -I {} perl processWHAMout.pl {} > INS-res.txt
done