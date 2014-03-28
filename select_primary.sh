#!/bin/bash

filename=$1

awk '/\G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0/{f=1;next}/\,/{f=0}f' $filename > tmp

mv tmp filename1

sed '/*/d' filename1 > filename2
sed '/Step#/d' filename2 > filename3


mv filename3 ${filename}_primary.txt