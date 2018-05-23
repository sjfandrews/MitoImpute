#!/bin/bash
cat ~/GitCode/MitoImpute/metadata/b37-StrandFiles.txt | while read line
do
  foo="http://www.well.ox.ac.uk/~wrayner/strand/"
  bar="platforms/"
  printf -v foo "%s$line" $foo
  printf -v bar "%s${line//.zip}" $bar
  wget $foo
  mkdir $bar
  mv $line $bar
done