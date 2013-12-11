#!/bin/bash
FILES=INPUTMCTUPLE/*
for f in $FILES
do
  filename=$(basename "$f")
  extension=${filename##*.}
  filename=${filename%.*}
  echo "Processing $filename file..."
  mkdir -p OUTPUTPLOTS
  mkdir -p STDOUT
  sed 's:INPUTMCTUPLE:INPUTMCTUPLE/'"$filename"':' <ThreePhotonPlotsv3.cpp >ThreePhotonPlotsv3TEMP.cpp
  sed 's:OUTPUTPLOTS:OUTPUTPLOTS/'"$filename"'_plots:' <ThreePhotonPlotsv3TEMP.cpp >ThreePhotonPlotsv3TEMP.cpp
  outputfile="STDOUT/"$filename"_stdout.txt"
  root -q -l ThreePhotonPlotsv3TEMP.cpp++ > $outputfile
# testing line below
#  sed 's:INPUTMCTUPLE:OUTPUTPLOTS/'"$filename"'_plots:' <test.txt >$outputfile
done