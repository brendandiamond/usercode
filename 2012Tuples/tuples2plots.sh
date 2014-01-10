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
  mkdir -p OUTPUTHIST
  # note that "" after -i is only for MAC OS!
  sed 's:INPUTMCTUPLE:INPUTMCTUPLE/'"$filename"':' <ThreePhotonPlotsv3.cpp >ThreePhotonPlotsv3TEMP.cpp
  sed -i "" 's:OUTPUTPLOTS:OUTPUTPLOTS/'"$filename"'_plots:' ThreePhotonPlotsv3TEMP.cpp
  sed -i "" 's:OUTPUTHIST:OUTPUTHIST/'"$filename"':' ThreePhotonPlotsv3TEMP.cpp
  sed -i "" 's:ThreePhotonPlotsv3:ThreePhotonPlotsv3TEMP:' ThreePhotonPlotsv3TEMP.cpp 
  outputfile="STDOUT/"$filename"_stdout.txt"
  root -q -l ThreePhotonPlotsv3TEMP.cpp++ > $outputfile
# testing line below
#  sed 's:INPUTMCTUPLE:OUTPUTPLOTS/'"$filename"'_plots:' <test.txt >$outputfile
done