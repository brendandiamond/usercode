#!/bin/bash
FILES=INPUTMCTUPLE/*
for f in $FILES
do
  echo "Processing $f file..."
  # take action on each file. $f store current file name
  # sed command
  cat $f
done