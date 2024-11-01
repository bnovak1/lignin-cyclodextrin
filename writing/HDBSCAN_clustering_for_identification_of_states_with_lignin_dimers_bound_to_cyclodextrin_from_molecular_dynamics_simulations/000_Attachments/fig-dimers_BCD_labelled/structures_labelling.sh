#!/bin/bash

# This script is used to label the chemical structures in the image

# Get the directory name from the command line
DIR=$1

# Dimers and BCD
convert -pointsize 22 \
        -fill blue -draw 'text 30,30 "1"' \
        -fill purple -draw 'text 10,170 "head"' \
        -fill purple -draw 'text 210,170 "tail"' \
        -fill orange -draw 'text 420,30 "2"' \
        -fill purple -draw 'text 385,80 "tail"' \
        -fill purple -draw 'text 610,80 "head"' \
        -fill green -draw 'text 780,30 "3"' \
        -fill black -draw 'text 450,300 "β-cyclodextrin"' \
        -fill purple -draw 'text 960,610 "primary face"' \
        -fill purple -draw 'text 960,370 "secondary face"' \
        $DIR/dimers_BCD.png \
        $DIR/dimers_BCD_labelled.png