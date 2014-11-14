#!/bin/bash

# Script to help automate the extraction of data from the Hobbs 2010 paper
# This wont run everything - as some manual labour is required it's here
# to remind me what to do in future

# First use pdfshuffler to save only the 4 pages of data called
source=HobbsLyneKramer_REDUCED.pdf

# Turn pdf into ppm
if [[ $@ == *ppm* ]]
then
    pdfimages ${source} images
fi

# Open ppm files in gimp and "clean" save as cleaned--

# Turn them into text files
if [[ $@ == *txt* ]]
then
    for file in cleaned*.ppm 
    do
    gocr -i $file -C "0-9.+-" -o ${file%.ppm}.txt -s 1000
    cat ${file%.ppm}.txt |tr -s " " > stripped-${file%.ppm}.txt
    done
fi

# Convert all stripped-cleaner*txt files to nice space separated file
if [[ $@ == *dat* ]]
then
    python CleanedTxtToDat.py
fi

# You now have a nice database of the data!

