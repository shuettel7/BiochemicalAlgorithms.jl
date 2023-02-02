#!/bin/bash
# Attention: Paste the script below into terminal if operating on Windows with WSL, errors with '\r' due to dos vs unix 
for file in ../huettel-msc/gaff_files/*.mol2; do
filename="$(basename "$file" .mol2).mol2"
echo "$filename"
antechamber -i $file -fi mol2 -o ../huettel-msc/gaff_files/antechamber_output/$filename -fo mol2 -at gaff
done
