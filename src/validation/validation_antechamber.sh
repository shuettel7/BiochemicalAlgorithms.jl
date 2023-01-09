#!/bin/bash
# Attention: Paste the script below into terminal if operating on Windows with WSL, errors with '\r' due to dos vs unix 
for file in ../huettel-msc/FDA_PubChem_category_mol2_files/*.mol2; do
filename=$(basename "$file")
echo "$filename"
antechamber -i $file -fi mol2 -o ../huettel-msc/FDA_PubChem_category_mol2_files/antechamber_output/$filename -fo mol2 -at gaff -s 0
done
