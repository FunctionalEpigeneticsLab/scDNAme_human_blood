#!/bin/bash

file1=$1
file2=$2
output_prefix=$3

sort -k1,1 -k2,2n "$file1" > "${file1}.sorted"
sort -k1,1 -k2,2n "$file2" > "${file2}.sorted"

# 1. Regions in file1 that OVERLAP with file2
bedtools intersect -a "${file1}.sorted" -b "${file2}.sorted" -u > "${output_prefix}.file1_overlapping.bed"

# 2. Get which regions in file2 overlap with file1
bedtools intersect -a "${file2}.sorted" -b "${file1}.sorted" -u > "${output_prefix}.file2_overlapping.bed"

# 3. Regions in file1 that do NOT overlap with file2
bedtools intersect -a "${file1}.sorted" -b "${file2}.sorted" -v > "${output_prefix}.file1_unique.bed"

# 4. Regions in file2 that do NOT overlap with file1
bedtools intersect -a "${file2}.sorted" -b "${file1}.sorted" -v > "${output_prefix}.file2_unique.bed"

# 5. Get the actual overlap regions (intersection coordinates)
bedtools intersect -a "${file1}.sorted" -b "${file2}.sorted" > "${output_prefix}.overlapping_regions.bed"

rm "${file1}.sorted" "${file2}.sorted"
