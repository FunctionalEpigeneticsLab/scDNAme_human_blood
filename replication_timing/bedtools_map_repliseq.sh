#!/bin/bash

dmr_bed=$1
dmr_bed2=$2
repliseq_bed=$3
output_dir=${4:-.}

mkdir -p "$output_dir"

dmr_prefix=$(basename "$dmr_bed" .bed)
prefix="${output_dir}/${dmr_prefix}"

sort -k1,1 -k2,2n "$dmr_bed" > "${prefix}.sorted.bed"
sort -k1,1 -k2,2n "$dmr_bed2" > "${dmr_bed2}.sorted"
sort -k1,1 -k2,2n "$repliseq_bed" > "${repliseq_bed}.sorted"

# Map replication timing to DMRs
bedtools map \
  -a "${prefix}.sorted.bed" \
  -b "${repliseq_bed}.sorted" \
  -c 4 \
  -o mean \
  > "${prefix}.replication_timing.bed"

# Get background
bedtools subtract \
  -a "${repliseq_bed}.sorted" \
  -b "${prefix}.sorted.bed" | \
  bedtools subtract \
  -a - \
  -b "${dmr_bed2}.sorted" | \
  sort -k1,1 -k2,2n \
  > "${prefix}.replication_timing_background.sorted.bed"

# Map replication timing to background regions
bedtools map \
  -a "${prefix}.replication_timing_background.sorted.bed" \
  -b "${repliseq_bed}.sorted" \
  -c 4 \
  -o mean \
  > "${prefix}.replication_timing_background.bed"

# cleanup
rm "${prefix}.sorted.bed" "${dmr_bed2}.sorted" "${repliseq_bed}.sorted" "${prefix}.replication_timing_background.sorted.bed"
