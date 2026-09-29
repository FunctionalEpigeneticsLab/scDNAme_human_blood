folder=pseudobulk 
bams=(${folder}/*.bam)
fa=genome.fa

# freebayes
for input_bam in "${bams[@]}"; do

tmp_vcf=$(basename "$input_bam" .sorted.bam).unfiltered.vcf.gz
freebayes-parallel <(python fasta_generate_regions.py ${fa}.fai 100000) 20 \
    -f $fa --min-coverage 3 --min-base-quality 30 $input_bam | bgzip > $tmp_vcf

done

# reheader
ls -1v ${sample}*chr*.vcf.gz > "$list"

while read file; do
    bcftools index "$file"
done < $list

bcftools concat -f $list -O z -o "${sample}.unfiltered.vcf.gz"
bcftools index "${sample}.unfiltered.vcf.gz"

echo ${sample} > ${sample}.txt
bcftools reheader -s ${sample}.txt "${sample}.unfiltered.vcf.gz" -o "${sample}.unfiltered.reheadered.vcf.gz"
bcftools index "${sample}.unfiltered.reheadered.vcf.gz" && rm ${sample}.txt

# bcftools filter vcf
for f in *.unfiltered.reheadered.vcf.gz; do
    prefix=${f%.unfiltered.reheadered.vcf.gz}
    bcftools view -i 'QUAL>=20 && GT!="0/0" && FORMAT/DP>=3' "$f" -Oz -o "${prefix}.filtered_AF0_QUAL20_DP3.vcf.gz"
done
