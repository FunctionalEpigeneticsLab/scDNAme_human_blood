# prepare data
methscan prepare --input-format allc ./*.allc.tsv.gz compact_data

# filter out the .npz files that are not in the list of allowed chromosomes
cp -r compact_data filtered_data
find filtered_data -name '*.npz' -type f -printf '%f\n' | sed 's/\.npz\$//' | grep -vFf "$chrom_methscan" | xargs -I{} rm -f "filtered_data/{}.npz"

# smooth
methscan smooth filtered_data

# call VMRs
methscan scan --bridge-gaps 1000 filtered_data VMRs.autosome.bed

# build matrix
matrix=VMR_matrix
bedfile=VMRs.autosome.bed
methscan_mtx_input_dir=filtered_data
methscan matrix $bedfile $methscan_mtx_input_dir $matrix

# convert into anndata object
python /path/to/csv_to_mm_ad.py --input_file ${matrix}/mean_shrunken_residuals.csv.gz --outdir ${matrix}/matrices --prefix MSR --write_anndata --write_cells_and_features
python /path/to/csv_to_mm_ad.py --input_file ${matrix}/total_sites.csv.gz --outdir ${matrix}/matrices --prefix cov --write_cells_and_features
python /path/to/csv_to_mm_ad.py --input_file ${matrix}/methylated_sites.csv.gz --outdir ${matrix}/matrices --prefix mc --write_cells_and_features
