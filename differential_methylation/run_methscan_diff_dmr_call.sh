set -euo pipefail

wd=`pwd`
outdir=${wd}/methscan_diff_outs
mkdir -p $outdir

for cell_groups_file in ${wd}/methscan_diff_input/cell_groups.*.csv; do

	methscan_dir=methscan
	datdir=${methscan_dir}/methscan_outs/filtered_data
	base=$(basename "$cell_groups_file")
	comparison="${base#cell_groups.}"
	comparison="${comparison%.csv}"
	outfile="${outdir}/methscan.dmr.${comparison}.bed"

	echo "Running methscan diff for $comparison"
	methscan diff --threads ${SLURM_CPUS_PER_TASK} "$datdir" "$cell_groups_file" "$outfile"

done

