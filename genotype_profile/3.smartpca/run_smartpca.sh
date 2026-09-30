convertf -p par.convertf.txt

awk 'NR==FNR {pop[$1]=$2; next} {print $1, "U", (pop[$1] ? pop[$1] : "TARGET")}' \
  integrated_call_samples_v3.20130502.ALL.panel \
  <(cut -d' ' -f1 merged.pruned.fam) \
  > merged.pruned.updated.ind

awk 'NR>1 {print $2}' integrated_call_samples_v3.20130502.ALL.panel \
  | sort -u > poplist.txt

smartpca -p par.smartpca.txt > smartpca.log
