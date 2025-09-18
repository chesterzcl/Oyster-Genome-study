#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=15G
#SBATCH --mail-type=ALL


# Input depth file (3 columns: scaffold, position, depth)
DEPTH_FILE=/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/wgs_depth.txt

if [ -z "$DEPTH_FILE" ]; then
  echo "Usage: $0 <depth.txt>"
  exit 1
fi

OUTDIR="/home/zl436/palmer_scratch/oyster_genome/Contigs_dedup_scfd_3ddna_mc/wgs_depth_summary"
mkdir -p $OUTDIR

# Output files
SUMMARY_FILE="${OUTDIR}/depth_summary_by_scaffold.tsv"
DIST_FILE="${OUTDIR}/depth_distribution_by_scaffold.tsv"

echo "Scaffold	Mean	Median	Min	Max	CoveredBases" > $SUMMARY_FILE

# Loop through each scaffold and compute stats
awk '
{
  scaffold=$1; depth=$3;
  depths[scaffold][++count[scaffold]]=depth;
  sum[scaffold]+=depth;
  if (depth > max[scaffold] || count[scaffold]==1) max[scaffold]=depth;
  if (depth < min[scaffold] || count[scaffold]==1) min[scaffold]=depth;
}
END {
  PROCINFO["sorted_in"] = "@ind_str_asc"
  for (sc in depths) {
    n = count[sc];
    split("", sorted);
    for (i = 1; i <= n; i++) sorted[i] = depths[sc][i];
    asort(sorted);
    if (n % 2 == 1) {
      median = sorted[(n + 1) / 2];
    } else {
      median = (sorted[n / 2] + sorted[n / 2 + 1]) / 2;
    }
    mean = sum[sc] / n;
    printf "%s\t%.2f\t%.2f\t%d\t%d\t%d\n", sc, mean, median, min[sc], max[sc], n >> "'$SUMMARY_FILE'"
  }
}' $DEPTH_FILE

# Generate depth distribution
awk '{count[$1"\t"$3]++} END {
  print "Scaffold\tDepth\tCount";
  for (i in count) print i "\t" count[i];
}' $DEPTH_FILE > $DIST_FILE

echo "✔ Depth summary written to: $SUMMARY_FILE"
echo "✔ Depth distribution written to: $DIST_FILE"

