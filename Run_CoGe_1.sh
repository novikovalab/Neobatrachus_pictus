# CoGe to circular synteny plot

# Set working directory (replace with your own CoGe analysis folder under home if needed)
cd ~/CoGe_analysis

# Step 1: Convert CoGe output into synteny blocks
#################################################

# CoGe Ks file downloaded from SynMap:
# "Results with synonymous/non-synonymous rate values" in the download section
file="~/CoGe_analysis/52485_59218.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.ks.txt"

# Remove CoGe header and comment lines
grep -v "##" "$file" > t
grep -v "#This" t > t2
mv t2 t

# Extract lines around "Ks" markers to define block starts and ends
grep -A 1 "Ks" t > t.starts
grep -B 1 "Ks" t > t.ends

# Clean up Ks lines for starts and ends
grep -v "#Ks" t.starts > x.starts
perl -i -p -e 's/-//g' x.starts

grep -v "#Ks" t.ends > x.ends
perl -i -p -e 's/-//g' x.ends

# Check lengths of start and end files
wc -l x.starts
wc -l x.ends

# Important: x.ends is typically one line shorter than x.starts,
# because grep -B 1 "#Ks" misses the last block end.
# We need to add the last line from the original CoGe file back to x.ends.
tail -1 "$file"

# Manually append the last line to x.ends or handle it in R (see below).
# After this, x.starts and x.ends should have the same number of lines.

# Next steps are done in R.
