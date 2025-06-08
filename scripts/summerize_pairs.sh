folder=$(pwd)
for name in ACEM ACET MBZ MIMD MIME MIMM MIND ETAM ETOH ETSH MGDM MSM NMA MPHE PRPA N1PA HOH; do
    count=$(find "$folder" -name "*${name}*" | wc -l)
    echo "Total files for $name: $count"
done
