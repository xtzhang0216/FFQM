for name in ACEM ACET MBZ MIMD MIME MIMM MIND ETAM ETOH ETSH MGDM MSM NMA MPHE PRPA N1PA HOH; do
    echo "Running $name"
    g16 < geometry_opt/gaussian_opt/${name}.gjf > geometry_opt/gaussian_results/${name}.log
done