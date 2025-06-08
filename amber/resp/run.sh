for name in ETAM MGDM MIMM; do
    echo "Running $name"
    bash RESP2_noopt.sh $name.pdb 1
done

# restrained calculation
bash gas_gaussian.sh ACET.pdb -1
Multiwfn_noGUI gau.fchk

mv gau.chg gas.chg

bash solv_gaussian.sh ACET.pdb -1
Multiwfn_noGUI gau.fchk
mv gau.chg solv.chg

bash average.sh ACET.pdb