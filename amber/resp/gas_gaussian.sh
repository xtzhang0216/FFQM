#!/bin/bash
delta=0.5
level_SP="B3LYP/6-311+G**"
Gaussian=g16

export inname=$1
filename=${inname%.*}
suffix=${inname##*.}

if [ $2 ];then
	echo "Net charge = $2"
	chg=$2
else
	echo "Net charge was not defined. Default to 0"
	chg=0
fi

if [ $3 ];then
	echo "Spin multiplicity = $3"
	multi=$3
else
	echo "Spin multiplicity was not defined. Default to 1"
	multi=1
fi

if [ $4 ];then
	echo Solvent is $4
	solvent="scrf(solvent="$4")"
else
	solvent="scrf(solvent=water)"
	echo "Solvent name was not defined. Default to water"
fi

echo delta parameter is $delta

keyword_SP_gas="# "$level_SP" pop=MK IOp(6/33=2,6/42=6)"
keyword_SP_solv="# "$level_SP" "$solvent" pop=MK IOp(6/33=2,6/42=6)"

#### Convert input file to .xyz file
Multiwfn_noGUI  $1 > /dev/null << EOF
100
2
2
tmp.xyz
0
q
EOF


#### Single point in gas
cat << EOF > gau.gjf
%nproc=32
%chk=gau.chk
$keyword_SP_gas

test

$chg $multi
EOF
awk '{if (NR>2) print }' tmp.xyz >> gau.gjf
cat << EOF >> gau.gjf


EOF
rm tmp.xyz

echo
echo Running single point task in gas phase via Gaussian...
$Gaussian < gau.gjf > gau.out

if grep -Fq "Normal termination" gau.out
then
	echo Done!
else
	echo The task has failed! Exit the script...
	exit 1
fi
echo Running formchk...
formchk gau.chk > /dev/null