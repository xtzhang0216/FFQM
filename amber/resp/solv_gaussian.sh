#!/bin/bash
delta=0.5
level_SP="B3LYP/6-311G**"
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


#### Single point in solvent
cat << EOF > gau.gjf
%nproc=32
%chk=gau.chk
$keyword_SP_solv geom=allcheck guess=read


EOF

echo
echo Running single point task in solvent phase via Gaussian...
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

