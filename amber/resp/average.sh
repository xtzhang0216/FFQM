#### Calculate RESP2
export inname=$1
filename=${inname%.*}
suffix=${inname##*.}
chgname=${1//$suffix/chg}
delta=0.5

paste gas.chg solv.chg |awk '{printf ("%-3s %12.6f %12.6f %12.6f %15.10f\n",$1,$2,$3,$4,(1-d)*$5+d*$10)}' d=$delta > $chgname
