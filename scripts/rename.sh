for ff in charmm amber opls amoeba drude;
do 
cd $ff/data
# for file in *ACEM*.csv; do mv "$file" "$file.bak"; done
# for file in *MIME*.csv; do mv "$file" "$file.bak"; done
# for file in *MIMM*.csv; do mv "$file" "$file.bak"; done
for file in *ACEM*.csv; do rm "$file"; done
for file in *MIME*.csv; do rm "$file"; done
for file in *MIMM*.csv; do rm "$file"; done
for file in *MIMD*.csv; do rm "$file"; done
cd ../../
done