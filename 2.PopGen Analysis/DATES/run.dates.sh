sed -r -e "s/CHG|Iran_GanjDareh_N|Iran_LN_SehGabi|Iran_TepeAbdulHosein_N.SG|Iran_Wezmeh_N.SG/Iran_Neo/g" data/MATRIX.ind > data/MATRIX_DATESgroup1.ind

sed -r -e "s/CHG|Iran_GanjDareh_N|Armenian.DG|Georgian.DG|Russia_Abkhasian.DG/Iran_Neo/g" data/MATRIX.ind > data/MATRIX_DATES_DATESgroup2.ind

sed -r -e "s/Barcin_N|Anatolia_TellKurdu_EC/Anatolia_Neo/g" data/MATRIX.ind > data/MATRIX_DATES_DATESgroup3.ind


for i in `seq 1 3` ; do dates -p par_dates_MenTesh_group$i > log_dates_MenTesh_group$i ; done
