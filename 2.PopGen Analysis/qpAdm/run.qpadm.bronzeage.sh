echo "Mbuti.DG
Ami.DG
Mixe.DG
Russia_Kostenki14.SG
Russia_MA1_HG.SG
EEHG
Italy_North_Villabruna_HG
Natufian
Iran_GanjDareh_N
Barcin_N" > right_base_BA

echo "CHG
Russia_Caucasus_Eneolithic
Caucasus_lowlands_LC
Mentesh
Russia_Steppe_Eneolithic" > set_BA_left

########## RUN FOR 2 POPS COMBINAISONS ##########

python3 rotation_combination_BA_left.py 2 > combination_left_MT_BA_2pop
sed -i -e "s/(//g" combination_left_MT_BA_2pop
sed -i -e "s/'//g" combination_left_MT_BA_2pop
sed -i -e "s/)//g" combination_left_MT_BA_2pop
sed -i -e "s/,//g" combination_left_MT_BA_2pop

COMB2=combination_left_MT_BA_2pop
LENGTH2=`wc -l ${COMB2}|awk '{print$1}'`

Rscript create_left_right_files.R combination_left_MT_BA_2pop set_BA_left right_base_BA Armenia_C

echo "DIR: ./
S1: MATRIX
indivname: MATRIX_qpadm_samples.ind
snpname: MATRIX.snp
genotypename: MATRIX.geno
popleft: left_Armenia_C_set1_2pop_1.txt
popright: right_Armenia_C_set1_2pop_1.txt
details: YES
allsnps: YES
inbreed: NO" > qpAdm_rotating_Armenia_C_set1_2pop_1.par

for i in $(seq 2 "$LENGTH2") ; do  sed -e "s/left_Armenia_C_set1_2pop_1/left_Armenia_C_set1_2pop_$i/g" qpAdm_rotating_Armenia_C_set1_2pop_1.par > qpAdm_rotating_Armenia_C_set1_2pop_$i.par ; done
for i in $(seq 2 "$LENGTH2") ; do  sed -i -e "s/right_Armenia_C_set1_2pop_1/right_Armenia_C_set1_2pop_$i/g" qpAdm_rotating_Armenia_C_set1_2pop_$i.par ; done

echo "Armenia_C
Armenia_EBA
Armenia_LBA.SG
Armenia_MBA
Armenia_Caucasus_EBA_KuraAraxes
ARM_Sarukhan_LBA
ARM_Pijut_LBA
ARM_Noratus_LBA
ARM_Nerkin_Getashen_LBA
ARM_Lori_Berd_LBA
ARM_Lchashen_LBA
ARM_Keti_Urartian
ARM_Keti_LBA
ARM_Karashamb_LBA
ARM_Dzori_Gekh_LBA
ARM_Black_Fortress_LBA
ARM_Tekhut_LBA
ARM_Tavshut_Trialeti_MBA
ARM_Shengavit_KuraAraxes_EBA
ARM_Karnut_KuraAraxes_EBA
ARM_Berkaber_KuraAraxes_EBA" > list.target.BA

for k in `cat list.target.BA | tail -n +2` ; do 
  for i in $(seq 1 "$LENGTH2") ; do
    sed -e "s/Armenia_C/${k}/g" left_Armenia_C_set1_2pop_$i.txt > left_${k}_set1_2pop_$i.txt
    sed -e "s/left_Armenia_C/left_${k}/g" qpAdm_rotating_Armenia_C_set1_2pop_$i.par > qpAdm_rotating_{k}_set1_2pop_$i.par
  done
done

for k in `cat list.target` ; do 
  for i in$(seq 1 "$LENGTH2") ; do 
    echo "qpAdm -p qpAdm_rotating_${k}_set1_2pop_$i.par > log_qpAdm_rotating_${k}_set1_2pop_$i.par" ; 
   done ; 
  done | parallel -j 20
  
########## REDO FOR 3POP COMBINAISONs ##########

python3 rotation_combination_BA_left.py 3 > combination_left_MT_BA_3pop
sed -i -e "s/(//g" combination_left_MT_BA_3pop
sed -i -e "s/'//g" combination_left_MT_BA_3pop
sed -i -e "s/)//g" combination_left_MT_BA_3pop
sed -i -e "s/,//g" combination_left_MT_BA_3pop

COMB3=combination_left_MT_BA_3pop
LENGTH3=`wc -l ${COMB3}|awk '{print$1}'`

Rscript create_left_right_files.R combination_left_MT_BA_3pop set_BA_left right_base_BA Armenia_C

echo "DIR: ./
S1: MATRIX
indivname: MATRIX_qpadm_samples.ind
snpname: MATRIX.snp
genotypename: MATRIX.geno
popleft: left_MT23_set1_3pop_1.txt
popright: right_MT23_set1_3pop_1.txt
details: YES
allsnps: YES
inbreed: NO" > qpAdm_rotating_MT23_set1_3pop_1.par

for i in $(seq 2 "$LENGTH3") ; do  sed -e "s/left_MT23_set1_3pop_1/left_MT23_set1_3pop_$i/g" qpAdm_rotating_MT23_set1_3pop_1.par > qpAdm_rotating_MT23_set1_3pop_$i.par ; done
for i in $(seq 2 "$LENGTH3") ; do  sed -i -e "s/right_MT23_set1_3pop_1/right_MT23_set1_3pop_$i/g" qpAdm_rotating_MT23_set1_3pop_$i.par ; done

for k in `cat list.target | tail -n 7` ; do 
  for i in $(seq 1 "$LENGTH3") ; do
    sed -e "s/MT23/${k}/g" left_MT23_set1_3pop_$i.txt > left_${k}_set1_3pop_$i.txt
    sed -e "s/left_MT23/left_${k}/g" qpAdm_rotating_MT23_set1_3pop_$i.par > qpAdm_rotating_{k}_set1_3pop_$i.par
  done
done

for k in `cat list.target` ; do 
  for i in $(seq 1 "$LENGTH3") ; do 
    echo "qpAdm -p qpAdm_rotating_${k}_set1_3pop_$i.par > log_qpAdm_rotating_${k}_set1_3pop_$i.par" ; 
   done ; 
  done | parallel -j 20
  
######### PARSE THE DATA ##########

COMB2=${DIR}combination_left_MT_BA_2pop
LENGTH2=`wc -l ${COMB2}|awk '{print$1}'`
COMB3=${DIR}combination_left_MT_BA_3pop
LENGTH3=`wc -l ${COMB3}|awk '{print$1}'`
LIST=${DIR}list.target.BA

echo "#ID POP1 POP2 p coeff_POP1 coeff_POP2 stderr1 stderr2 nested_model nested_model_p " > ${DIR}results.qpAdm_BA_2pop.txt
for NAME in `cat ${LIST}`; do
    for i in $(seq 1 "$LENGTH2") ; do
        POP1=`sed -n ${i}p $COMB2 | awk '{print $1}'`
        POP2=`sed -n ${i}p $COMB2 | awk '{print $2}'`
        echo -en $NAME"\t"
        echo "${POP1} ${POP2}" | tr "\n" "\t"
        grep summ ${DIR}log_qpAdm_rotating_${NAME}_set1_2pop_${i}.par | awk '{print $4" " $5" " $6}' | tr "\n" "\t"
        grep std ${DIR}log_qpAdm_rotating_${NAME}_set1_2pop_${i}.par| awk '{print $3" " $4}' | tr "\n" "\t"
        grep nested ${DIR}log_qpAdm_rotating_${NAME}_set1_2pop_${i}.par| head -n 1| awk '{print $3" " $11}' | tr "\n" "\t"
        echo "";
    done |sort -k4 >> ${DIR}results.qpAdm_BA_2pop.txt ;
done	

echo "#ID POP1 POP2 POP3 p coeff_POP1 coeff_POP2 coeff_POP3 stderr1 stderr2 stderr3 nested_model nested_model_p " > ${DIR}results.qpAdm_BA_3pop.txt
for NAME in `cat ${LIST}`; do
for i in $(seq 1 "$LENGTH3") ; do
POP1=`sed -n ${i}p $COMB3 | awk '{print $1}'`
POP2=`sed -n ${i}p $COMB3 | awk '{print $2}'`
POP3=`sed -n ${i}p $COMB3 | awk '{print $3}'`
echo -en $NAME"\t"
echo "${POP1} ${POP2} ${POP3}" | tr "\n" "\t"
grep summ ${DIR}log_qpAdm_rotating_${NAME}_set1_3pop_${i}.par | awk '{print $4" " $5" " $6" "$7}' | tr "\n" "\t"
grep std ${DIR}log_qpAdm_rotating_${NAME}_set1_3pop_${i}.par| awk '{print $3" " $4" "$5}' | tr "\n" "\t"
grep nested ${DIR}log_qpAdm_rotating_${NAME}_set1_3pop_${i}.par| head -n 1| awk '{print $3" " $11}' | tr "\n" "\t"
echo "";
done |sort -k4 >> ${DIR}results.qpAdm_BA_3pop.txt ;
done	
