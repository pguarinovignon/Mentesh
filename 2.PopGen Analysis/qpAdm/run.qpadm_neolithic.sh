echo "Mbuti.DG
Ami.DG
Mixe.DG
Russia_Kostenki14.SG
Russia_MA1_HG.SG
EEHG
Italy_North_Villabruna_HG
Natufian" > right_base

echo "Iran_GanjDareh_N
CHG
PPN
Barcin_N
Anatolia_TellKurdu_EC
IRQ_Nemrik9_PPN
TUR_SE_Mardin_PPN" > set_1_left

########## RUN FOR 2 POPS COMBINAISONS ##########

python3 rotation_combination_Neo_left.py 2 > combination_left_MT_2pop
sed -i -e "s/(//g" combination_left_MT_2pop
sed -i -e "s/'//g" combination_left_MT_2pop
sed -i -e "s/)//g" combination_left_MT_2pop
sed -i -e "s/,//g" combination_left_MT_2pop

N=`wc -l combination_left_MT_2pop`

Rscript create_left_right_files.R

echo "DIR: ./
S1: MATRIX
indivname: MATRIX_qpadm_samples.ind
snpname: MATRIX.snp
genotypename: MATRIX.geno
popleft: left_MT23_set1_2pop_1.txt
popright: right_MT23_set1_2pop_1.txt
details: YES
allsnps: YES
inbreed: NO" > qpAdm_rotating_MT23_set1_2pop_1.par

for i in `seq 2 $N` ; do  sed -e "s/left_MT23_set1_2pop_1/left_MT23_set1_2pop_$i/g" qpAdm_rotating_MT23_set1_2pop_1.par > qpAdm_rotating_MT23_set1_2pop_$i.par ; done
for i in `seq 2 $N` ; do  sed -i -e "s/right_MT23_set1_2pop_1/right_MT23_set1_2pop_$i/g" qpAdm_rotating_MT23_set1_2pop_$i.par ; done

echo "MT23
MT26
MT7
MT1
PLT2
ARM_Akn
ARM_MB
Mentesh" > list.target
for k in `cat list.target | tail -n 7` ; do 
  for i in `seq 1 21` ; do
    sed -e "s/MT23/${k}/g" left_MT23_set1_2pop_$i.txt > left_${k}_set1_2pop_$i.txt
    sed -e "s/left_MT23/left_${k}/g" qpAdm_rotating_MT23_set1_2pop_$i.par > qpAdm_rotating_{k}_set1_2pop_$i.par
  done
done

for k in `cat list.target` ; do 
  for i in `seq 1 21` ; do 
    echo "qpAdm -p qpAdm_rotating_${k}_set1_2pop_$i.par > log_qpAdm_rotating_${k}_set1_2pop_$i.par" ; 
   done ; 
  done | parallel -j 20
  
########## REDO FOR 3POP COMBINAISONs ##########

python3 rotation_combination_Neo_left.py 3 > combination_left_MT_3pop
sed -i -e "s/(//g" combination_left_MT_3pop
sed -i -e "s/'//g" combination_left_MT_3pop
sed -i -e "s/)//g" combination_left_MT_3pop
sed -i -e "s/,//g" combination_left_MT_3pop

N=`wc -l combination_left_MT_3pop`

Rscript create_left_right_files_3pop.R

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

for i in `seq 2 $N` ; do  sed -e "s/left_MT23_set1_3pop_1/left_MT23_set1_3pop_$i/g" qpAdm_rotating_MT23_set1_3pop_1.par > qpAdm_rotating_MT23_set1_3pop_$i.par ; done
for i in `seq 2 $N` ; do  sed -i -e "s/right_MT23_set1_3pop_1/right_MT23_set1_3pop_$i/g" qpAdm_rotating_MT23_set1_3pop_$i.par ; done

for k in `cat list.target | tail -n 7` ; do 
  for i in `seq 1 21` ; do
    sed -e "s/MT23/${k}/g" left_MT23_set1_3pop_$i.txt > left_${k}_set1_3pop_$i.txt
    sed -e "s/left_MT23/left_${k}/g" qpAdm_rotating_MT23_set1_3pop_$i.par > qpAdm_rotating_{k}_set1_3pop_$i.par
  done
done

for k in `cat list.target` ; do 
  for i in `seq 1 21` ; do 
    echo "qpAdm -p qpAdm_rotating_${k}_set1_3pop_$i.par > log_qpAdm_rotating_${k}_set1_3pop_$i.par" ; 
   done ; 
  done | parallel -j 20
  
######### PARSE THE DATA ##########

COMB2=${DIR}combination_left_MT_2pop
LENGTH2=`wc -l ${COMB2}|awk '{print$1}'`
COMB3=${DIR}combination_left_MT_3pop
LENGTH3=`wc -l ${COMB3}|awk '{print$1}'`
LIST=${DIR}list.target
echo "#ID POP1 POP2 p coeff_POP1 coeff_POP2 stderr1 stderr2 nested_model nested_model_p " > ${DIR}results.qpAdm_2pop.txt
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
    done |sort -k4 >> ${DIR}results.qpAdm_2pop.txt ;
done	

echo "#ID POP1 POP2 POP3 p coeff_POP1 coeff_POP2 coeff_POP3 stderr1 stderr2 stderr3 nested_model nested_model_p " > ${DIR}results.qpAdm_3pop.txt
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
done |sort -k4 >> ${DIR}results.qpAdm_3pop.txt ;
done	
