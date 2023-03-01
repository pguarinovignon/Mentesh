echo "genotypename:    1240v2_HO_NP_final_MDH_laza2022.geno
snpname:         1240v2_HO_NP_final_MDH_laza2022.snp
indivname:       1240v2_HO_NP_final_MDH_laza2022.ind
outputformat:    PED
genotypeoutname: 1240v2_HO_NP_final_MDH_laza2022.ped
snpoutname:      1240v2_HO_NP_final_MDH_laza2022.map
indivoutname:    1240v2_HO_NP_final_MDH_laza2022.pedind
outputgroup : YES" > convertf_ped.txt

./convertf -p convertf_ped.txt > log_convertf_ped.txt

Rscript subset_for_admixture.r

plink1.9 --file 1240v2_HO_NP_final_MDH_laza2022 --keep tokeep_admixture.txt --recode --out MATRIX_adm
plink1.9 --file MATRIX_adm --indep-pairwise 200 25 0.4
plink1.9 --file MATRIX_adm --extract plink.prune.in --make-bed --out MATRIX_adm_pruned --allow-no-sex

for i in {01..20}; do
        cd {root}/ADMIXTURE/rep$i
        for k in {2..14}; do
                echo "admixture --cv -s time /MATRIX_adm_pruned.bed $k -j10 > log_Eurasie_NP_ULGSK_K${k}_rep${i}.log";
        done | parallel ;
        for k in {2..14}; do
                mv MATRIX_adm_pruned.$k.P pruned_Eurasie_HO_NP_laza2022.$k.rep$i.P;
                mv MATRIX_adm_pruned.$k.Q pruned_Eurasie_HO_NP_laza2022.$k.rep$i.Q;
        done ;
done
