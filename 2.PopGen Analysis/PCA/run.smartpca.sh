echo "genotypename: ../intermediaire_merge/1240v2_HO_NP_final_MDH_laza2022.geno
snpname: ../intermediaire_merge/1240v2_HO_NP_final_MDH_laza2022.snp
indivname: ../intermediaire_merge/1240v2_HO_NP_final_MDH_laza2022.ind
evecoutname:1240v2_HO_NP_final_MDH_laza2022_eurasie_focus_caucase.evec
evaloutname: 1240v2_HO_NP_final_MDH_laza2022_eurasie_focus_caucase.eval
poplistname: toinclude_eurasie_focus_caucase_noEA_HO
lsqproject: YES
outliermode: 2" > par_pca_eurasie.txt

./smartpca -p par_pca_eurasie.txt > log_pca_eurasie.txt
