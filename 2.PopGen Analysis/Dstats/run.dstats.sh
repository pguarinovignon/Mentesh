 cat list_select_MT23 > list_select_MT_samples
 sed 's/MT23/MT26/g' list_select_MT23 >> list_select_MT_samples
 sed 's/MT23/MT7/g' list_select_MT23 >> list_select_MT_samples
 sed 's/MT23/Mentesh_001/g' list_select_MT23 >> list_select_MT_samples
 sed 's/MT23/Mentesh/g' list_select_MT23 > list_select_Mentesh


for i in {1..301..50} ; do echo "qpDstat -p Dstat_Mentesh_Paleo_Neo_select.par -L $i -H $(($i+49)) > log_Dstat_Mentesh_Paleo_Neo_select_$i" ; done | parallel

for i in {1..1501..50} ; do echo "qpDstat -p Dstat_MT_Caucase_Paleo_Neo_select.par -L $i -H $(($i+49)) > log_Dstat_MT_Caucase_Paleo_Neo_select_$i" ; done | parallel

grep result log_Dstat_Mentesh_Paleo_Neo_select_* > out_Dstat_Mentesh_Paleo_Neo_select
grep result log_Dstat_MT_Paleo_Neo_select_* > out_Dstat_MT_Paleo_Neo_select

for i in {1..151..50} ; do echo "qpDstat -p Dstat_Mentesh_caucase_BA.par -L $i -H $(($i+49)) > log_Dstat_Mentesh_caucase_BA_$i" ; done | parallel
for i in {1..101..50} ; do echo "qpDstat -p Dstat_Mentesh_anat_BA.par -L $i -H $(($i+49)) > log_Dstat_Mentesh_anat_BA_$i" ; done | parallel

grep result log_Dstat_Mentesh_caucase_BA_* > out_Dstat_Mentesh_caucase_BA
grep result log_Dstat_Mentesh_anat_BA_* > out_Dstat_Mentesh_anat_BA
