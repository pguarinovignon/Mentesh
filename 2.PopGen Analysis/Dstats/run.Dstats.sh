for i in {1..301..50} ; do echo "qpDstat -p Dstat_Mentesh_Paleo_Neo_select.par -L $i -H $(($i+49)) > log_Dstat_Mentesh_Paleo_Neo_select_$i" ; done | parallel

for i in {1..1501..50} ; do echo "qpDstat -p Dstat_MT_Caucase_Paleo_Neo_select.par -L $i -H $(($i+49)) > log_Dstat_MT_Caucase_Paleo_Neo_select_samples_$i" ; done | parallel -j 8

for i in {1..151..50} ; do echo "qpDstat -p Dstat_Mentesh_caucase_BA.par -L $i -H $(($i+49)) > log_Dstat_Mentesh_caucase_BA_$i" ; done | parallel

or i in {1..101..50} ; do echo "qpDstat -p Dstat_Mentesh_anat_BA.par -L $i -H $(($i+49)) > log_Dstat_Mentesh_anat_BA_$i" ; done | parallel
