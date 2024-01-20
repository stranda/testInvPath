####
#### this is a series of function calls that create abc simulation environments for each
#### taxon.  The source material is in the gigas-popgen github repo
####
library(testInvPath)



##Undaria
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Undaria_pinnatifida/",fas="Uwai_141inds_2loc.fas",indmeta="Uwai_141inds.meta.csv",mname="upinn_meta_updated.csv",newdir="/home/astrand/tmp/upinn",species="Undaria_pinnatifida",dataType="sequence",fsc_exec="fsc27")

## Battilaria
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Batillaria_attramentaria",genofile=NULL,fas="Battr_180inds.fas",indmeta="Battr_180inds_meta.csv",mname="battr_meta.csv",newdir="/home/astrand/tmp/battr",species="Batillaria_attramentaria",dataType="sequence",fsc_exec="fsc27")

##Cercaria_batillaria_HL1
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Cercaria_batillaria_HL1",genofile=NULL,fas="HL1_230inds.fas",indmeta="HL1_230inds_meta.csv",mname="battrHL1_meta.csv",newdir="/home/astrand/tmp/battrHL1",species="Cercaria_batillaria_HL1",dataType="sequence",fsc_exec="fsc27")

##Cercaria_batillaria_HL6
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Cercaria_batillaria_HL6",genofile=NULL,fas="HL6_428inds.fas",indmeta="HL6_428inds_meta.csv",mname="battrHL6_meta.csv",newdir="/home/astrand/tmp/battrHL6",species="Cercaria_batillaria_HL6",dataType="sequence",fsc_exec="fsc27")

##Haminoea_japonica
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Haminoea_japonica/",fas="Haminoea_140ind_aligned.fas",indmeta="Table1_IndAcc.csv",mname="hami_meta.csv",newdir="/home/astrand/tmp/hami",species="Haminoea_japonica",dataType="sequence",fsc_exec="fsc27")

##Ulva_pertusa
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Ulva_pertusa/",fas="Ulva_118ind.fas",indmeta="Ulva_118ind_meta.csv",mname="ulva_meta.csv",newdir="/home/astrand/tmp/ulva",species="Ulva_pertusa",dataType="sequence",fsc_exec="fsc27")


##Mutimo_cylindricus 
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Mutimo_cylindricus/",fas="Mutimo_91ind.fas",indmeta="Mutimo_91ind_meta.csv",mname="mutimo_meta.csv",newdir="/home/astrand/tmp/muti",species="Mutimo_cylindricus",dataType="sequence",fsc_exec="fsc27")

##Acanthogobius_flavimanus
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Acanthogobius_flavimanus/",fas="Acanthogobius_flavimanus_605ind.fas",indmeta="Acanthogobius_flavimanus_605ind_meta.csv",mname="acanthogobius_meta.csv",newdir="/home/astrand/tmp/goby",species="Acanthogobius_flavimanus",dataType="sequence",fsc_exec="fsc27")

##Didemnum_vexillum
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Didemnum_vexillum/",fas="Dide_189inds.fas",indmeta="Dide_189inds_meta.csv",mname="dide_meta.csv",newdir="/home/astrand/tmp/dide",species="Didemnum_vexillum",dataType="sequence",fsc_exec="fsc27")


##Hemigrapsus_takanoi
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Hemigrapsus_takanoi/",fas="Hemigrapsus_takanoi_624ind_aligned.fas",indmeta="Hemigrapsus_takanoi_624ind_meta.csv",mname="htakanoi_meta.csv",newdir="/home/astrand/tmp/htakanoi",species="Hemigrapsus_takanoi",dataType="sequence",fsc_exec="fsc27")

##Hemigrapsus_sanguineus
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Hemigrapsus_sanguineus/",fas="Hemigrapsus_sang_413ind.fas",indmeta="Hemigrapsus_sang_413ind_meta.csv",mname="hsanguineus_meta.csv",newdir="/home/astrand/tmp/hsanguineus",species="Hemigrapsus_sanguineus",dataType="sequence",fsc_exec="fsc27")

##Gracillaria_vermiculophylla_SSR

species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Gracilaria_vermiculophylla_SSR/",genofile="gverm_usat_byPop_gtype",mname="gverm_meta.csv",newdir="/opt/data1/oyster/introduced/gvermSSR",species="Gracilaria_vermiculophylla_SSR",dataType="microsatellite",fsc_exec="fsc27",popPairwise=TRUE)


##Gracillaria_vermiculophylla_SNP

species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Gracilaria_vermiculophylla_SNP/",genofile="gverm_snps_229ind_1000loci.RDS",mname="meta-SNP.pops_edited.V2.csv",newdir="/home/astrand/tmp/gvermSNP",species="Gracilaria_vermiculophylla_SNP",dataType="snp",fsc_exec="fsc27",nativeTopology= matrix(c(0,0,0,0,  #kag
                                      0,0,0,1,  #hon
                                      2,3,0,0,  #tok
                                      0,0,0,0), #hok
                                          byrow=T,nrow=4,
                                          dimnames=list(c("kag","hon","tok","hok"),
                                                        c("kag","hon","tok","hok"))
                                          ))

##Gracillaria_vermiculophylla_mtDNA
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Gracilaria_vermiculophylla_mtDNA/",fas="gvermMTDNA.606ind.fas",indmeta="gvermMTDNA.606ind.csv",mname="gverm_meta.csv",newdir="/opt/data1/oyster/introduced/gvermMTDNA",species="Gracilaria_vermiculophylla_mtDNA",dataType="sequence",fsc_exec="fsc27",popPairwise=TRUE)


##

species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Palaemon_macrodactylus/",fas="Palaemon_238inds.fas",indmeta="Palaemon_238inds_meta.csv",mname="pala_meta.csv",newdir="/home/astrand/tmp/pala",species="Palaemon_macrodactylus",dataType="sequence",fsc_exec="fsc27")

##

species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Polydora_hoplura/",fas="Polydora_final.fas",indmeta="Polydora_inds_final.csv",mname="Polydora_meta.csv",newdir="/home/astrand/tmp/polyd",species="Polydora_hoplura",dataType="sequence",fsc_exec="fsc27")

species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Test_species/",fas="Test_final.fas",indmeta="Test_inds_final.csv",mname="Test_meta.csv",newdir="/home/astrand/tmp/testsp",species="Test_species",dataType="sequence",fsc_exec="fsc27",popPairwise=T)



