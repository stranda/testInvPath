####
#### this is a series of function calls that create abc simulation environments for each
#### taxon.  The source material is in the gigas-popgen github repo
####

##Undaria
species_setup(root="~/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Undaria_pinnatifida/",fas="Uwai_141inds_2loc.fas",indmeta="Uwai_141inds.meta.csv",mname="upinn_meta_updated.csv",newdir="~/tmp/upinn",species="Undaria_pinnatifida",dataType="sequence",fsc_exec="fsc27")

## Battilaria
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Batillaria_attramentaria",genofile=NULL,fas="Battr_180inds.fas",indmeta="Battr_180inds_meta.csv",mname="battr_meta.csv",newdir="~/tmp/battr",species="Batillaria_attramentaria",dataType="sequence",fsc_exec="fsc27")

##Cercaria_batillaria_HL1
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Cercaria_batillaria_HL1",genofile=NULL,fas="HL1_230inds.fas",indmeta="HL1_230inds_meta.csv",mname="battrHL1_meta.csv",newdir="~/tmp/battrHL1",species="Cercaria_batillaria_HL1",dataType="sequence",fsc_exec="fsc27")

##Cercaria_batillaria_HL6
species_setup(root="/home/astrand/GoogleDrive/data/Oyster/gigas-popgen/ABC_input/Cercaria_batillaria_HL6",genofile=NULL,fas="HL6_428inds.fas",indmeta="HL6_428inds_meta.csv",mname="battrHL6_meta.csv",newdir="~/tmp/battrHL6",species="Cercaria_batillaria_HL6",dataType="sequence",fsc_exec="fsc27")

