###
### a function that creates a list of parameters for the simulations
###  MUST be edited for each species
###

#' creates the species and problem-level parameters
#' @param cores number of cores (usually set to 1)
#' @param nreps number of reps to run per computation
#' @param species name of species under consideration
#' @param demeN controls ghost pops in native regions
#' @param demeI controls ghost pops in introduced regions
#' @param sources the possible source regions for this species (must have a sample in each)
#' @param intros the introduced regions (must have a sample in each of these regions)
#' @param exclude source populations to exclude from analysis (ususally nonsource)
#' @param dataType the type of genetic data (found in datafile.R)
#' @param fsc_exec name of the fastsimcoal executable
#'
#' @description The sources and intros are calculated when you run the
#'     'species_setup()' function and when the simulation hierarchy is
#'     created, these are stored in a file called datafiles.R in the
#'     'src'subdirectory.  This function is really only used in the call to setupReps()
#' @return list containing parameters
#' @export
mainparams <-function(cores=1, 
                      nreps=1, 
                      species="",
                      demeN = 0,  #number of populations in each native region (demeN - #sampled) pops are ghosts in each native region)
                      demeI = 0,  #number of populations in each intro region (demeI - #sampled) pops are ghosts in each intro region)
                      sources=c("kag","hok","tok","hon","sea"), #there should be samples from everything in this vector
                      intros=c("PNW","Cali","NZ","Arg","EU"), #there should be samples from everything in this vector
                      dataType="",
                      fsc_exec="",
                      exclude="nonSource",
                      popPairwise=popPairwise
                      )
{
    list(cores=cores, #cores 
         nreps=nreps, #numerof replicated sims
         species=species,
         demeN = demeN,  #number of populations in each native region (demeN - #sampled) pops are ghosts in each native region)
         demeI = demeI,  #number of populations in each intro region (demeI - #sampled) pops are ghosts in each intro region)
         sources=sources, #there should be samples from everything in this vector
         intros=intros, #there should be samples from everything in this vector
         dataType=dataType,
         fsc_exec=fsc_exec,
         exclude=exclude,
         popPairwise=popPairwise
         )
}
