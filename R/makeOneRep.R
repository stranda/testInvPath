#' Create a single rep of a coalescent simulation based on pulls from prior distributions
#'
#' @param params these are the parameters produced by setupReps()
#' @param priors pulls from prior distribution to parameterize this
#'     rep (from getPriors())
#' @param demeN the total number of demes to be simulated in the
#'     native (source range)
#' @param demeI the total number of demes to be simulated in the
#'     introduced (source range)
#' @description this function currently uses strataG to parameterize a
#'     fastsimcoal model that depends on the structure determined in
#'     params (primarily source population, native topology,
#'     introduced populations).  the priors determine choices of
#'     timing, actual introduction model and rates of mutation and
#'     migration.  demeN and demeI try to simulate ghost populations
#'     that are ultimately not sampled.  If the values for these
#'     variables are less that the number of sampled demes in each
#'     region then just the sampled populations are simulated.
#'     Otherwise additional unsampled populations are simulated up to
#'     the number of demeN or demeI for each of the native or
#'     introduced regions, respectively
#' @return a fscWrite() object that can be used in strataG's fscRun() function
