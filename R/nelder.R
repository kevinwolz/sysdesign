#' Create a Nelder Fan experimental design
#' @description Creates a Nelder Fan experimental design.
#' @details The Nelder Fan or Nelder Wheel experimental design (Nelder 1962) is an experimental design that
#' systematically varies plant desnity within a single plot. This function creates Nelder Fan Type Ia (Nelder 1962),
#' where the growing area around each plant has a constant shape throughout the design but increases as radius increases.
#' The terminology and calculations used here follow Parrott, Brinks, and Lhotka (2012).
#' @return An object of class sysd and class nelder. This is a list of 3 elements:
#' \itemize{
#'  \item{"plants"}{ - A data frame (tibble) containing one row for each for each plant in the design.}
#'  \item{"plot"}{ - A data frame (tibble) containing plot charateristics.}
#'  \item{"even.optim"}{ - If \code{even = TRUE}, a data frame (tibble) comparing plot characteristics between
#'  the scenario where theta is even throughout the design and the scenario where theta is uneven.
#'  If \code{even = FALSE}, then \code{FALSE}.}
#' }
#' @param DN Plant density within the last experimental arc (plants ha-1)
#' (i.e. lower extreme of experimental plant density range).
#' @param D1 Plant density within the first experimental arc (plants ha-1)
#' (i.e. upper extreme of experimental plant density range).
#' @param N Number of experimental arcs (i.e. number of densities to be tested within D1 to DN).
#' @param tau The "rectangularity" proportion. "Rectangularity" is the proportional relationship between the arc length
#' between spokes and the radial length between arcs, where the numerator represents the arc length and the denominator
#' represents radial distance. This proportion has been referred to as "rectangularity" in the historical literature
#' and it remains constant throughout the design. Given that the inner and outer borders of the growing space shape
#' surrounding a plant in a Nelder design are arcs, and that the shape is not truly rectangular or trapezoidal in
#' nature, the term "rectangularity" can be confusing.
#' @param even Logical indicated whether or not the design should be adjusted so that the angle between spokes goes
#' evenly into \code{max.angle} (i.e. so that there are no spokes that must be removed from the experiment as border
#' spokes).
#' @param max.angle The maximum rotation (in degrees) of the design. If 360, then a full circle design will be created.
#' @param arc.borders Number of border arcs on either extreme.
#' @param spoke.borders Number of border spokes on either extreme (only used if max.angle < 360).
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @references
#' \itemize{
#'  \item Nelder JA (1962) New kinds of systematic designs for spacing experiments.
#'  Biometrics 18:283-307.
#'  \url{http://www.jstor.org/stable/2527473}
#'  \item Parrott DL, Brinks JS, Lhotka JM (2011) Designing Nelder wheel plots for tree density experiments.
#'  New Forests 43:245–254.
#'  \url{http://link.springer.com/10.1007/s11056-011-9278-4}
#' }
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- nelder(DN         = 1000,
#'                D1         = 3000,
#'                N          = 5,
#'                tau        = 1,
#'                even       = TRUE,
#'                max.angle  = 360)
nelder <- function(DN, D1, N, tau = 1, even = FALSE, max.angle = 360, arc.borders = 1, spoke.borders = 1) {

  if(!(is.numeric(D1) & (D1 %% 1 == 0) & length(D1) == 1))         stop("D1 must an integer", call. = FALSE)
  if(!(is.numeric(DN) & (DN %% 1 == 0) & length(DN) == 1))         stop("DN must an integer", call. = FALSE)
  if(!(is.numeric(N)  & (N  %% 1 == 0) & length(N)  == 1 & N > 2)) stop("N must be an integer greater than 2", call. = FALSE)
  if(D1 <= DN)                              stop("D1 must be greater than DN",                    call. = FALSE)
  if(!(is.numeric(tau) & length(tau) == 1)) stop("tau must be numeric and of length 1",           call. = FALSE)
  if(!is.logical(even))                     stop("even must be a logical",                        call. = FALSE)
  if(max.angle > 360 | max.angle <= 0)      stop("max.angle must be between 0 and 360",           call. = FALSE)
  if(!(is.numeric(arc.borders) & (arc.borders %% 1 == 0) & length(arc.borders) == 1 & arc.borders >= 1))
    stop("arc.borders must be an integer greater than 0", call. = FALSE)
  if(!(is.numeric(spoke.borders) & (spoke.borders %% 1 == 0) & length(spoke.borders) == 1 & spoke.borders >= 0))
    stop("spoke.borders must be a positive integer", call. = FALSE)
  if(max.angle < 360 & spoke.borders == 0) stop("spoke.borders must be greater than 0 when max.angle is less than 360", call. = FALSE)

  max.angle.rad <- max.angle * pi / 180  # radians

  ## ORIGINAL VALUES WITHOUT EVEN SPOKE OPTIMIZATION
  tau.orig      <- tau
  alpha.orig    <- exp((log(D1) - log(DN)) / (2 * N - 2))            # Parrot et al. (2012) - Eq. 1
  theta.orig    <- tau.orig * (alpha.orig ^ 0.5 - alpha.orig ^ -0.5) # radians - Parrot et al. 2012 - (Eq. 3)
  n.spokes.orig <- ceiling(max.angle.rad / theta.orig)

  orig.dat <- nelder_calc(alpha         = alpha.orig,
                          theta         = theta.orig,
                          tau           = tau.orig,
                          D1            = D1,
                          N             = N,
                          max.angle     = max.angle,
                          n.spokes      = n.spokes.orig,
                          arc.borders   = arc.borders,
                          spoke.borders = spoke.borders)

  ## NEW VALUES WITH EVEN SPOKE OPTIMIZATION
  n.spokes <- round(max.angle.rad / theta.orig)
  theta    <- max.angle.rad / n.spokes

  get_alpha <- function(alpha) abs((alpha ^ 0.5 - alpha ^ -0.5) - theta / tau.orig)
  alpha <- optimize(get_alpha, c(0, alpha.orig * 10000))$minimum
  tau   <- theta / (alpha ^ 0.5 - alpha ^ -0.5)

  even.dat <- nelder_calc(alpha         = alpha,
                          theta         = theta,
                          tau           = tau,
                          D1            = D1,
                          N             = N,
                          max.angle     = max.angle,
                          n.spokes      = n.spokes,
                          arc.borders   = arc.borders,
                          spoke.borders = spoke.borders)

  if(even) {
    even.spoke.dat <- dplyr::bind_rows(orig.dat$plot, even.dat$plot) %>%
      dplyr::mutate(even = c(FALSE, TRUE)) %>%
      dplyr::select(even, dplyr::everything())

    out <- c(even.dat, list(even.optim = even.spoke.dat))

  } else {
    out <- c(orig.dat, list(even.optim = FALSE))
  }

  class(out) <- c(class(out), "nelder", "sysd")
  return(out)
}

#' Create a Nelder Fan experimental design decision table
#' @description Creates a Nelder Fan experimental design decision table.
#' @details This function helps explore Nelder Fan design options and select a design that meets external constraints
#' (e.g. plant or space availability).
#' Function inputs are identical to \code{\link{nelder}}, but inputs of any length are allowed.
#' All possible combinations of inputs are created using \code{expand.grid}, and then each of these cases is passed to
#' \code{\link{nelder}} for evaluation. Inputs and outputs are all combined and returned for evaluation.
#' @return A tibble containing a wide range of traits of the experimental designs.
#' @param DN Plant density within the last experimental arc (plants ha-1)
#' (i.e. lower extreme of experimental plant density range).
#' @param D1 Plant density within the first experimental arc (plants ha-1)
#' (i.e. upper extreme of experimental plant density range).
#' @param N Number of experimental arcs (i.e. number of densities to be tested within D1 to DN).
#' @param tau The "rectangularity" proportion. See \code{\link{nelder}} for details.
#' @param even Logical indicated whether or not the design should be adjusted so that the angle between spokes goes
#' evenly into \code{max.angle} (i.e. so that there are no spokes that must be removed from the experiment as border
#' spokes).
#' @param max.angle The maximum rotation (in degrees) of the design. If 360, then a full circle design will be created.
#' @param arc.borders Number of border arcs on either extreme.
#' @param spoke.borders Number of border spokes on either extreme (only used if max.angle < 360).
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @family definition functions
#' @examples
#' dat <- nelder_decision(DN         = seq(1000, 2000, 50),
#'                        D1         = 3000,
#'                        N          = 5,
#'                        tau        = 1,
#'                        even       = TRUE,
#'                        max.angle  = c(90, 180, 360))
nelder_decision <- function(DN, D1, N, tau = 1, even = FALSE, max.angle = 360, arc.borders = 1, spoke.borders = 1) {
  test.cases <- expand.grid(DN = DN, D1 = D1, N = N, tau = tau, even = even, max.angle = max.angle,
                            arc.borders = arc.borders, spoke.borders = spoke.borders)

  OUT <- dplyr::tibble()
  for(i in 1:nrow(test.cases)) {
    out <- nelder(DN            = test.cases$DN[i],
                  D1            = test.cases$D1[i],
                  N             = test.cases$N[i],
                  tau           = test.cases$tau[i],
                  even          = test.cases$even[i],
                  max.angle     = test.cases$max.angle[i],
                  arc.borders   = test.cases$arc.borders[i],
                  spoke.borders = test.cases$spoke.borders[i])
    OUT <- dplyr::bind_rows(OUT, out$plot)
  }

  OUT <- dplyr::select(OUT, -max.angle)

  return(dplyr::bind_cols(test.cases, OUT))
}

#' Create a Biculture Nelder Fan experimental design
#' @description Creates a Biculture Nelder Fan experimental design.
#' @details The Nelder Fan or Nelder Wheel Type Ia experimental design (Nelder 1962) is an experimental design that
#' systematically varies plant desnity within a single plot, where the growing area around each plant has a constant
#' shape throughout the design but increases as radius increases. Goelz (2001) adapted this design to simultaneously
#' study the effect of species composition by superimposing a species gradient along the arc (Figure 6 of Goelz, 2001).
#' This function takes a Nelder Fan design from \code{\link{nelder}} and adds species identities to create
#' the Goelz (2001) biculture version. For Nelder Fan designs with \code{max.angle = 360}, species monoculture are
#' set to opposite poles of the circle, and the composition gradient occurs in two direction along either side of the
#' circle in between. For Nelder Fan designs with \code{max angle < 360}, the composition gradient occurs in one
#' direction between monoculture extremes at either edge of the design.
#' This function does NOT robustly impose any conformity criterion to ensure that competition levels are tested evenly across
#' all density levels. Species are assigned randomly to positions within each spoke of the design based on the probability of each
#' species in that spoke. For a robust implementation of the competition-density conformity criterion that optimizes the rough
#' initial approach of this function using an evolutionary algorithm, use \code{\link{nelder_biculture_optim}}.
#' @return An object of classes sysd, nelder, and nelder-biculture. This is a list of 5 elements,
#' the first 3 of which are the same as for \code{\link{nelder}}, and the last 2 of which are:
#' \itemize{
#'  \item{"species.counts"}{ - An abject of class "table" containing the total counts of each species in the design.}
#'  \item{"spoke.composition"}{ - A data frame (tibble) containing the ratio and counts of each species by spoke.}
#' }
#' @param data An object of class nelder.
#' @param comps An optional numeric vector containing the ratios of one species in each spoke.
#' This can effectively be used to create non-standard bi-culture designs that deviate from the Goelz (2001) approach.
#' If \code{NULL}, the default, then a linear
#' gradient of speies composition is used between monoculture extremes for each species.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @references
#' \itemize{
#'  \item Nelder JA (1962) New kinds of systematic designs for spacing experiments.
#'  Biometrics 18:283-307.
#'  \url{http://www.jstor.org/stable/2527473}
#'  \item Goelz J (2001) Systematic experimental designs for mixed species plantings.
#'  Native Plants Journal 2:90–96.
#'  \url{http://npj.uwpress.org/content/2/2/90.short}
#' }
#' @export
#' @family definition functions
#' @examples
#' nelder.design <- nelder(DN = 1000, D1 = 3000, N = 5)
#' dat <- nelder_biculture(data = nelder.design)
nelder_biculture <- function(data, comps = NULL) {

  nelder_class_check(data)

  if(is.null(comps)) { # generate the sequence of species A composition in each spoke
    if(data$plot$max.angle == 360) {
      if(data$plot$spokes %% 2 != 0) stop("when max.angle = 360,
                                            data must have an even number of spokes", call. = FALSE)
      inc.length <- data$plot$spokes / 2 + 1
      inc.seq <- seq(0, 1, length.out = inc.length)
      dec.seq <- rev(inc.seq[c(-1, -length(inc.seq))])
      comps <- c(inc.seq, dec.seq)
    } else {
      if(!(is.numeric(comps) & length(comps) == data$plot$spokes))  stop("comps must be NULL or a numeric vector
                                                                          of length equal to the number of spokes
                                                                          in the data", call. = FALSE)
      comps <- seq(0, 1, length.out = data$plot$spokes)
    }
  }

  ## ASSIGN SPECIES
  data <- assign_species(data = data, comps = comps)

  ## STATS
  data$plants$species    <- factor(data$plants$species)
  data$species.counts    <- table(data$plants$species)
  data$spoke.composition <- dplyr::tibble(spoke   = 1:data$plot$spokes,
                                          A.ratio = comps,
                                          B.ratio = 1 - A.ratio,
                                          A.count = round(comps * data$plot$arcs),
                                          B.count = data$plot$arcs - A.count)

  class(data) <- c(class(data), "nelder-biculture")
  return(data)
}

#' Create a compeition-conformity-optimized Biculture Nelder Fan experimental design
#' @description Creates a competition-conformity-optimized Biculture Nelder Fan experimental design.
#' @details While \code{\link{nelder_biculture}} creates Biculture Nelder Fan designs that follow a composition-conformity criterion
#' for species composition across spokes, this function creates designs that are optimized to also conform to a compeition-conformity
#' criterion  for species competition environments across arcs using an evolutionary algorithm. Function parameters other than
#' \code{data} are all controls on the evolutionary algorithm.
#' The \code{\link{ecr}} package is used for the evolutionary algorithm.
#' @return An list containing:
#' \itemize{
#'  \item{"layout"}{ - A data frame (tibble) containing the base nelder design.}
#'  \item{"stats"}{ - A data frame (tibble) containing statistics on each generation in the evolutionary algorithm.}
#'  \item{"data"}{ - A data frame (tibble) containing the actual data (designs) of each Biculture Nelder Fan in eac generation.}
#' }
#' @param data An object of class nelder-biculture.
#' @param save.path A character string indicating the path for where to save data after each generation.
#' @param MU The population size.
#' @param LAMBDA The number of offspring to produce in each generation.
#' @param MAX.GEN The number of generations to run the evolutionary algorithm.
#' @param P.RECOMB The probability of two parents to perform crossover.
#' @param RECOMB The crossover probability for each gene.
#' @param P.MUT The probability to apply mutation to a child.
#' @param MUT The mutation probability for each gene.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- nelder()
#' dat.bi <- nelder_biculture(data = dat)
#' dat.bi.optim <- nelder_biculture_optim(data = dat.bi)
nelder_biculture_optim <- function(data,
                                   save.path = ".",
                                   MU        = 100,
                                   LAMBDA    = MU,
                                   MAX.GEN   = 150,
                                   P.RECOMB  = 1,
                                   RECOMB    = 0.1,
                                   P.MUT     = 1,
                                   MUT       = 0.05) {

  nelder_biculture_class_check(data)
  if(!requireNamespace("ecr", quietly = TRUE)) stop("The package 'ecr' is required for nelder_biculture() if
                                                      optim = TRUE. Please install and load it", call. = FALSE)

  ### GA CONTROLS
  control <- ecr::initECRControl(fitness.fun  = nelder_fitness,
                                 n.objectives = 1,
                                 minimize     = TRUE) %>%
    ecr::registerECROperator(slot = "mutate",
                             fun  = ecr::setup(ecr::makeMutator(mutShuffle, supported = "float"),
                                               p = MUT))

  ## INITIAL POPULATION
  data$plants$species <- NULL
  FAN <- data
  FAN$even.optim <- FAN$species.counts <- FAN$spoke.composition <- NULL

  make_spoke_chromosomes <- function(x) {
    x$plants %>%
      dplyr::select(spoke, species) %>%
      dplyr::mutate(species = match(species, LETTERS)) %>%
      dplyr::group_by(spoke) %>%
      tidyr::nest() %>%
      .$data %>%
      purrr::map("species")
  }

  population <- purrr::map(rep(list(FAN), MU),
                           assign_species,
                           comps = data$spoke.composition$A.ratio) %>%
    purrr::map(make_spoke_chromosomes)

  ## GA
  ga.out <- run_ga(control     = control,
                   population  = population,
                   layout      = FAN,
                   layout.comp = FAN$plants,
                   LAMBDA      = LAMBDA,
                   MAX.GEN     = MAX.GEN,
                   P.RECOMB    = P.RECOMB,
                   RECOMB      = RECOMB,
                   P.MUT       = P.MUT,
                   save.path   = save.path)

  OUT <- list(layout = data, stats = ga.out$pop.stats, data = ga.out$pop.data)
  class(OUT) <- c(class(OUT), "nelder-optim")
  return(OUT)
}

#' Select optimal Biculture Nelder Fan designs from the results of nelder_biculture_optim
#' @description Selects optimal Biculture Nelder Fan designs from the results of \code{\link{nelder_biculture_optim}}.
#' @return A list of objects of class nelder-biculture.
#' @param optim.results An object of class nelder-optim.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' optim.dat <- nelder_biculture_optim()
#' my_nelder <- select_optimal_nelder_biculture(optim.results = optim.dat)
select_optimal_nelder_biculture <- function(optim.results) {

  if(!("nelder-optim" %in% class(optim.results))) stop("optim.results must be of class nelder-optim", call. = FALSE)

  best.designs <- select_optimal(optim.results)

  OUT <- list()
  for(i in 1:length(best.designs)) {
    optim.results$layout$plants <- optim.results$layout$plants %>%
      dplyr::arrange(spoke, arc) %>%
      dplyr::mutate(species = factor(LETTERS[best.designs[[i]]$species]))
    OUT <- c(OUT, list(optim.results$layout))
  }

  return(OUT)
}

#' Calculate competition indices for a Nelder Fan experimental design
#' @description Calculates competition indices for a Nelder Fan experimental design.
#' @return A data frame (tibble), which is the \code{plants} element of the nelder object passed to \code{data},
#' but with 4 additioanl columns:
#' \itemize{
#'  \item{"A.inv.dist"}{ - The inverse distance weighted competition felt by that individual due to nearby individuals
#'  of species A}
#'  \item{"B.inv.dist"}{ - The inverse distance weighted competition felt by that individual due to nearby individuals
#'  of species B}
#'  \item{"A.neighbors"}{ - The number of directly adjacent individuals of species A}
#'  \item{"B.neighbors"}{ - The number of directly adjacent individuals of species B}
#' }
#' @param data An object of class nelder-biculture.
#' @param search.radius Search radius to use for calculating inverse distance weighted competition (m). If \code{NULL},
#' the default, then inverse distance weighted competition is not calculated.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- nelder(DN         = 1000,
#'               D1         = 3000,
#'               N          = 5,
#'               tau        = 1,
#'               even       = TRUE,
#'               max.angle  = 360)
#' dat.bi <- nelder_biculture(dat)
#' nelder_biculture_competition(dat.bi)
nelder_biculture_competition <- function(data, search.radius = NULL) {

  nelder_biculture_class_check(data)
  if(!(is.null(search.radius) | (is.numeric(search.radius) & length(search.radius) == 1)))
    stop("search.readius must be NULL or a numeric vector of length 1")

  data$plants$A.inv.dist  <- NA_real_
  data$plants$B.inv.dist  <- NA_real_
  data$plants$A.neighbors <- NA_real_
  data$plants$B.neighbors <- NA_real_

  data$plants$search.radius <- search.radius

  for(i in 1:nrow(data$plants)) {
    focal.indiv <- data$plants[i, ]

    if(!is.null(search.radius)) {
      ## Inverse Distance
      X <- focal.indiv$x.field
      Y <- focal.indiv$y.field

      inverse.distance <- data$plants[-i, ] %>%
        dplyr::mutate(d = sqrt((x.field - X)^2 + (y.field - Y)^2)) %>%
        dplyr::filter(d <= search.radius) %>%
        dplyr::mutate(id = 1 / (1 + d)) %>%
        dplyr::group_by(species) %>%
        dplyr::summarize(id = sum(id)) %>%
        tidyr::spread(key = "species", value = "id")

      if(!("A" %in% names(inverse.distance))) inverse.distance$A <- 0
      if(!("B" %in% names(inverse.distance))) inverse.distance$B <- 0

      data$plants$A.inv.dist[i] <- inverse.distance$A
      data$plants$B.inv.dist[i] <- inverse.distance$B
    }

    ## Neighbor Counts
    s <- focal.indiv$spoke
    a <- focal.indiv$arc

    neighbor.spokes <- c(s - 1, s, s + 1)
    neighbor.arcs   <- c(a - 1, a, a + 1)

    if(data$plot$max.angle == 360) {
      neighbor.spokes[neighbor.spokes < 1] <- neighbor.spokes[neighbor.spokes < 1] + data$plot$spokes
      neighbor.spokes[neighbor.spokes > data$plot$spokes] <- neighbor.spokes[neighbor.spokes > data$plot$spokes] - data$plot$spokes
    } else {
      neighbor.spokes[neighbor.spokes < 1] <- NA
      neighbor.spokes[neighbor.spokes > data$plot$spokes] <- NA
    }

    neighbor.arcs[neighbor.arcs < 0] <- NA
    neighbor.arcs[neighbor.arcs > (data$plot$arcs - 1)] <- NA

    neighbors <- data$plants[-i, ] %>%
      dplyr::filter(arc   %in% neighbor.arcs) %>%
      dplyr::filter(spoke %in% neighbor.spokes) %>%
      dplyr::summarize(A = sum(species == "A"),
                       B = sum(species == "B"))

    if(!(a %in% c(0, data$plot$arcs - 1))) {
      data$plants$A.neighbors[i] <- neighbors$A
      data$plants$B.neighbors[i] <- neighbors$B
    }
  }

  return(data$plants)
}

#' Calculate distance between two individuals in the same arc and adjacent spokes of a Nelder Fan experimental design
#' @description Calculates distance between two individuals in the same arc and adjacent spokes of a Nelder Fan experimental design.
#' @return A numeric value of the distance (m).
#' @param data An object of class nelder.
#' @param arc An numeric value indicating which arc of the design to calculate distance for.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- nelder(DN         = 1000,
#'               D1         = 3000,
#'               N          = 5,
#'               tau        = 1,
#'               even       = TRUE,
#'               max.angle  = 360)
#' nelder_interspoke_distance(dat)
nelder_interspoke_distance <- function(data, a = data$plot$arcs) {

  nelder_class_check(data)

  two.points <- data$plants %>%
    dplyr::filter(arc == a) %>%
    dplyr::filter(spoke %in% 1:2)

  out <- sqrt((two.points$x.field[1] - two.points$x.field[2])^2 + (two.points$y.field[1] - two.points$y.field[2])^2)

  return(out)
}

#' Calculate Nelder Fan design
#' @description Calculates Nelder Fan design.
#' Used within \code{\link{nelder}}.
#' @return A list containing the plant and plot data.
#' @param alpha Rate of change along the spokes.
#' @param theta The angle between the spokes.
#' @param tau The "rectangularity" proportion
#' @param D1 Plant density within the first experimental arc (plants ha-1)
#' (i.e. upper extreme of experimental plant density range).
#' @param N Number of experimental arcs (i.e. number of densities to be tested within D1 to DN).
#' @param max.angle The maximum rotation (in degrees) of the design. If 360, then a full circle design will be created.
#' @param n.spokes The number of spokes.
#' @param arc.borders Number of border arcs on either extreme.
#' @param spoke.borders Number of border spokes on either extreme (only used if max.angle < 360).
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @importFrom dplyr %>%
#' @keywords internal
nelder_calc <- function(alpha, theta, tau, D1, N, max.angle, n.spokes, arc.borders, spoke.borders) {

  n.arcs <- N + 2 * arc.borders
  spoke.borders <- ifelse(max.angle < 360, spoke.borders, 0)

  if((2 * spoke.borders) > n.spokes) stop("there are too many spoke borders for this size design", call. = FALSE)

  spoke.border.ids <- c(head(1:n.spokes, spoke.borders), tail(1:n.spokes, spoke.borders))
  arc.border.ids   <- c(head(1:n.arcs,    arc.borders),  tail(1:n.arcs,   arc.borders))

  Co <- (2 - (alpha + alpha ^ -1)) / (2 * (alpha - alpha ^ -1)) * 100 # % - Parrot et al. (2012) - Eq. 2
  r0 <- sqrt(20000 / (D1 * theta * (alpha ^ 3 - alpha)))              # m - Parrot et al. (2012) - Eq. 4

  arc.dat <- dplyr::tibble(arc     = 1:n.arcs,
                           r       = r0 * alpha ^ (arc - 1 - (arc.borders - 1)), # m           - Parrot et al. (2012) - Eq. 5 (modified **)
                           space   = theta * r ^ 2 * (alpha - alpha ^ -1) / 2,   # m2          - Parrot et al. (2012) - Eq. 6
                           density = 10000 / space)                              # plants ha-1 - Parrot et al. (2012) - Eq. 7
  # ** The original Eq. 5 from Parrot et al. (2012) reads: r = r0 * alpha ^ arc
  # This equation needs to be modified to adjust for two things:
  #     (1) arc ids are specified here starting with 1, rather than 0.
  #     (2) D1 is the density of the innermost *experimental* arc. Parrot et al. (2012) assumes that the design has *one* border
  #     arc interior to this innermost experimental arc and calculates r0 assuming it is the radious to this single border arc.
  #     Here, flexibility in the number of border arcs must be allowed via the input parameter arc.borders.
  # To account for #1, 1 is substracted in the expondnet.
  # To account for #2, (arc.borders - 1) is subsracted in the exponent (-1 is included because Eq. 5 already assumes 1 border arc).

  arc.dat[arc.border.ids, c("space", "density")] <- NA

  plot.dat <- dplyr::tibble(area            = pi * max(arc.dat$r) ^ 2 / 10000 * max.angle / 360, # ha - Parrot et al. (2012) - Eq. 8 (modified **)
                            arcs            = n.arcs,
                            exp.arcs        = N,
                            spokes          = n.spokes,
                            exp.spokes      = n.spokes - 2 * spoke.borders,
                            rmin            = min(arc.dat$r),
                            rmax            = max(arc.dat$r),
                            alpha           = alpha,
                            angle           = theta * 180 / pi, # convert to degrees
                            max.angle       = max.angle,
                            rectangularity  = tau,
                            non.centrality  = Co,
                            min.exp.density = min(arc.dat$density, na.rm = TRUE),
                            max.exp.density = max(arc.dat$density, na.rm = TRUE),
                            plants          = arcs * spokes,
                            exp.plants      = exp.arcs * exp.spokes,
                            perc.exp        = round(exp.plants / plants * 100))
  # ** The original Eq. 8 from Parrot et al. (2012) reads: area = pi * r(N+1) ^ 2 / 10000
  # This equation needs to be modified to adjust for two things:
  #     (1) Parrot et al. (2012) assumes that the design has *one* border arc exterior to the outermost experimental arc (the N+1 in the Eq).
  #     Here, flexibility in the number of border arcs must be allowed via the input parameter arc.borders.
  #     (2) Parrot et al. (2012) assumes a full 360-degree design.
  #     Here, flexiblity in the shape of the design must be allowed via the input parameter max.angle
  # To account for #1, max(arc.dat$r) is used to replace r(N+1).
  # To account for #2, (max.angle / 360) is multiplied as well.

  plant.dat <- dplyr::tibble(spoke  = rep(1:plot.dat$spokes, each  = plot.dat$arcs),
                             arc    = rep(arc.dat$arc,       times = plot.dat$spokes),
                             theta  = rep(plot.dat$angle * 0:(plot.dat$spokes - 1), each = plot.dat$arcs)) %>%
    dplyr::left_join(arc.dat, by = "arc") %>%
    dplyr::mutate(x.field = r * cos(theta * pi / 180)) %>%
    dplyr::mutate(y.field = r * sin(theta * pi / 180)) %>%
    dplyr::mutate(exp     = (spoke %in% spoke.border.ids) | (arc %in% arc.border.ids)) %>%
    dplyr::mutate(exp     = factor(exp, levels = c(FALSE, TRUE), labels = c("Experimental", "Border")))

  temp.dat <- list(plants = plant.dat)
  class(temp.dat) <- "nelder"
  plot.dat <- dplyr::mutate(plot.dat, outer.interspoke.distance = nelder_interspoke_distance(data = temp.dat, a = arcs))

  return(list(plants = plant.dat, plot = plot.dat))
}

#' Calculate the fitness of a Nelder Fan biculture design
#' @description Calculates the fitness of a Nelder Fan biculture design. Fitness is calculated to reflect the
#' uniformity of compeition environments present across experimental arcs.
#' Used within \code{\link{nelder_biculture}}.
#' @return A numeric value.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @importFrom dplyr %>%
#' @keywords internal
nelder_fitness <- function(design, layout) {

  layout$plants <- layout$plants %>%
    dplyr::mutate(species = LETTERS[unlist(design)])

  complete.combos <- expand.grid(species     = c("A", "B"),
                                 arc         = 1:layout$plot$exp.arcs,
                                 A.neighbors = 0:8) %>%
    dplyr::as_tibble()

  out <- layout %>%
    nelder_biculture_competition() %>%
    dplyr::group_by(species, arc, A.neighbors) %>%
    dplyr::summarize(dens = dplyr::n()) %>%
    dplyr::filter(!is.na(A.neighbors)) %>%
    dplyr::right_join(complete.combos, by = c("species", "arc", "A.neighbors")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dens = tidyr::replace_na(dens, 0)) %>%
    dplyr::group_by(species, A.neighbors) %>%
    #dplyr::summarize(stdev = diff(range(dens))) %>%
    dplyr::summarize(stdev = sd(dens)) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(fitness = mean(stdev)) %>%
    as.numeric()

  return(out)
}

#' Assign speices identities in a Nelder Fan biculture design
#' @description Assigns speices identities in a Nelder Fan biculture design.
#' Used within \code{\link{nelder_biculture}}.
#' @return An object of the same class as \code{data}.
#' @param data An object of class nelder or class nelder-biculture.
#' @param comps An optional numeric vector containing the ratios of one species in each spoke.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @keywords internal
assign_species <- function(data, comps) {
  data$plants$species <- "B"              # initialize all spots as species B
  quants <- round(comps * data$plot$arcs) # determine # of spots in each spoke that whould be species A

  for(s in 1:data$plot$spokes) {          # randomly select species A spots
    q <- quants[s]
    a <- sample(x       = 1:data$plot$arcs,
                size    = q,
                replace = FALSE)
    data$plants$species[which(data$plants$spoke == s & data$plants$arc %in% a)] <- "A"
  }

  return(data)
}

#' Check if an object if of class nelder
#' @description Checks if an object if of class nelder. Used by many nelder definition functions.
#' @return An error.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @keywords internal
nelder_class_check <- function(data) {
  if(!("nelder" %in% class(data))) stop("data must be of class 'nelder'", call. = FALSE)
}

#' Check if an object if of class nelder-biculture
#' @description Checks if an object if of class nelder-biculture. Used by many nelder functions.
#' @return An error.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @keywords internal
nelder_biculture_class_check <- function(data) {
  if(!("nelder-biculture" %in% class(data))) stop("data must be of class 'nelder-biculture'", call. = FALSE)
}
