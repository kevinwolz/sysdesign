#' Convert yards to feet with fractional inches
#' @description Converts yards to feet with fractional inches.
#' @return A character string
#' @param y A numeric value of length 1.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
yard_to_feet_frac <- function(y) {

  if(length(y) > 1) stop("y must be of length 1", call. = FALSE)

  dec_to_frac <- dplyr::tibble(dec  = (0:15)/16,
                               frac = c("", "1/16", "1/8", "3/16",  "1/4", "5/16",  "3/8", "7/16", "1/2",
                                        "9/16", "5/8", "11/16", "3/4", "13/16", "7/8", "15/16"))
  orig <- y * 3
  feet <- floor(orig)
  inch.orig <- (orig - feet) * 12
  inches    <- floor(inch.orig)
  dec <- inch.orig - inches
  frac <- ""
  if(dec != 0) frac <- paste0(" ", dec_to_frac$frac[which.min(abs(dec - dec_to_frac$dec))])

  inch <- paste0(" ", inches, frac, "in")
  if(inches == 0 & frac == "") inch <- ""

  return(paste0(feet, "ft", inch))
}

#' Run genetic algorithm
#' @description Runs genetic algorithm for \code{\link{goelz_optim}} & \code{\link{nelder_biculture_optim}}.
#' @return A list containing fitness stats and population data.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @keywords internal
run_ga <- function(control, population, layout, layout.comp, LAMBDA, MAX.GEN, P.RECOMB, RECOMB, P.MUT, save.path, start.gen) {

  parallelMap::parallelRegisterLevels(package = "ecr",
                                      levels  = c("evaluateFitness", "generateOffspring",
                                                  "computeDominanceRanking", "computeIndicators"))

  control <- control %>%
    ecr::registerECROperator(slot = "recombine",
                             fun  = ecr::setup(ecr::recUnifCrossover, p = RECOMB)) %>%
    ecr::registerECROperator(slot = "selectForMating",
                             fun  = ecr::setup(ecr::selTournament,
                                               k = 2L)) %>%
    ecr::registerECROperator(slot = "selectForSurvival",
                             fun  = ecr::selGreedy)

  GEN <- start.gen
  print(paste("STARTING GENERATION:", GEN))
  startGen <- proc.time()[3]

  fitness <- ecr::evaluateFitness(control = control,
                                  inds    = population,
                                  layout  = layout)
  population <- annotate_pop(population = population,
                             fitness    = fitness,
                             gens       = rep(GEN, length(population)),
                             gen.ids    = paste(GEN, 1:length(population), sep = "-"))

  GAout <- compile_pop(pop = population, fit = fitness, layout = layout.comp, GEN = GEN)
  pop.stats <- GAout$stats
  pop.data  <- GAout$data

  readr::write_csv(pop.stats, paste0(save.path, "/sysdesign_ga_stats.csv"), append = (start.gen != 1))
  readr::write_csv(pop.data,  paste0(save.path, "/sysdesign_ga_data.csv"),  append = (start.gen != 1))

  elapsedGen = round((proc.time()[3] - startGen) / 60, 2)
  print(paste0("DONE WITH  GENERATION: ", GEN, " (", elapsedGen, " minutes)"))

  for(GEN in (start.gen + 1):MAX.GEN) {
    print(paste("STARTING GENERATION:", GEN))
    startGen <- proc.time()[3]

    ## GENERATE NEW POPULATION
    offspring <- ecr::recombinate(control  = control,
                                  inds     = population,
                                  fitness  = fitness,
                                  lambda   = LAMBDA,
                                  p.recomb = P.RECOMB) %>%
      ecr::mutate(control = control,
                  inds    = .,
                  p.mut   = P.MUT)

    ## EVALUATE FITNESS OF OFFSPRING
    fitness <- ecr::evaluateFitness(control = control,
                                    inds    = offspring,
                                    layout  = layout)
    offspring <- annotate_pop(population = offspring,
                              fitness    = fitness,
                              gens       = rep(GEN, length(offspring)),
                              gen.ids    = paste(GEN, 1:length(offspring), sep = "-"))

    ## SELECT FOR SURVIVAL
    sel <- ecr::replaceMuPlusLambda(control    = control,
                                    population = population,
                                    offspring  = offspring)
    population <- sel$population
    fitness    <- sel$fitness

    ## SAVE DATA
    GAout <- compile_pop(pop = population, fit = fitness, layout = layout.comp, GEN = GEN)
    pop.stats <- dplyr::bind_rows(pop.stats, GAout$stats)
    pop.data  <- dplyr::bind_rows(pop.data,  GAout$data)

    readr::write_csv(GAout$stats, paste0(save.path, "/sysdesign_ga_stats.csv"), append = TRUE)
    readr::write_csv(GAout$data,  paste0(save.path, "/sysdesign_ga_data.csv"),  append = TRUE)

    elapsedGen = round((proc.time()[3] - startGen) / 60, 2)
    print(paste0("DONE WITH  GENERATION: ", GEN, " (", elapsedGen, " minutes)"))
  }

  return(list(pop.stats = pop.stats, pop.data = pop.data))
}

#' Integer mutation
#' @description Performs integer mutation between a lower and upper bound. Used by \code{\link{goelz_optim}}.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @keywords internal
mutInteger <- function(ind, p, lower, upper) {
  if(!(is.integer(lower) & length(lower) == 1)) stop("lower must be an integer of length 1",     call. = FALSE)
  if(!(is.integer(upper) & length(upper) == 1)) stop("upper must be an integer of length 1",     call. = FALSE)
  if(length(lower) != length(upper)) stop("Length of lower and upper bounds need to be equal!",  call. = FALSE)

  for(idx in 1:length(ind)) {
    available <- seq(from = lower, to = upper, by = 1)
    available <- available[available != ind[idx]]
    if(runif(1L) < p) ind[idx] <- sample(available, size = 1)
  }
  return(ind)
}

#' Shuffle mutation
#' @description Performs mutation by shuffling elements within a vector. Used by \code{\link{nelder_biculture_optim}}.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @keywords internal
mutShuffle <- function(ind, p) {
  for(idx in 1:length(ind)) if(runif(1L) < p) ind[[idx]] <- sample(ind[[idx]])
  return(ind)
}

#' Annotate population with fitness and stats
#' @description Annotates population with fitness and stats.
#' Used by \code{\link{goelz_optim}} & \code{\link{nelder_biculture_optim}}.
#' @return A popluation matrix.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @keywords internal
annotate_pop <- function(population, fitness, gens, gen.ids) {

  for(i in seq_along(population)) {
    attr(population[[i]], "fitness")    <- fitness[, i]
    attr(population[[i]], "generation") <- gens[i]
    attr(population[[i]], "gen.id")     <- gen.ids[i]
  }

  return(population)
}

#' Compile popluation data
#' @description Compiles popluation data.
#' Used by \code{\link{goelz_optim}} & \code{\link{nelder_biculture_optim}}.
#' @return A list.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @importFrom dplyr %>%
#' @keywords internal
compile_pop <- function(pop, fit, layout, GEN) {

  orig.gen <- pop %>%
    purrr::map_dbl(attr, which = "generation") %>%
    dplyr::as_tibble() %>%
    dplyr::rename(orig.gen = value)

  gen.id <- pop %>%
    purrr::map_chr(attr, which = "gen.id") %>%
    dplyr::as_tibble() %>%
    dplyr::rename(gen.id = value)

  fit <- fit %>%
    t() %>%
    as.data.frame() %>%
    dplyr::as_tibble() %>%
    dplyr::rename(fitness = V1)

  stats <- bind_cols(dplyr::tibble(gen = rep(GEN, nrow(orig.gen))), orig.gen, gen.id, fit)

  bind_layout <- function(x, lay) dplyr::bind_cols(lay, x)

  data <- pop %>%
    purrr::map(matrix, ncol = 1) %>%
    purrr::map(dplyr::as_tibble, .name_repair = "minimal") %>%
    purrr::map(setNames, nm = "species") %>%
    purrr::map(tidyr::unnest, cols = c(species)) %>% # addition for Nelder
    purrr::map(bind_layout, lay = layout) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(orig.gen = rep(orig.gen$orig.gen, each = nrow(layout))) %>%
    dplyr::mutate(gen.id   = rep(gen.id$gen.id,     each = nrow(layout))) %>%
    dplyr::mutate(gen = GEN) %>%
    dplyr::select(gen, orig.gen, gen.id, dplyr::everything()) # goelz: zone, zone.id, id, x, y, z, species

  return(list(stats = stats, data = data))
}

#' Select optimal individual
#' @description Selects the optimal individual from an object of class goelz-optim or nelder-optim.
#' Used by \code{\link{select_optimal_goelz}} & \code{\link{select_optimal_nelder_biculture}}.
#' @return An numeric vector of the optimal species composition.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @importFrom dplyr %>%
#' @keywords internal
select_optimal <- function(x) {
  best.design.ids <- x$stats %>%
    dplyr::filter(fitness == min(fitness)) %>%
    dplyr::select(gen.id, fitness) %>%
    dplyr::distinct() %>%
    .$gen.id

  best.designs <- list()
  for(i in best.design.ids) {
    best.design <- x$data %>%
      dplyr::filter(gen.id == i) %>%
      dplyr::filter(gen == orig.gen) %>%
      dplyr::select(-gen, -orig.gen)

    if("zone" %in% names(best.design)) {
      best.design <- dplyr::arrange(best.design, zone, zone.id)
    } else {
      best.design <- dplyr::arrange(best.design, spoke, arc)
    }
    best.designs <- c(best.designs, list(best.design))
  }

  return(best.designs)
}

#' Add a column with species names
#' @description Adds a column with species names.
#' @return An object of the same class as data.
#' @param data An object of class nelder-biculture or goelz.
#' @param sp.names A vector of species names. Must be of length 2 when data is of class nelder-biculture and length 3 when data is of class goelz.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
name_species <- function(data, sp.names) {
  if("nelder-biculture" %in% class(data)) {
    if(length(sp.names) != 2) stop("sp.names must be of length 2 when data is of class nelder-biculture", call. = FALSE)
    data$plants <- dplyr::mutate(data$plants, species.name = as.character(factor(species, labels = sp.names)))
  } else if("goelz" %in% class(data)) {
    if(length(sp.names) != 3) stop("sp.names must be of length 3 when data is of class nelder-biculture", call. = FALSE)
    data <- dplyr::mutate(data, species.name = as.character(factor(species, labels = sp.names)))
  } else {
    stop("data must be of class nelder-biculture or goelz", call. = FALSE)
  }

  return(data)
}
