
create_goelz <- function(N = 35, reps = 1, split = FALSE) {

  L       <- seq(from = 0, to = 100, length.out = N)
  SPECIES <- 1:3

  triangle <- expand.grid(x = L, y = L, z = L) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(s = x + y + z) %>%
    dplyr::filter(abs(s - 100) < 1) %>%
    dplyr::arrange(y, z, x) %>%
    dplyr::mutate(id = 1:nrow(.)) %>%
    dplyr::select(id, x, y, z) %>%
    dplyr::mutate(y.pos = match(y, L)) %>%
    dplyr::mutate(y.field = sqrt(3) / 2 * (y.pos - 1)) %>%
    dplyr::group_by(y.pos) %>%
    dplyr::mutate(x.pos = match(z, L)) %>%
    dplyr::mutate(x.field = x.pos + 0.5 * (y.pos - 1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(border = FALSE)

  A <- triangle %>%
    dplyr::filter(x > 35 & y < 35) %>%
    dplyr::arrange(y, z) %>%
    dplyr::mutate(trident = "A") %>%
    dplyr::mutate(trident.id = 1:nrow(.))

  B <- triangle %>%
    dplyr::filter(y > 35 & z < 35) %>%
    dplyr::arrange(z, x) %>%
    dplyr::mutate(trident = "B") %>%
    dplyr::mutate(trident.id = 1:nrow(.))

  C <- triangle %>%
    dplyr::filter(z > 35 & x < 35) %>%
    dplyr::arrange(x, y) %>%
    dplyr::mutate(trident = "C") %>%
    dplyr::mutate(trident.id = 1:nrow(.))

  triangle <- dplyr::bind_rows(A, B, C)

  init_species <- function(x, y, z) {
    OUT <- NULL
    for(i in 1:length(x)) {
      out <- sample(SPECIES, size = 1, replace = TRUE, prob = c(x[i], y[i], z[i]))
      OUT <- c(OUT, out)
    }
    return(OUT)
  }

  out <- list()
  A.out <- list()
  index <- 1
  if(split) {
    while(index <= reps) {
      A.design <- triangle %>%
        dplyr::filter(trident == "A") %>%
        dplyr::mutate(species = init_species(x, y, z)) %>%
        .$species

      B.design <- add_one(A.design)
      C.design <- add_one(B.design)

      out <- c(out, list(c(A.design, B.design, C.design)))
      A.out <- c(A.out, list(A.design))
      index <- index + 1
    }
    out <- list(triangle = triangle, design = out, A.design = A.out)
  } else {
    while(index <= reps) {
      A.design <- triangle %>%
        dplyr::filter(trident == "A") %>%
        dplyr::mutate(species = init_species(x, y, z))

      B.design <- triangle %>%
        dplyr::filter(trident == "B") %>%
        dplyr::mutate(species = add_one(A.design$species))

      C.design <- triangle %>%
        dplyr::filter(trident == "C") %>%
        dplyr::mutate(species = add_one(B.design$species))

      out <- c(out, list(dplyr::bind_rows(A.design, B.design, C.design)))
      index <- index + 1
    }
  }

  return(out)
}

add_goelz_border <- function(goelz, n) {
  i <- 0
  while(i < n) {
    bottom.start <- goelz %>%
      dplyr::filter(y.pos == min(y.pos)) %>%
      dplyr::select(x.pos, y.pos, x.field, y.field, trident, species) %>%
      dplyr::mutate(trident = "buffer-bottom") %>%
      dplyr::mutate(x.field = x.field - 1.5) %>%
      dplyr::mutate(y.field = y.field - (sqrt(3) / 2)) %>%
      dplyr::mutate(x.pos   = x.pos - 1) %>%
      dplyr::mutate(y.pos   = y.pos - 1)

    bottom.end <- bottom.start %>%
      dplyr::arrange(x.pos) %>%
      tail(2) %>%
      dplyr::mutate(x.field = x.field + 2) %>%
      dplyr::mutate(x.pos   = x.pos   + 2)

    bottom <- dplyr::bind_rows(bottom.start, bottom.end) %>%
      dplyr::mutate(border = TRUE)


    left.start <- goelz %>%
      dplyr::filter(x.pos == min(x.pos)) %>%
      dplyr::select(x.pos, y.pos, x.field, y.field, trident, species) %>%
      dplyr::mutate(trident = "buffer-left") %>%
      dplyr::mutate(x.field = x.field - 1) %>%
      dplyr::mutate(x.pos = x.pos - 1)

    left.end <- left.start %>%
      dplyr::arrange(y.pos) %>%
      tail(2) %>%
      dplyr::mutate(x.field = x.field + 1.0) %>%
      dplyr::mutate(y.field = y.field + 2 * (sqrt(3) / 2)) %>%
      dplyr::mutate(y.pos = y.pos + 2)

    left <- dplyr::bind_rows(left.start, left.end) %>%
      dplyr::mutate(border = TRUE)


    right.start <- goelz %>%
      dplyr::group_by(y.pos) %>%
      dplyr::filter(x.pos == max(x.pos)) %>%
      dplyr::ungroup() %>%
      dplyr::select(x.pos, y.pos, x.field, y.field, trident, species) %>%
      dplyr::mutate(trident = "buffer-right") %>%
      dplyr::mutate(x.field = x.field + 0.5) %>%
      #dplyr::mutate(x.pos = x.pos + 1) %>%
      dplyr::mutate(y.field = y.field + (sqrt(3) / 2)) %>%
      dplyr::mutate(y.pos = y.pos + 1)

    right.end <- right.start %>%
      dplyr::arrange(dplyr::desc(y.pos)) %>%
      tail(2) %>%
      dplyr::mutate(x.field = x.field + 1) %>%
      dplyr::mutate(x.pos = x.pos + 2) %>%
      dplyr::mutate(y.field = y.field - 2 * (sqrt(3) / 2)) %>%
      dplyr::mutate(y.pos = y.pos - 2)

    right <- dplyr::bind_rows(right.start, right.end) %>%
      dplyr::mutate(border = TRUE)

    goelz <- dplyr::bind_rows(goelz, bottom, left, right)
    i <- i + 1
  }

  return(goelz)
}

plot_goelz <- function(data,
                       species.colors = c("black", "white", "red"),
                       border.colors  = c("black", "green")){
  out <- ggplot(data, aes(x     = x.field,
                          y     = y.field,
                          fill  = factor(species),
                          color = border)) +
    labs(fill = "Species", color = "Border") +
    geom_point(size = 4, shape = 21) +
    scale_fill_manual(values  = species.colors) +
    scale_color_manual(values = border.colors) +
    theme_void() +
    coord_equal()

  # ggtern(g, aes(x = x, y = y, z = z, fill = factor(species))) +
  #   geom_point(size = 2, shape = 21) +
  #   #geom_text(aes(label = x.pos), size = 2) +
  #   scale_fill_manual(values = c("grey50", "white", "red ")) +
  #   theme_void() +
  #   theme(tern.axis.line = element_blank(),
  #         panel.background = element_blank())

  return(out)
}

optim_goelz <- function(N        = 35,
                        MU       = 100, # Population size
                        LAMBDA   = MU,  # Number of offspring to produce in each generation
                        MAX.GEN  = 150, # Number of generations to run GA
                        P.RECOMB = 1,   # Probability of two parents to perform crossover
                        RECOMB   = 0.1, # Crossover probability for each gene
                        P.MUT    = 1,   # Probability to apply mutation to a child
                        MUT      = 0.005) {# Probability of mutation of each gene

  ### GA CONTROLS
  mutInteger <- function (ind, p, lower, upper) {
    assertInteger(lower, any.missing = FALSE, all.missing = FALSE)
    assertInteger(upper, any.missing = FALSE, all.missing = FALSE)
    if (length(lower) != length(upper)) {
      stopf("Uniform mutator: length of lower and upper bounds need to be equal!")
    }

    for(idx in 1:length(ind)) {
      available <- seq(from = lower, to = upper, by = 1)
      available <- available[available != ind[idx]]
      if(runif(1L) < p) ind[idx] <- sample(available, size = 1)
    }
    return(ind)
  }
  mutInteger <- ecr::makeMutator(mutInteger, supported = "float")

  control <- ecr::initECRControl(fitness.fun  = goelz_fitness,
                                 n.objectives = 1,
                                 minimize     = TRUE) %>%
    ecr::registerECROperator(slot = "mutate",
                             fun  = ecr::setup(mutInteger,
                                               p     = MUT,
                                               lower = 1L,
                                               upper = 3L)) %>%
    ecr::registerECROperator(slot = "recombine",
                             fun  = ecr::setup(ecr::recUnifCrossover, p = RECOMB)) %>%
    ecr::registerECROperator(slot = "selectForMating",
                             fun  = ecr::setup(ecr::selTournament,
                                               k = 2L)) %>%
    ecr::registerECROperator(slot = "selectForSurvival",
                             fun  = ecr::selGreedy)

  ### GA
  INITIAL.POP <- create_goelz(N = N, reps = MU, split = TRUE)
  TRIANGLE    <- INITIAL.POP$triangle
  population  <- INITIAL.POP$A.design
  GEN <- 1
  init.fitness <- fitness <- ecr::evaluateFitness(control  = control,
                                                  inds     = population,
                                                  triangle = TRIANGLE)
  init.population <- population <- annotate_pop(population = population,
                                                fitness    = fitness,
                                                gens       = rep(GEN, length(population)),
                                                gen.ids    = paste(GEN, 1:length(population), sep = "-"))

  GAout <- compile_pop(population, fitness, triangle = TRIANGLE, GEN = GEN)
  pop.stats <- GAout$stats
  pop.data  <- GAout$data

  for(GEN in 2:MAX.GEN) {
    print(paste("\nSTARTING GENERATION:", GEN))
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
    fitness <- ecr::evaluateFitness(control  = control,
                                    inds     = offspring,
                                    triangle = TRIANGLE)
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
    GAout <- compile_pop(population, fitness, triangle = TRIANGLE, GEN = GEN)
    pop.stats <- dplyr::bind_rows(pop.stats, GAout$stats)
    pop.data  <- dplyr::bind_rows(pop.data, GAout$data)

    elapsedGen = round((proc.time()[3] - startGen) / 60)
    print(paste0("DONE WITH  GENERATION: ", GEN, " (", elapsedGen, " minutes)"))
  }

  return(list(stats = pop.stats, data = pop.data))
}

goelz_fitness <- function(design, triangle) {
  design <- triangle %>%
    dplyr::filter(trident == "A") %>%
    dplyr::mutate(species = design)
  N <- length(unique(triangle$x))
  L <- seq(from = 0, to = 100, length.out = N)
  se <- NULL
  for(i in 1:N) {
    x.se <- design %>%
      dplyr::filter(x == L[i]) %>%
      dplyr::summarize(species = (sum(species == 1) / nrow(.) * 100 - L[i])^2) %>%
      .$species

    y.se <- design %>%
      dplyr::filter(y == L[i]) %>%
      dplyr::summarize(species = (sum(species == 2) / nrow(.) * 100 - L[i])^2) %>%
      .$species

    z.se <- design %>%
      dplyr::filter(z == L[i]) %>%
      dplyr::summarize(species = (sum(species == 3) / nrow(.) * 100 - L[i])^2) %>%
      .$species

    se <- c(se, x.se[!is.nan(x.se)], y.se[!is.nan(y.se)], z.se[!is.nan(z.se)])
  }
  return(sum(se))
}

annotate_pop <- function(population, fitness, gens, gen.ids) {
  for(i in seq_along(population)) {
    attr(population[[i]], "fitness")    <- fitness[, i]
    attr(population[[i]], "generation") <- gens[i]
    attr(population[[i]], "gen.id")     <- gen.ids[i]
  }
  return(population)
}

compile_pop <- function(pop, fit, triangle, GEN) {

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

  bind_triangle <- function(x, tri) dplyr::bind_cols(tri, x)

  data <- pop %>%
    purrr::map(matrix, ncol = 1) %>%
    purrr::map(dplyr::as_tibble) %>%
    purrr::map(setNames, nm = "species") %>%
    purrr::map(bind_triangle, tri = dplyr::filter(triangle, trident == "A")) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(orig.gen = rep(orig.gen$orig.gen, each = nrow(dplyr::filter(triangle, trident == "A")))) %>%
    dplyr::mutate(gen.id = rep(gen.id$gen.id, each = nrow(dplyr::filter(triangle, trident == "A")))) %>%
    dplyr::mutate(gen = GEN) %>%
    dplyr::select(gen, orig.gen, gen.id, trident, trident.id, id, x, y, z, species)

  return(list(stats = stats, data = data))
}

add_one <- function(x) {
  x[x == 3] <- 4
  x[x == 2] <- 3
  x[x == 1] <- 2
  x[x == 4] <- 1
  return(x)
}
