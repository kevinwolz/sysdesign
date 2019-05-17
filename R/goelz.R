
goelz <- function(N     = 35,
                  reps  = 1,
                  split = FALSE) {

  L <- seq(from = 0, to = 100, length.out = N)
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
    dplyr::mutate(x.field = x.pos + 0.5 * (y.pos - 1) - 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(border = FALSE)

  A <- triangle %>%
    dplyr::filter(x > 35 & y <= 35) %>%
    dplyr::arrange(y, z) %>%
    dplyr::mutate(zone = "A") %>%
    dplyr::mutate(zone.id = 1:nrow(.))

  B <- triangle %>%
    dplyr::filter(y > 35 & z <= 35) %>%
    dplyr::arrange(z, x) %>%
    dplyr::mutate(zone = "B") %>%
    dplyr::mutate(zone.id = 1:nrow(.))

  C <- triangle %>%
    dplyr::filter(z > 35 & x <= 35) %>%
    dplyr::arrange(x, y) %>%
    dplyr::mutate(zone = "C") %>%
    dplyr::mutate(zone.id = 1:nrow(.))

  triangle <- dplyr::bind_rows(A, B, C) %>%
    dplyr::mutate(zone = factor(zone, levels = c("A", "B", "C"))) %>%
    dplyr::arrange(y.pos, x.pos) %>%
    dplyr::select(id, x.pos, y.pos, x.field, y.field, zone, zone.id, x, y, z, border)

  init_species <- function(x, y, z) {
    OUT <- NULL
    for(i in 1:length(x)) {
      out <- sample(SPECIES, size = 1, replace = TRUE, prob = c(x[i], y[i], z[i]))
      OUT <- c(OUT, out)
    }
    return(OUT)
  }

  out <- OUT <- A.out <- list()
  index <- 1
  if(split) {
    while(index <= reps) {
      A.design <- triangle %>%
        dplyr::filter(zone == "A") %>%
        dplyr::mutate(species = init_species(x, y, z)) %>%
        .$species

      B.design <- add_one(A.design)
      C.design <- add_one(B.design)

      out <- c(out, list(c(A.design, B.design, C.design)))
      A.out <- c(A.out, list(A.design))
      index <- index + 1
    }
    OUT <- list(triangle = triangle, design = out, A.design = A.out)
  } else {
    while(index <= reps) {
      A.design <- triangle %>%
        dplyr::filter(zone == "A")
      out <- A_to_triangle(triangle  = triangle,
                           A.species = init_species(A.design$x, A.design$y, A.design$z))
      out[[1]] <- out[[1]] %>% dplyr::select(id, x.pos, y.pos, x.field, y.field, species, zone, zone.id, x, y, z, border)
      OUT <- c(OUT, out)
      index <- index + 1
    }
  }

  g.class <- ifelse(split, "goelz-split", "goelz")
  class(OUT) <- c(class(OUT), g.class, "sysd")
  return(OUT)
}

goelz_add_border <- function(data, n) {

  resample <- function(x, ...) x[sample.int(length(x), ...)]

  data <- data %>%
    dplyr::mutate(zone = as.character(zone))

  ## CREATE BORDER LOCATIONS
  border.max <- 0
  border.row <- 1
  while(border.row <= n) {
    bottom.start <- data %>%
      dplyr::filter(y.pos == min(y.pos)) %>%
      dplyr::select(x.pos, y.pos, x.field, y.field, zone, species) %>%
      dplyr::mutate(zone    = "border-bottom") %>%
      dplyr::mutate(x.field = x.field - 1.5) %>%
      dplyr::mutate(y.field = y.field - (sqrt(3) / 2)) %>%
      dplyr::mutate(x.pos   = x.pos   - 1) %>%
      dplyr::mutate(y.pos   = y.pos   - 1)

    bottom.end <- bottom.start %>%
      dplyr::arrange(x.pos) %>%
      tail(2) %>%
      dplyr::mutate(x.field = x.field + 2) %>%
      dplyr::mutate(x.pos   = x.pos   + 2)

    bottom.full <- dplyr::bind_rows(bottom.start, bottom.end) %>%
      dplyr::arrange(x.pos) %>%
      dplyr::mutate(border     = TRUE) %>%
      dplyr::mutate(border.num = border.row) %>%
      dplyr::mutate(zone.id    = 1:nrow(.) + border.max) %>%
      dplyr::arrange(zone.id)

    for(i in 1:nrow(bottom.full)) {
      x.search <- bottom.full$x.pos[i] + c(-1, 0)
      y.search <- bottom.full$y.pos[i] + 1

      search.window <- data %>%
        dplyr::filter(x.pos %in% x.search & y.pos %in% y.search)

      if(nrow(search.window) != 0) {
        bottom.full$species[i] <- resample(search.window$species, 1)
      }
    }

    left.start <- data %>%
      dplyr::filter(x.pos == min(x.pos)) %>%
      dplyr::select(x.pos, y.pos, x.field, y.field, zone, species) %>%
      dplyr::mutate(zone    = "border-left") %>%
      dplyr::mutate(x.field = x.field - 1) %>%
      dplyr::mutate(x.pos   = x.pos   - 1)

    left.end <- left.start %>%
      dplyr::arrange(dplyr::desc(y.pos)) %>%
      head(2) %>%
      dplyr::mutate(x.field = x.field + 1.0) %>%
      dplyr::mutate(y.field = y.field + 2 * (sqrt(3) / 2)) %>%
      dplyr::mutate(y.pos   = y.pos   + 2)

    left.full <- dplyr::bind_rows(left.start, left.end) %>%
      dplyr::arrange(dplyr::desc(y.pos)) %>%
      dplyr::mutate(border     = TRUE) %>%
      dplyr::mutate(border.num = border.row) %>%
      dplyr::mutate(zone.id    = 1:nrow(.) + border.max) %>%
      dplyr::arrange(zone.id)

    left.full$species <- add_one(bottom.full$species)

    right.start <- data %>%
      dplyr::group_by(y.pos) %>%
      dplyr::filter(x.pos == max(x.pos)) %>%
      dplyr::ungroup() %>%
      dplyr::select(x.pos, y.pos, x.field, y.field, zone, species) %>%
      dplyr::mutate(zone    = "border-right") %>%
      dplyr::mutate(x.field = x.field + 0.5) %>%
      dplyr::mutate(y.field = y.field + (sqrt(3) / 2)) %>%
      dplyr::mutate(y.pos   = y.pos   + 1)

    right.end <- right.start %>%
      dplyr::arrange(y.pos) %>%
      head(2) %>%
      dplyr::mutate(x.field = x.field + 1) %>%
      dplyr::mutate(x.pos   = x.pos   + 2) %>%
      dplyr::mutate(y.field = y.field - 2 * (sqrt(3) / 2)) %>%
      dplyr::mutate(y.pos   = y.pos   - 2)

    right.full <- dplyr::bind_rows(right.start, right.end) %>%
      dplyr::arrange(y.pos) %>%
      dplyr::mutate(border     = TRUE) %>%
      dplyr::mutate(border.num = border.row) %>%
      dplyr::mutate(zone.id    = 1:nrow(.) + border.max) %>%
      dplyr::arrange(zone.id)

    right.full$species <- add_one(left.full$species)

    data <- dplyr::bind_rows(data,
                              bottom.full,
                              left.full,
                              right.full)

    border.row <- border.row + 1
    border.max <- border.max + nrow(bottom.full)
  }

  data <- data %>%
    dplyr::mutate(x.pos = x.pos + n) %>%
    dplyr::mutate(y.pos = y.pos + n) %>%
    dplyr::mutate(x.field = x.field + abs(min(x.field))) %>%
    dplyr::mutate(y.field = y.field + abs(min(y.field))) %>%
    dplyr::mutate(zone = factor(zone, levels = c("A", "B", "C", "border-bottom", "border-left", "border-right"))) %>%
    dplyr::arrange(y.pos, x.pos)

  class(data) <- unique(c(class(data), "goelz", "sysd"))
  return(data)
}

goelz_optim <- function(N        = 35,
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
  INITIAL.POP <- goelz(N = N, reps = MU, split = TRUE)
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

goelz_mirror <- function(data, joining.borders = max(data$border.num, na.rm = TRUE)) {

  # Create & rotate second triangle
  data.mirrored <- data

  theta <- -pi / 3
  x.field.new <- data$x.field * cos(theta) - data$y.field * sin(theta) + max(data$x.field) / 2 + 1
  x.field.new <- x.field.new - 2 * (x.field.new - mean(x.field.new))
  y.field.new <- data$y.field * cos(theta) + data$x.field * sin(theta) + max(data$y.field)

  data.mirrored$x.pos   <- -data$y.pos + max(data$x.pos) + 2
  data.mirrored$y.pos   <- -data$x.pos + max(data$y.pos) + 1
  data.mirrored$x.field <- x.field.new
  data.mirrored$y.field <- y.field.new

  # Modify border rows if any
  n.borders <- max(data$border.num, na.rm = TRUE)
  borders.to.remove <- n.borders * 2 - joining.borders

  if(borders.to.remove > 0) {
    if(borders.to.remove %% 2 == 0) { # borders.to.remove is even
      data          <- remove_edge(data = data,          side = "right", n = borders.to.remove / 2)
      data.mirrored <- remove_edge(data = data.mirrored, side = "left",  n = borders.to.remove / 2)
      y.pos.shift <- 0
    } else { # borders.to.remove is odd
      data          <- remove_edge(data = data,          side = "right", n = floor(borders.to.remove / 2))
      data.mirrored <- remove_edge(data = data.mirrored, side = "left",  n = floor(borders.to.remove / 2) + 1)
      y.pos.shift <- -1
    }
  }

  data.mirrored <- data.mirrored %>%
    dplyr::mutate(triangle = "B") %>%
    dplyr::mutate(x.pos   = x.pos - (borders.to.remove - 1)) %>%
    dplyr::mutate(y.pos   = y.pos + y.pos.shift) %>%
    dplyr::mutate(x.field = x.field - borders.to.remove + 0.5) %>%
    dplyr::mutate(y.field = y.field + (y.pos.shift * sqrt(3) / 2))

  data.combine <- data %>%
    dplyr::mutate(triangle = "A") %>%
    dplyr::bind_rows(data.mirrored) %>% # Combine two triangles
    dplyr::arrange(y.pos, x.pos) %>% # Reassign ids
    dplyr::mutate(id = 1:nrow(.))

  class(data.combine) <- unique(c(class(data.combine), "goelz", "sysd"))
  return(data.combine)
}

goelz_corners <- function(data) {
  x.range <- range(round(data$x.field, 2))
  y.range <- range(round(data$y.field, 2))
  out <- expand.grid(x.field = x.range, y.field = y.range) %>%
    dplyr::arrange(x.field, y.field)
  return(out)
}

goelz_guides <- function(data) {
  x.range  <- range(round(data$x.field, 2))
  x.range <- c(x.range, x.range[1] + diff(x.range) / 3, x.range[2] - diff(x.range) / 3)
  y.unique <- unique(round(data$y.field, 2))
  out <- expand.grid(x.field = x.range, y.field = y.unique) %>%
    dplyr::arrange(x.field, y.field)
  return(out)
}

goelz_starts <- function(data) {
  out <- data %>%
    dplyr::filter(x.pos == 1) %>%
    dplyr::select(x.pos, y.pos, x.field, y.field) %>%
    dplyr::mutate(y.field = round(y.field, 2))
  return(out)
}

goelz_fitness <- function(design, triangle) {
  design <- triangle %>%
    dplyr::filter(zone == "A") %>%
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
    purrr::map(bind_triangle, tri = dplyr::filter(triangle, zone == "A")) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(orig.gen = rep(orig.gen$orig.gen, each = nrow(dplyr::filter(triangle, zone == "A")))) %>%
    dplyr::mutate(gen.id = rep(gen.id$gen.id, each = nrow(dplyr::filter(triangle, zone == "A")))) %>%
    dplyr::mutate(gen = GEN) %>%
    dplyr::select(gen, orig.gen, gen.id, zone, zone.id, id, x, y, z, species)

  return(list(stats = stats, data = data))
}

add_one <- function(x) {
  x[x == 3] <- 4
  x[x == 2] <- 3
  x[x == 1] <- 2
  x[x == 4] <- 1
  return(x)
}

A_to_triangle <- function(triangle, A.species) {
  A.design <- triangle %>%
    dplyr::filter(zone == "A") %>%
    dplyr::arrange(zone.id) %>%
    dplyr::mutate(species = A.species)

  B.design <- triangle %>%
    dplyr::filter(zone == "B") %>%
    dplyr::arrange(zone.id) %>%
    dplyr::mutate(species = add_one(A.design$species))

  C.design <- triangle %>%
    dplyr::filter(zone == "C") %>%
    dplyr::arrange(zone.id) %>%
    dplyr::mutate(species = add_one(B.design$species))

  return(list(dplyr::bind_rows(A.design, B.design, C.design)))
}

remove_edge <- function(data, side, n) {
  if(side == "right") func <- max else func <- min
  for(i in 1:n) {
    data <- data %>%
      dplyr::group_by(y.pos) %>%
      dplyr::filter(x.pos != func(x.pos)) %>%
      dplyr::ungroup()
  }
  return(data)
}
