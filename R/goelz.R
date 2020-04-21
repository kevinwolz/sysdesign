#' Create a Goelz Triangle experimental design
#' @description Creates a Goelz Triangle experimental design.
#' @details The Goelz Triangle experimental design (Goelz 2001) is an experimental design that systematically
#' varies species composition of three plant species within a single plot, while maintaining a constant plant density.
#' This function creates a Goelz Triangle design following several criteria set forth by Goelz (2001):
#' \itemize{
#'  \item{"symmetry"}{ - Symmetry requires that, for every planting spot assigned, the 2 corresponding planting
#'  spots be assigned corresponding species.}
#'  \item{"equality"}{ - Equality merely requires that equal numbers of each species be assigned to each plot.
#'  This will be achieved if the number of planting spots per plot is a multiple of three and if symmetry is imposed.}
#'  \item{"conformity"}{ - Conformity to the intended pattern requires that the species proportion in any subsection
#'  of the triangular plot is close to the expectations.}
#' }
#' This function does NOT robustly impose the conformity criterion. Species are assigned randomly to each position
#' by sampling the three species using weights based on the theoretical probability of each species at that point.
#' For a robust implementation of the conformity criterion that optimizes the rough initial approach of this function
#' using an evolutionary algorithm, use \code{\link{goelz_optim}}.
#' @return A data frame (tibble) of class sysd and class goelz.
#' @param N The number of plants to be on each edge of the design.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @references
#' \itemize{
#'  \item Goelz JCG (2001) Systematic experimental designs for mixed species plantings. Native Plants Journal 2:90–96.
#'  \url{http://npj.uwpress.org/content/2/2/90.short}
#' }
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- goelz()
goelz <- function(N = 35) {

  if(!(is.numeric(N) & length(N) == 1))       stop("N must be numeric and of length 1",    call. = FALSE)
  if(goelz_count(N = N)$remainder.3 != 0)     stop(paste0("A triangle with N = ", N,
                                                          " does not have a number of points divisible by 3. ",
                                                          "Please use a differnet N."), call. = FALSE)
  if(N < 5) stop("N must be greater than 5 to create a useful Goelz Triangle",             call. = FALSE)

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

  threshold <- 1 / 3 * 100
  A <- triangle %>%
    dplyr::filter(x > threshold & y <= threshold) %>%
    dplyr::arrange(y, z) %>%
    dplyr::mutate(zone = "A") %>%
    dplyr::mutate(zone.id = 1:nrow(.))

  B <- triangle %>%
    dplyr::filter(y > threshold & z <= threshold) %>%
    dplyr::arrange(z, x) %>%
    dplyr::mutate(zone = "B") %>%
    dplyr::mutate(zone.id = 1:nrow(.))

  C <- triangle %>%
    dplyr::filter(z > threshold & x <= threshold) %>%
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

  A.design <- triangle %>%
    dplyr::filter(zone == "A")
  out <- A_to_triangle(triangle  = triangle,
                       A.species = init_species(A.design$x, A.design$y, A.design$z)) %>%
    dplyr::select(id, x.pos, y.pos, x.field, y.field, species, zone, zone.id, x, y, z, border)
  class(out) <- c(class(out), "sysd", "goelz")

  return(out)
}

#' Add border rows to a Goelz Triangle experimental design
#' @description Adds border rows to a Goelz Triangle experimental design.
#' @details Goelz (2001) suggests that all Goelz Triangle experimental designs be enclosed within a number of
#' border/buffer rows. These border rows reduce edge effects and thereby increase the number of plants with useable
#' data. The presence of border rows will also increases the number of individuals whose nearest neighbors are all
#' conspecifics. Goelz (2001) suggests that the species in each planting spot in the border rows be assigned as a
#' 50:50 probability of the species in the two adjacent spots towards the interior of the triangle. This is the method
#' used here, with each additional border row being determined successively. The border rows are also held to the same
#' standard of symmetry across the three zones in the triangle.
#' @return A data frame (tibble) of class sysd, goelz, and goelz-border.
#' @param data An object of class goelz.
#' @param n The number of border rows to add on each side of the design.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @references
#' \itemize{
#'  \item Goelz JCG (2001) Systematic experimental designs for mixed species plantings. Native Plants Journal 2:90–96.
#'  \url{http://npj.uwpress.org/content/2/2/90.short}
#' }
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- goelz()
#' dat.border <- goelz_add_border(data = dat, n = 3)
goelz_add_border <- function(data, n) {

  goelz_class_check(data)
  if(!(is.numeric(n) & length(n) == 1)) stop("n must be numeric and of length 1", call. = FALSE)

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

  class(data) <- unique(c(class(data), "sysd", "goelz", "goelz-border"))
  return(data)
}

#' Mirror a Goelz Triangle experimental design
#' @description Mirrors a Goelz Triangle experimental design and combines it with the original design to create a
#' parallelogram design.
#' @details When establishing two or more replicate Goelz Triangle plots, the plots can be stacked next to each other,
#' with species matched, to further reduce edge effects. This function takes a single Goelz Triangle, mirrors it to
#' create a second triangle, and then combines the two together into a parallelogram design. When merging the two
#' designs, the border rows between the two triangles, if there are any, "double up". This may or may not be desired.
#' Border rows between the two triangles is not really necessary for reducing edge effects, since this edge is now on
#' the interior of the parallelogram. However, at least some borders between triangles may be desired to "separate"
#' the two replicates and "maintain" independence. Any approach can be created by specifiying the exact number of
#' border rows to maintain betwen the two triangles via the \code{joining.borders} argument. The default is to maintain
#' only the number of border rows that each triangle originally had (i.e. half the number of borders initially
#' between the two triangles when merging them.)
#' @return A data frame (tibble) of class sysd, goelz, and goelz-mirror.
#' @param data An object of class goelz.
#' @param joining.borders The number of border rows to maintain between the two triangles.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @references
#' \itemize{
#'  \item Goelz JCG (2001) Systematic experimental designs for mixed species plantings. Native Plants Journal 2:90–96.
#'  \url{http://npj.uwpress.org/content/2/2/90.short}
#' }
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- goelz()
#' dat.border <- goelz_add_border(data = dat, n = 3)
#' dat.border.mirror <- goelz_mirror(data = dat.border)
goelz_mirror <- function(data, joining.borders = max(data$border.num, na.rm = TRUE)) {

  goelz_class_check(data)
  if(!(is.numeric(joining.borders) & length(joining.borders) == 1))
    stop("joining.borders must be a numeric vector of length 1", call. = FALSE)

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
    dplyr::mutate(x.pos    = x.pos - (borders.to.remove - 1)) %>%
    dplyr::mutate(y.pos    = y.pos + y.pos.shift) %>%
    dplyr::mutate(x.field  = x.field - borders.to.remove + 0.5) %>%
    dplyr::mutate(y.field  = y.field + (y.pos.shift * sqrt(3) / 2))

  data.combine <- data %>%
    dplyr::mutate(triangle = "A") %>%
    dplyr::bind_rows(data.mirrored) %>% # Combine two triangles
    dplyr::arrange(y.pos, x.pos) %>% # Reassign ids
    dplyr::mutate(id = 1:nrow(.))

  class(data.combine) <- unique(c(class(data.combine), "sysd", "goelz", "goelz-mirror"))
  return(data.combine)
}

#' Create conformity-optimized Goelz Triangle experimental designs
#' @description Creates conformity-optimized Goelz Triangle experimental designs
#' @details While \code{\link{goelz}} creates Goelz Triangle designs that roughly follow the conformity criterion set
#' forth by Goelz (2001), this function optimizes designs for conformity using an evolutionary algorithm. Function
#' parameters other than \code{data} are all controls on the evolutionary algorithm.
#' The \code{\link{ecr}} package is used for the evolutionary algorithm.
#' @return An list containing:
#' \itemize{
#'  \item{"layout"}{ - A data frame (tibble) containing the base goelz design.}
#'  \item{"stats"}{ - A data frame (tibble) containing statistics on each generation in the evolutionary algorithm.}
#'  \item{"data"}{ - A data frame (tibble) containing the actual data (designs) of each Goelz Triangle in each
#'  generation.}
#' }
#' @param data An object of class goelz.
#' @param save.path A character string indicating the path for where to save data after each generation.
#' @param MU The population size.
#' @param LAMBDA The number of offspring to produce in each generation.
#' @param MAX.GEN The number of generations to run the evolutionary algorithm.
#' @param P.RECOMB The probability of two parents to perform crossover.
#' @param RECOMB The crossover probability for each gene.
#' @param P.MUT The probability to apply mutation to a child.
#' @param MUT The mutation probability for each gene.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @references
#' \itemize{
#'  \item Goelz JCG (2001) Systematic experimental designs for mixed species plantings. Native Plants Journal 2:90–96.
#'  \url{http://npj.uwpress.org/content/2/2/90.short}
#' }
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- goelz_optim()
goelz_optim <- function(data,
                        save.path = ".",
                        MU        = 100,
                        LAMBDA    = MU,
                        MAX.GEN   = 150,
                        P.RECOMB  = 1,
                        RECOMB    = 0.1,
                        P.MUT     = 1,
                        MUT       = 0.005) {

  goelz_class_check(data)
  if(!requireNamespace("ecr", quietly = TRUE)) stop("The package 'ecr' is required for goelz_optim().
                                                    Please install and load it", call. = FALSE)

  ### GA CONTROLS
  control <- ecr::initECRControl(fitness.fun  = goelz_fitness,
                                 n.objectives = 1,
                                 minimize     = TRUE) %>%
    ecr::registerECROperator(slot = "mutate",
                             fun  = ecr::setup(ecr::makeMutator(mutInteger, supported = "float"),
                                               p     = MUT,
                                               lower = 1L,
                                               upper = 3L))

  ## INITIAL POPULATION
  TRIANGLE <- data %>%
    dplyr::select(id, x.pos, y.pos, x.field, y.field, zone, zone.id, x, y, z, border)

  INITIAL.POP <- list(data)
  for(i in 1:(MU - 1)) {
    g <- list(goelz(N = max(data$x.pos)))
    INITIAL.POP <- c(INITIAL.POP, g)
  }

  pull.A.design <- function(triangle) {
    triangle %>%
      dplyr::filter(zone == "A") %>%
      .$species
  }

  population  <- purrr::map(INITIAL.POP, pull.A.design)

  ## GA
  ga.out <- run_ga(control     = control,
                   population  = population,
                   layout      = dplyr::filter(TRIANGLE, zone == "A"),
                   layout.comp = dplyr::filter(TRIANGLE, zone == "A"),
                   LAMBDA      = LAMBDA,
                   MAX.GEN     = MAX.GEN,
                   P.RECOMB    = P.RECOMB,
                   RECOMB      = RECOMB,
                   P.MUT       = P.MUT,
                   save.path   = save.path)

  OUT <- list(layout = TRIANGLE, stats = ga.out$pop.stats, data = ga.out$pop.data)
  class(OUT) <- c(class(OUT), "goelz-optim")
  return(OUT)
}

#' Select optimal Goelz Triangle design sfrom the results of goelz_optim
#' @description Selects optimal Goelz Triangle designs from the results of \code{\link{goelz_optim}}.
#' @return A list of objects of class goelz.
#' @param optim.results An object of class goelz-optim.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' optim.dat <- goelz_optim()
#' my_goelz <- select_optimal_goelz(optim.results = optim.dat)
select_optimal_goelz <- function(optim.results) {

  if(!("goelz-optim" %in% class(optim.results))) stop("optim.results must be of class goelz-optim", call. = FALSE)

  N <- max(optim.results$layout$x.pos)
  best.designs <- select_optimal(optim.results)

  OUT <- list()
  for(i in 1:length(best.designs)) {
    out <- goelz(N = N) %>%
      dplyr::arrange(zone, zone.id) %>%
      A_to_triangle(A.species = best.designs[[i]]$species)
    OUT <- c(OUT, list(out))
  }

  return(OUT)
}

#' Get coordinates of corners of Goelz Triangle experimental design
#' @description Gets coordinates of corners of Goelz Triangle experimental design.
#' @return A data frame (tibble).
#' @param data An object of class goelz.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- goelz()
#' goelz_corners(data = dat)
goelz_corners <- function(data) {

  goelz_class_check(data)

  x.range <- range(round(data$x.field, 2))
  y.range <- range(round(data$y.field, 2))

  out <- expand.grid(x.field = x.range, y.field = y.range) %>%
    dplyr::arrange(x.field, y.field) %>%
    dplyr::as_tibble()

  return(out)
}

#' Get coordinates of guides for Goelz Triangle experimental design
#' @description Gets coordinates of guides for Goelz Triangle experimental design.
#' @return A data frame (tibble).
#' @param data An object of class goelz.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- goelz()
#' goelz_guides(data = dat)
goelz_guides <- function(data) {

  goelz_class_check(data)

  x.range  <- range(round(data$x.field, 2))
  x.range  <- c(x.range, x.range[1] + diff(x.range) / 3, x.range[2] - diff(x.range) / 3)
  y.unique <- unique(round(data$y.field, 2))

  out <- expand.grid(x.field = x.range, y.field = y.unique) %>%
    dplyr::arrange(x.field, y.field)

  return(out)
}

#' Get coordinates of row starts for Goelz Triangle experimental design
#' @description Gets coordinates of row starts for Goelz Triangle experimental design.
#' @return A data frame (tibble).
#' @param data An object of class goelz.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @importFrom dplyr %>%
#' @family definition functions
#' @examples
#' dat <- goelz()
#' goelz_starts(data = dat)
goelz_starts <- function(data) {

  goelz_class_check(data)

  out <- data %>%
    dplyr::filter(x.pos == 1) %>%
    dplyr::select(x.pos, y.pos, x.field, y.field) %>%
    dplyr::mutate(y.field = round(y.field, 2))

  return(out)
}

#' Calculate fitness of Goelz Trianlge design
#' @description Calculates fitness of Goelz Trianlge design.
#' Used by \code{\link{goelz_optim}} as the fitness function.
#' @return A numeric value.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @importFrom dplyr %>%
#' @keywords internal
goelz_fitness <- function(design, layout) {

  design <- layout %>%
    dplyr::mutate(species = design)
  N <- length(unique(layout$x))
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

#' Add one
#' @description Adds one.
#' Used by \code{\link{goelz}}, \code{\link{goelz_add_border}}, and \code{\link{A_to_triangle}}.
#' @return A numeric vector.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @keywords internal
add_one <- function(x) {
  x[x == 3] <- 4
  x[x == 2] <- 3
  x[x == 1] <- 2
  x[x == 4] <- 1
  return(x)
}

#' Apply symmetry from one trianlge region across other two regions.
#' @description Applies symmetry from one trianlge region across other two regions.
#' Used by \code{\link{goelz}}.
#' @return A list.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @importFrom dplyr %>%
#' @keywords internal
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

  return(dplyr::bind_rows(A.design, B.design, C.design))
}

#' Remove edge
#' @description Removes edge.
#' Used by \code{\link{goelz_mirror}}.
#' @return An object of class goelz.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @importFrom dplyr %>%
#' @keywords internal
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

#' Count total plants in Goelz Triangle design.
#' @description Counts total plants in Goelz Triangle design. Used by \code{\link{goelz}}.
#' @return A data frame (tibble).
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @keywords internal
goelz_count <- function(N) {
  count_func <- function(x) sum(1:x)
  out <- dplyr::tibble(N            = N,
                       total.points = purrr::map_dbl(N, count_func),
                       remainder.3  = total.points %% 3)
  return(out)
}

#' Check if an object if of class goelz
#' @description Checks if an object if of class goelz. Used by many goelz definition functions.
#' @return An error.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @keywords internal
goelz_class_check <- function(data) {
  if(!("goelz" %in% class(data))) {
    if("goelz-mirror" %in% class(data)) {
      stop("data must be of class 'goelz', not 'goelz-mirror'", call. = FALSE)
    } else stop("data must be of class 'goelz'", call. = FALSE)
  }
}
