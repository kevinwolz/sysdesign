#' Create a Nelder Fan experimental design
#' @description Creates a Nelder Fan experimental design.
#' @details The Nelder Fan or Nelder Wheel experimental design (Nelder 1962) is an experimental design that systematically
#' varies plant desnity within a single plot. This function creates Nelder Fan Type Ia (Nelder 1962), where the growing area
#' around each plant has a constant shape throughout the design but increases as radius increases.
#' The terminology and calculations used here follow Parrott, Brinks, and Lhotka (2012).
#' @return An object of class "sysd" and class "nelder". This is a list of 3 elements:
#' \itemize{
#'  \item{"plants"}{ - A data frame (tibble) containing one row for each for each plant in the design.}
#'  \item{"plot"}{ - A data frame (tibble) containing plot charateristics.}
#'  \item{"optim"}{ - If \code{even = TRUE}, a data frame (tibble) comparing plot characteristics between the scenario
#'  where theta is even throughout the design and the scenario where theta is uneven.
#'  If \code{even = FALSE}, then \code{FALSE}.}
#' }
#' @param DN Plant density within the last experimental arc (plants ha-1) (i.e. lower extreme of experimental plant density range)
#' @param D1 Plant density within the first experimental arc (plants ha-1) (i.e. upper extreme of experimental plant density range)
#' @param N Number of experimental arcs (i.e. number of densities to be tested within D1 to DN)
#' @param tau The "rectangularity" proportion. "Rectangularity" is the proportional relationship between the arc length between
#'  spokes and the radial length between arcs, where the numerator represents the arc length and the denominator represents
#'  radial distance. This proportion has been referred to as "rectangularity" in the historical literature and it remains constant
#'  throughout the design. Given that the inner and outer borders of the growing space shape surrounding a plant in a Nelder
#'  design are arcs, and that the shape is not truly rectangular or trapezoidal in nature, the term "rectangularity" can be
#'  confusing.
#' @param even Logical indicated whether or not the design should be adjusted so that the angle between spokes goes evenly
#' into \code{max.angle} (i.e. so that there are no spokes that must be removed from the experiment as border spokes).
#' @param max.angle The maximum rotation (in degrees) of the design. If 360, then a full circle design will be created.
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
nelder <- function(DN, D1, N, tau = 1, even = FALSE, max.angle = 360) {

  if(!(is.numeric(D1) & (D1 %% 1 == 0) & length(D1) == 1))      stop("D1 must an integer and of length 1",           call. = FALSE)
  if(!(is.numeric(DN) & (DN %% 1 == 0) & length(DN) == 1))      stop("DN must an integer and of length 1",           call. = FALSE)
  if(!(is.numeric(N) & (N %% 1 == 0) & length(N) == 1 & N > 2)) stop("N must be an integer greater than 2",          call. = FALSE)
  if(D1 <= DN)                                                  stop("D1 must be greater than DN",                   call. = FALSE)
  if(!(is.numeric(tau) & length(tau) == 1))                     stop("tau argument must be numeric and of length 1", call. = FALSE)
  if(!is.logical(even))                                         stop("even argument must be a logical",              call. = FALSE)
  if(max.angle > 360 | max.angle <= 0)                          stop("max.angle argument must be between 0 and 360", call. = FALSE)

  max.angle.rad <- max.angle * pi / 180  # radians

  ## ORIGINAL VALUES WITHOUT EVEN SPOKE OPTIMIZATION
  tau.orig      <- tau
  alpha.orig    <- exp((log(D1) - log(DN)) / (2 * N - 2))
  theta.orig    <- tau.orig * (alpha.orig ^ 0.5 - alpha.orig ^ -0.5) # radians
  n.spokes.orig <- ceiling(max.angle.rad / theta.orig)

  orig.dat <- nelder_calc(alpha         = alpha.orig,
                          theta         = theta.orig,
                          tau           = tau.orig,
                          D1            = D1,
                          N             = N,
                          max.angle     = max.angle,
                          n.spokes      = n.spokes.orig,
                          spoke.borders = TRUE)

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
                          spoke.borders = max.angle != 360)

  if(even) {
    even.spoke.dat <- dplyr::bind_rows(orig.dat$plot, even.dat$plot) %>%
      dplyr::mutate(even = c(FALSE, TRUE)) %>%
      dplyr::select(even, dplyr::everything())

    out <- c(even.dat, list(optim = even.spoke.dat))

  } else {
    out <- c(orig.dat, list(optim = FALSE))
  }

  class(out) <- c(class(out), "nelder", "sysd")
  return(out)
}

#' Create a Nelder Fan experimental design decision table
#' @description Creates a Nelder Fan experimental design decision table.
#' @details This function helps explore Nelder Fan design options and select a design that meets external constraints
#' (e.g. plant or space availability). Function inputs are identical to \code{\link{nelder}}, but inputs of any length are allowed.
#' All possible combinations of inputs are created using \code{expand.grid}, and then each of these cases is passed to
#' \code{\link{nelder}} for evaluation. Inputs and outputs are all combined and returned for evaluation.
#' @return A tibble containing a wide range of traits of the experimental designs.
#' @param DN Plant density within the last experimental arc (plants ha-1) (i.e. lower extreme of experimental plant density range)
#' @param D1 Plant density within the first experimental arc (plants ha-1) (i.e. upper extreme of experimental plant density range)
#' @param N Number of experimental arcs (i.e. number of densities to be tested within D1 to DN)
#' @param tau The "rectangularity" proportion. See \code{\link{nelder}} for details.
#' @param even Logical indicated whether or not the design should be adjusted so that the angle between spokes goes evenly
#' into \code{max.angle} (i.e. so that there are no spokes that must be removed from the experiment as border spokes).
#' @param max.angle The maximum rotation (in degrees) of the design. If 360, then a full circle design will be created.
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
nelder_decision <- function(DN, D1, N, tau = 1, even = FALSE, max.angle = 360) {
  test.cases <- expand.grid(DN = DN, D1 = D1, N = N, tau = tau, even = even, max.angle = max.angle)

  OUT <- dplyr::tibble()
  for(i in 1:nrow(test.cases)) {
    out <- nelder(DN        = test.cases$DN[i],
                  D1        = test.cases$D1[i],
                  N         = test.cases$N[i],
                  tau       = test.cases$tau[i],
                  even      = test.cases$even[i],
                  max.angle = test.cases$max.angle[i])
    OUT <- dplyr::bind_rows(OUT, out$plot)
  }
  return(dplyr::bind_cols(test.cases, OUT))
}


#' Calculate Nelder Fan design
#' @description Calculates Nelder Fan design
#' Used within \code{\link{nelder}}.
#' @return A list containing the plant and plot data.
#' @param alpha Rate of change along the spokes.
#' @param theta The angle between the spokes.
#' @param tau The "rectangularity" proportion
#' @param D1 Plant density within the first experimental arc (plants ha-1) (i.e. upper extreme of experimental plant density range)
#' @param N Number of experimental arcs (i.e. number of densities to be tested within D1 to DN)
#' @param max.angle The maximum rotation (in degrees) of the design. If 360, then a full circle design will be created.
#' @param n.spokes The number of spokes.
#' @param spoke.borders Logical indicating whether or not there are spokes required to be border spokes.
#' @importFrom dplyr %>%
#' @keywords internal
nelder_calc <- function(alpha, theta, tau, D1, N, max.angle, n.spokes, spoke.borders) {

  Co <- (2 - (alpha + alpha ^ -1)) / (2 * (alpha - alpha ^ -1)) * 100
  r0 <- sqrt(20000 / (D1 * theta * (alpha ^ 3 - alpha))) # m

  arc.dat <- dplyr::tibble(arc     = c(0:(N + 1)),
                           r       = r0 * alpha ^ arc,                         # m
                           space   = theta * r ^ 2 * (alpha - alpha ^ -1) / 2, # m2
                           density = 10000 / space)                            # plants ha-1

  arc.dat$space[c(1, N + 2)]   <- NA
  arc.dat$density[c(1, N + 2)] <- NA

  plot.dat <- dplyr::tibble(area           = pi * tail(arc.dat$r, 1) ^ 2 / 10000 * max.angle / 360, # ha
                            arcs           = N + 2,
                            exp.arcs       = N,
                            spokes         = n.spokes,
                            exp.spokes     = ifelse(spoke.borders, n.spokes - 2, n.spokes),
                            r0             = r0,
                            rmax           = max(arc.dat$r),
                            alpha          = alpha,
                            angle          = theta * 180 / pi,
                            rectangularity = tau,
                            non.centrality = Co,
                            min.density    = min(arc.dat$density, na.rm = TRUE),
                            max.density    = D1,
                            plants         = (N + 2) * n.spokes)

  if(max.angle < 360) plot.dat$spokes <- plot.dat$spokes + 1

  plant.dat <- dplyr::tibble(spoke  = rep(1:plot.dat$spokes, each  = plot.dat$arcs),
                             arc    = rep(arc.dat$arc,       times = plot.dat$spokes),
                             theta  = rep(plot.dat$angle * 0:(plot.dat$spokes - 1), each = plot.dat$arcs)) %>%
    dplyr::left_join(arc.dat, by = "arc") %>%
    dplyr::mutate(x = r * cos(theta * pi / 180)) %>%
    dplyr::mutate(y = r * sin(theta * pi / 180))

  ## DESIGNATE EXPERIMENTAL VS BORDER ROWS
  if(spoke.borders) {
    plant.dat <- dplyr::mutate(plant.dat, exp = (theta %in% range(theta)) | (arc %in% range(arc)))
  } else {
    plant.dat <- dplyr::mutate(plant.dat, exp = arc %in% range(arc))
  }
  plant.dat <- dplyr::mutate(plant.dat, exp = factor(exp, levels = c(FALSE, TRUE), labels = c("Experimental", "Border")))

  plot.dat <- plot.dat %>%
    dplyr::mutate(exp.plants = sum(plant.dat$exp == "Experimental")) %>%
    dplyr::mutate(perc.exp   = round(exp.plants / plants * 100))

  return(list(plants = plant.dat, plot = plot.dat))
}
