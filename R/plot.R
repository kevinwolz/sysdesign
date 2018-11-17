#' Plot a Nelder Fan experimental design
#' @description Plots a Nelder Fan experimental design.
#' @return A ggplot2 object.
#' @param data An object of class 'nelder' created via \code{\link{nelder}}.
#' @param size A numeric value indicating point size.
#' @param fill.colors A character vector of length 2 indicating the fill colors of experimental and border plants, respectively.
#' @param legend Logical indicating whether or not to include a legend for fill colors.
#' @param ex.area Logical indicating whether or not and example growing area of one plant should be shown.
#' @param caption Logical indicating whether or not to include a caption with generic data about the design.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @import ggplot2
#' @family plot functions
#' @examples
#' data <- nelder(DN         = 1000,
#'                D1         = 3000,
#'                N          = 5,
#'                tau        = 1,
#'                even       = TRUE,
#'                max.angle  = 360)
#' plot_nelder(data)
plot_nelder <- function(data,
                        size        = 3,
                        fill.colors = c("black", "white"),
                        legend      = FALSE,
                        ex.area     = FALSE,
                        caption     = FALSE) {

  if(!("nelder" %in% class(data)))                            stop("data must be of class 'nelder'",                     call. = FALSE)
  if(!(is.numeric(size) & length(size) == 1))                 stop("size must be numeric and of length 1",               call. = FALSE)
  if(!(is.character(fill.colors) & length(fill.colors) == 2)) stop("fill.colors must be a character vector of length 2", call. = FALSE)
  if(!is.logical(ex.area))                                    stop("ex.area must be a logical",                          call. = FALSE)
  if(!is.logical(caption))                                    stop("caption must be a logical",                          call. = FALSE)

  plot.caption <- NULL
  if(caption) {
    plot.caption <- paste0("# Spokes: ",        data$plot$spokes, " (", data$plot$exp.spokes, ")",
                           "\n# Arcs: ",        data$plot$arcs,   " (", data$plot$exp.arcs, ")",
                           "\n# Angle: ",       round(data$plot$angle, 2),  " degrees",
                           "\n# Plants: ",      data$plot$plants, " (", data$plot$exp.plants, ")",
                           "\nDensity range: ", round(data$plot$min.density), " - ", data$plot$max.density, " plants/ha",
                           "\nr0-rmax: ",       round(data$plot$r0, 2),       " - ", round(data$plot$rmax, 2), " m",
                           "\nPlot area: ",     round(data$plot$area, 2), " ha")
  }

  area.geom <- geom_blank()
  if(ex.area) {
    r <- unique(data$plants$r)
    i <- length(r) - 1
    radius.max <- (r[i] + r[i + 1]) / 2
    radius.min <- (r[i] + r[i - 1]) / 2

    area.shade <- dplyr::tibble(ymin = rep(data$plot$angle - data$plot$angle / 2, 2),
                                ymax = rep(data$plot$angle + data$plot$angle / 2, 2),
                                x    = c(radius.min, radius.max))

    area.geom <- geom_ribbon(aes(ymin = ymin, ymax = ymax, x = x),
                             data = area.shade,
                             fill = "grey50")
  }

  plot.obj <- ggplot(data$plants) +
    area.geom +
    geom_point(aes(x = r, y = theta, fill = exp), shape = 21, size = size) +
    scale_fill_manual(values = fill.colors) +
    lims(x = c(0, NA),
         y = c(0, 360)) +
    coord_polar(theta = "y", start = -pi / 2, direction = -1) +
    theme_void() +
    labs(fill = NULL, caption = plot.caption) +
    theme(plot.caption    = element_text(size = 20, hjust = 0),
          legend.position = ifelse(legend, "bottom", "none"),
          legend.text     = element_text(size = 18))

  return(plot.obj)
}
