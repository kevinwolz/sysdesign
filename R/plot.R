#' Plot a Nelder Fan experimental design
#' @description Plots a Nelder Fan experimental design.
#' @return A ggplot2 object.
#' @param data An object of class 'nelder' created via \code{\link{nelder}}.
#' @param size A numeric value indicating point size.
#' @param fill.palette A character vector of length 2 indicating the fill colors of different species.
#' @param species Logical indicating whether or not to plot species as differnet colors
#' (see \code{\link{nelder_biculture}}).
#' @param legend Logical indicating whether or not to include a legend for fill colors.
#' @param ex.area Logical indicating whether or not and example growing area of one plant should be shown.
#' @param caption Logical indicating whether or not to include a caption with generic data about the design.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @import ggplot2
#' @family plot functions
#' @examples
#' dat <- nelder(DN         = 1000,
#'               D1         = 3000,
#'               N          = 5,
#'               tau        = 1,
#'               even       = TRUE,
#'               max.angle  = 360)
#' plot_nelder(data = dat)
plot_nelder <- function(data,
                        size         = 3,
                        fill.palette = c("black", "white"),
                        species      = TRUE,
                        legend       = FALSE,
                        ex.area      = FALSE,
                        caption      = FALSE) {

  nelder_class_check(data)
  if(!(is.numeric(size) & length(size) == 1)) stop("size must be numeric and of length 1",                call. = FALSE)
  if(!(is.character(fill.palette) &
       length(fill.palette) == 2))            stop("fill.palette must be a character vector of length 2", call. = FALSE)
  if(!is.logical(ex.area))                    stop("ex.area must be a logical",                           call. = FALSE)
  if(!is.logical(caption))                    stop("caption must be a logical",                           call. = FALSE)
  if(!is.logical(legend))                     stop("legend must be a logical",                            call. = FALSE)

  plot.caption <- NULL
  if(caption) {
    plot.caption <- paste0("# Spokes: ",        data$plot$spokes, " (", data$plot$exp.spokes, ")",
                           "\n# Arcs: ",        data$plot$arcs,   " (", data$plot$exp.arcs, ")",
                           "\n# Angle: ",       round(data$plot$angle, 2),  " degrees",
                           "\n# Plants: ",      data$plot$plants, " (", data$plot$exp.plants, ")",
                           "\nDensity range: ", round(data$plot$min.exp.density), " - ", data$plot$max.exp.density, " plants/ha",
                           "\nr0-rmax: ",       round(data$plot$rmin, 2), " - ", round(data$plot$rmax, 2), " m",
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
    lims(x = c(0, NA),
         y = c(0, 360)) +
    coord_polar(theta = "y", start = -pi / 2, direction = -1) +
    theme_void() +
    labs(shape = NULL, caption = plot.caption) +
    guides(fill = FALSE) +
    theme(plot.caption    = element_text(size = 20, hjust = 0),
          legend.position = ifelse(legend, "bottom", "none"),
          legend.text     = element_text(size = 18)) +
    scale_fill_manual(values  = fill.palette) +
    scale_shape_manual(values = c(21, 24))

  if(species & "species" %in% names(data$plants)) {
    plot.obj <- plot.obj +
      geom_point(aes(x = r, y = theta, fill = species, shape = border), size = size)

  } else {
    plot.obj <- plot.obj +
      geom_point(aes(x = r, y = theta, shape = border), size = size, fill = "black")
  }

  return(plot.obj)
}

#' Plot a Goelz Triangle experimental design
#' @description Plots a Goelz Triangle experimental design.
#' @return A ggplot2 object.
#' @param data An object of class 'goelz' created via \code{\link{goelz}} with \code{split} = FALSE.
#' @param size A numeric value indicating point size.
#' @param fill One of "species", "zone", or "border", indicating the column of \code{data} that specifies the point fill color.
#' @param color One of "species", "zone", or "border", indicating the column of \code{data} that specifies the point border color.
#' Use "none" for no variation in point borders.
#' @param label A character string indicating the column name of \code{data} to use for labels on each point.
#' Use "none" for no labels.
#' @param fill.discrete.palette A character vector indicating the color palette of used for point fill.
#' @param color.discrete.palette A character vector indicating the color palette of used for point borders.
#' @param fill.continuous.palette ggplot2 function for a continuous color palette used for point fill.
#' @param color.continuous.palette ggplot2 function for a continuous color palette used for point color
#' @param corners Logical indicating whether or not to plot red X's at each corner of the rectangle containing the design.
#' @param guides Logical indicating whether or not to plot red +'s as guides at the extremes of each row in the design.
#' @param legend Logical indicating whether or not to include a legend for fill and color values.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @import ggplot2
#' @family plot functions
#' @examples
#' dat <- goelz(N = 35, reps = 1, split = FALSE)
#' plot_goelz(data = dat)
plot_goelz <- function(data,
                       size    = 4,
                       fill    = "species",
                       color   = "border",
                       label   = "none",
                       fill.discrete.palette    = c("grey50", "white", "#E69F00", "#56B4E9",
                                                    "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                       color.discrete.palette   = c("black", "green", "red", "blue", "orange", "purple"),
                       fill.continuous.palette  = ggplot2::scale_fill_viridis_c(option = "magma"),
                       color.continuous.palette = ggplot2::scale_color_viridis_c(option = "magma"),
                       corners = FALSE,
                       guides  = FALSE,
                       legend  = TRUE){

  goelz_class_check(data)
  if(!(is.numeric(size) & length(size) == 1))     stop("size must be numeric and of length 1",         call. = FALSE)
  if(!(is.character(label) & length(label) == 1)) stop("label must be a character vector of length 1", call. = FALSE)
  if(!is.character(fill.discrete.palette))        stop("fill.palette must be a character vector",      call. = FALSE)
  if(!is.character(color.discrete.palette))       stop("color.palette must be a character vector",     call. = FALSE)
  if(!is.logical(corners))                        stop("corners must be a logical",                    call. = FALSE)
  if(!is.logical(guides))                         stop("guides must be a logical",                     call. = FALSE)
  if(!(fill  %in% names(data)))                   stop("fill must be a column name in data",            call. = FALSE)
  if(!(color %in% c("none", names(data))))        stop("color must be a column name in data or 'none'", call. = FALSE)

  data$species <- factor(data$species)

  plot.obj <- ggplot(data, aes(x = x.field, y = y.field)) +
    theme_void() +
    coord_equal() +
    theme(legend.position = ifelse(legend, "right", "none"))

  if(is.numeric(data[[fill]])) {
    plot.obj <- plot.obj +
      fill.continuous.palette
  } else {
    plot.obj <- plot.obj +
      scale_fill_manual(values = fill.discrete.palette)
  }

  if(color == "none") {
    plot.obj <- plot.obj +
      geom_point(size = size, shape = 21, aes_string(fill = fill), color = color.discrete.palette[1])
  } else if(is.numeric(data[[color]])) {
    plot.obj <- plot.obj +
      geom_point(size = size, shape = 21, aes_string(fill = fill, color = color)) +
      color.continuous.palette
  } else {
    plot.obj <- plot.obj +
      geom_point(size = size, shape = 21, aes_string(fill = fill, color = color)) +
      scale_color_manual(values = color.discrete.palette)
  }

  if(label != "none") plot.obj <- plot.obj + geom_text(aes_string(label = label))
  if(corners)         plot.obj <- plot.obj + geom_point(data = goelz_corners(data), color = "red", shape = 4, size = size / 2)
  if(guides)          plot.obj <- plot.obj + geom_point(data = goelz_guides(data),  color = "red", shape = 3, size = size / 4)

  return(plot.obj)
}

#' Plot a diagnositc display of compeitition in a Nelder Fan biculture experimental design
#' @description Plots a diagnositc display of compeitition in a Nelder Fan biculture experimental design.
#' The optimal competition scenario is to have counts consistent across arcs. See \code{\link{nelder_optim}} for more information.
#' @return If \code{plot = TRUE}, returns a ggplot object, otherwise the data that would create the plot is returned.
#' @param data An object of class 'nelder-biculture' created via \code{\link{goelz}} with \code{split} = FALSE.
#' @param plot If \code{TRUE}, the default, a ggplot object is returned.
#' If \code{FALSE}, the data that would create the plot is returned.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @import ggplot2
#' @importFrom dplyr %>%
#' @family plot functions
#' @examples
#' dat <- nelder()
#' dat.bi <- nelder_biculture(data = dat)
#' plot_nelder_biculture_competition(data = dat.bi)
plot_nelder_biculture_competition <- function(data, plot = TRUE) {

  nelder_biculture_class_check(data)

  plot.data <- data %>%
    nelder_biculture_competition() %>%
    dplyr::mutate(common.neighbors = ifelse(species == "A", A.neighbors, B.neighbors)) %>%
    dplyr::group_by(species, arc, common.neighbors) %>%
    dplyr::filter(!is.na(common.neighbors)) %>%
    dplyr::summarize(dens = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(species = factor(species,
                                   levels = c("A", "B"),
                                   labels = c("Species A", "Species B")))

  base_size <- 18
  top.dens <- plot.data$dens %>%
    unique() %>%
    sort() %>%
    rev() %>%
    .[1:ceiling(length(.) * 0.33)]

  plot.obj <- ggplot(plot.data, aes(x     = arc,
                                    y     = common.neighbors,
                                    fill  = dens,
                                    label = dens)) +
    labs(x    = "Arc #",
         y    = "Number of same-\nspecies neighbors",
         fill = "Number of\nindividuals") +
    scale_x_continuous(breaks = unique(plot.data$arc), expand = c(0,0)) +
    scale_y_continuous(breaks = 0:8, expand = c(0,0)) +
    geom_raster(na.rm = TRUE) +
    geom_text(aes(color = dens %in% top.dens)) +
    facet_wrap(~species) +
    scale_fill_viridis_c(option = "magma") +
    scale_color_manual(values = c("white", "black")) +
    guides(color = FALSE) +
    coord_equal() +
    theme_bw(base_size = base_size) +
    theme(plot.margin       = unit(base_size * c(1,1,1,1), "points"),
          panel.border      = element_rect(size = 2, color = "black"),
          axis.text         = element_text(color = "black"),
          axis.title.x      = element_text(vjust = -1),
          axis.title.y      = element_text(vjust = 2),
          legend.key.width  = unit(1.5, "cm"),
          panel.grid.major  = element_blank(),
          axis.ticks        = element_line(color = "black"),
          axis.ticks.length = unit(base_size * 0.25, "points"))

  return(if(plot) plot.obj else plot.data)
}

#' Plot a diagnositc display of the fitness trajectory
#' @description Plots a diagnositc display of the fitness trajectory from the genetic algorithm run by
#' \code{\link{goelz_optim}} or \code{\link{nelder_biculture_optim}}.
#' @return If \code{plot = TRUE}, returns a ggplot object, otherwise the data that would create the plot is returned.
#' @param data An object of class 'goelz-optim' created via \code{\link{goelz_optim}} or
#' 'nelder-optim' created via \code{\link{nelder_biculture_optim}}.
#' @param plot If \code{TRUE}, the default, a ggplot object is returned.
#' If \code{FALSE}, the data that would create the plot is returned.
#' @author Kevin J Wolz, \email{kevin@@savannainstitute.org}
#' @export
#' @import ggplot2
#' @importFrom dplyr %>%
#' @family plot functions
#' @examples
#' dat <- nelder()
#' dat.bi <- nelder_biculture(data = dat)
#' dat.bi.optim <- nelder_biculture_optimr(data = dat.bi)
#' plot_fitness_trajectory(data = dat.bi.optim)
plot_fitness_trajectory <- function(data, plot = TRUE) {

  if(!any(c("nelder-optim", "goelz-optim") %in% class(data))) {
    stop("data must be of class nelder-optim or goelz-optim", call. = FALSE)
  }

    plot.data <- data$stats %>%
      dplyr::group_by(gen) %>%
      dplyr::summarize(fitness.mean = mean(fitness),
                       fitness.min  = min(fitness),
                       fitness.max  = max(fitness))

    base_size <- 18

    plot.obj <- ggplot(plot.data, aes(x = gen)) +
      labs(x = "Generation",
           y = "Fitness") +
      geom_line(aes(y = fitness.mean), size = 2) +
      geom_line(aes(y = fitness.min), linetype = "dashed") +
      geom_line(aes(y = fitness.max), linetype = "dashed") +
      theme_bw(base_size = base_size) +
      theme(plot.margin       = unit(base_size * c(1,1,1,1), "points"),
            panel.border      = element_rect(size = 2, color = "black"),
            axis.text         = element_text(color = "black"),
            axis.title.x      = element_text(vjust = -1),
            axis.title.y      = element_text(vjust = 2),
            axis.ticks        = element_line(color = "black"),
            axis.ticks.length = unit(base_size * 0.25, "points"),
            aspect.ratio      = 0.75)

    return(if(plot) plot.obj else plot.data)
  }
