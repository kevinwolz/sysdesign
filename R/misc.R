yard_to_feet_frac <- function(y) {

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
