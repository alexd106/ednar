
#' ggplot2 Theme For Calibration Curve Plots
#'
#' @param base_size base size
#' @param base_family base family
#' @param base_line_size base_line_size
#' @param base_rect_size base_rect_size
#'
#' @return A \code{\link{ggplot2}} object
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' theme_calib
#' }
theme_calib <- function (base_size = 11, base_family = "", base_line_size = base_size/22,
												 base_rect_size = base_size/22)
{
	half_line <- base_size/2
	theme_grey(base_size = base_size, base_family = base_family,
	  base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
	theme(panel.background = element_rect(fill = "white", colour = NA),
	  panel.border = element_rect(fill = NA, colour = "grey70", size = rel(1)),
	  panel.grid = element_line(colour = "grey87"),
		panel.grid.major = element_line(size = rel(0.5)),
		panel.grid.minor = element_line(size = rel(0.25)),
		axis.ticks = element_line(colour = "grey70", size = rel(0.5)),
		axis.title = element_text(size = 13),
		axis.text.x = element_text(size = 11),
		axis.text.y = element_text(size = 11),
		legend.key = element_rect(fill = "white", colour = NA),
		strip.background = element_rect(fill = "grey70", colour = NA),
		strip.text = element_text(colour = "white", size = rel(0.8),
		  margin = margin(0.8 * half_line, 0.8 * half_line, 0.8 * half_line, 0.8 * half_line)),
		plot.margin = unit(c(5.5, 15, 5.5, 5.5), "points"),
	complete = TRUE)
}

