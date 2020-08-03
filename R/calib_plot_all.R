#' @title Plot Multiple qPCR Calibration Plots
#'
#' @description Create and plot multiple qPCR calibration plots.
#'
#' @param calib_df A \code{\link{tibble}} or \code{\link{data.frame}} containing calibration curve data.
#' @param ... Further arguments to be passed to the \code{\link{calib_plot}} function.
#'
#' @details The \code{data} object contains data for at least one qPCR calibration curve
#'    usually presented as Cq (Ct) value and corresponding copy number from from a series of serial dilutions.
#'    The \code{data.frame} must contain the headers \code{Target}, \code{Cq} and \code{SQ}. The \code{Target}
#'    column must contain unique identifiers for each calibration curve. The \code{Cq} column contains the
#'    Cq (Ct) values and \code{SQ} contains the copy number data. Additional columns will be ignored.
#'
#'    Non-detections in \code{data} should be represented as \code{NA}.
#'
#'    LOD and LOQ values must be supplied as a \code{\link{data.frame}} or \code{\link{tibble}} object.
#'
#' @return A \code{\link{list}} of \code{\link{ggplot2}} objects.
#'
#' @export
#'
#' @importFrom patchwork wrap_plots
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' require(tibble)
#'
#' # basic plots
#' calib_plot_all(calib_data)
#'
#' # lod and loq as tibble
#' lod_data <- tibble(Target = c("706", "10720201", "90720202" ),
#'               lod = c(1.5, 4, 10), loq = c(4.3, 10, 20))
#' calib_plot_all(calib_data, lod = lod_data)
#'
#' # Robust model fit
#' calib_plot_all(calib_data, lod = lod_data, robust = TRUE)
#'}
#'

calib_plot_all <- function(calib_df, ...){

	stopifnot(!missing(calib_df))
	calib_df_name <- deparse(substitute(calib_df))

	# check columns names supplied
	if(!any(colnames(calib_df) == "Target")) {
		stop(paste0(calib_df, " does not contain Targets column."))
	}
	if(!any(colnames(calib_df) == "Cq")) {
		stop(paste0(calib_df, " does not contain Cq column."))
	}
	if(!any(colnames(calib_df) == "SQ")) {
		stop(paste0(calib_df, " does not contain SQ column."))
	}
	if(length(unique(calib_df$Target)) > 8) {
		stop(paste0("More than 8 Targets detected in ", calib_df_name, ". Consider subsetting your data before plotting."))
	}

	#list object to store plots
	plot_list <- vector(mode = "list", length = length(unique(calib_df$Target)))
	names(plot_list) <- as.character(unique(calib_df$Target))

	# plots for each target
	for(i in seq_along(unique(calib_df$Target))){
	  plot_tmp <- calib_plot(calib_df, target = unique(calib_df$Target)[i], ...)
	  plot_list[[i]] <- plot_tmp
	}
	patchwork::wrap_plots(plot_list, ncol = 2)
}



