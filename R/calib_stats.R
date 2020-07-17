#' @title Fit and Summarise qPCR Calibration Curves
#'
#' @description Fit and summarise qPCR calibration curves.
#'
#' @param data A \code{\link{tibble}} or \code{\link{data.frame}} containing data for
#'     one or more calibration curves.
#' @param type The level of detail to return. \code{"effects"} returns slope and intercept estimates, \code{"r2"}
#'     returns the \code{R^2} of the model fits.
#' @param ... Placeholder for further arguments that might be needed by future implementations.
#'
#' @details The calibration curves are fitted using a linear model with the \code{\link{lm}} function. \code{SQ}
#'     values are \code{\link{log10}} transformed prior to model fitting.
#'
#'     The \code{data} object contains data for at least one qPCR calibration curve
#'     usually presented as Cq (Ct) value and corresponding copy number from from a series of serial dilutions.
#'     The \code{data} object must contain the headers \code{Target}, \code{Cq} and \code{SQ}. The \code{Target}
#'     column must contain unique identifiers for each calibration curve. The \code{Cq} column contains the
#'     Cq (Ct) values and \code{SQ} contains the copy number data. Additional columns will be ignored.
#'
#'     The \code{"type"} argument can only take the character strings \code{"effects"} or \code{"r2"}.
#'     Specifying \code{"effects"} will return the slope and intercept estimates, standard error for each estimate,
#'     the test statistic and p.value for each calibration curve. Specifying \code{"r2"} returns the \code{R^2}
#'     for each model.
#'
#'     Non-detections in \code{data} should be represented as \code{NA}.
#'
#' @return A \code{tibble} containing summary statistics of model fit for each calibration curve.
#'
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom broom tidy glance
#' @importFrom dplyr group_by mutate
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom tidyr nest unnest
#'
#' @examples
#' \dontrun{
#' sum_stats <- calib_stats(calib_data, type = "effects")
#' sum_stats <- calib_stats(calib_data, type = "r2")
#' }
calib_stats <- function(data, type = "effects", ...){
	stopifnot(!missing(data) || class(data == "data.frame") || class(data == "tibble"))

	data <- tibble::tibble(data)
	data$Target <- factor(data$Target)

	out <- data %>%
		tidyr::nest(data = -(.data$Target)) %>%
		dplyr::mutate(
			fit = purrr::map(data, ~ lm(Cq ~ log10(SQ), data = .x)),
			tidied = purrr::map(.data$fit, broom::tidy),
			glanced = purrr::map(.data$fit, broom::glance)
		)

	if(type == "effects"){
		out %>%
			tidyr::unnest(.data$tidied) %>%
			dplyr::select(.data$Target, .data$term, .data$estimate, .data$std.error, .data$statistic, .data$p.value)
	} else if(type == "r2"){
		out %>%
			tidyr::unnest(.data$glanced) %>%
			dplyr::select(.data$Target, .data$r.squared)
	} else {
		stop(paste0("type = must contain one of either 'effects' or 'r2'."))
	}
}
