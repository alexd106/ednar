
#' @title qPCR Calibration Plot
#'
#' @description Create qPCR calibation plots.
#'
#' @param data A \code{\link{tibble}} or \code{\link{data.frame}} containing calibration curve data.
#' @param target A character string containing a unique identifier for the target calibration curve to be plotted.
#' @param lod A vector or \code{tibble} or \code{data.frame} specifying LOD and optionally LOQ values to plot.
#' @param ... Placeholder for further arguments that might be needed by future implementations.
#'
#' @details The \code{data} object contains data for at least one qPCR calibration curve
#'    usually presented as Cq (Ct) value and corresponding copy number from from a series of serial dilutions.
#'    The \code{data.frame} must contain the headers \code{Target}, \code{Cq} and \code{SQ}. The \code{Target}
#'    column must contain unique identifiers for each calibration curve. The \code{Cq} column contains the
#'    Cq (Ct) values and \code{SQ} contains the copy number data. Additional columns will be ignored.
#'
#'    The \code{target} argument takes a single character string specifying the target calibration curve
#'    to plot.
#'
#'    The \code{lod} argument can be supplied as a vector or as a \code{data.frame} or \code{tibble}. If
#'    supplied as a vector then a single LOD value or optionally an LOD and LOQ value can be included. If both LOD
#'    and LOQ values are specified then LOD much be supplied first in the vector (i.e. \code{c(1.5, 4.1)}).
#'    If supplied as a \code{data.frame} or \code{tibble} the object must contain the headers \code{Targets},
#'    and \code{lod} and with an optional \code{loq} header specifying the calibration curve target, the lod and loq
#'    values respectively.
#'
#'    Non-detections in \code{data} should be represented as \code{NA}.
#'
#' @return A \code{\link{ggplot2}} object
#'
#' @export
#'
#' @importFrom scales comma
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' require(tibble)
#' # lod as tibble
#' lod_data <- tibble(Target = "706", lod = 1.5, loq = 4.3)
#' calib_plot(calib_data, target = "706", lod = lod_data)
#'
#' # lod as vector
#' lod_data <- c(1.5, 4.3)
#' calib_plot(calib_data, target = "706", lod = lod_data)
#'
#' # loq missing
#' lod_data <- tibble(Target = "706", lod = 1.5)
#' calib_plot(calib_data, target = "706", lod = lod_data)
#'}
calib_plot <- function(data, target, lod = NULL, ...){

  stopifnot(!missing(target) || !missing(data))
  sub.data <- subset(data, data$Target %in% c(target))

  if(!is.null(lod) && is.vector(lod)){
    if(length(lod) == 1){
      ggplot(sub.data, aes(x = .data$SQ, y = .data$Cq)) +
        geom_point(alpha = 0.65, size = 2) +
        scale_x_continuous(trans = 'log10') +
        geom_smooth(method = "lm", se = FALSE, colour = "black") +
        geom_vline(xintercept = lod[1], color = "black", size = 0.7) +
        ggtitle(paste0("Target: ", target)) +
        xlab(expression(standard~concentrations~(copies~reaction^{-1}))) +
        ylab("Cq - value") +
        theme_calib()
    } else if(length(lod) == 2){
      ggplot(sub.data, aes(x = .data$SQ, y = .data$Cq)) +
      	geom_point(alpha = 0.65, size = 2) +
      	scale_x_continuous(trans = 'log10') +
      	geom_smooth(method = "lm", se = FALSE, colour = "black") +
      	geom_vline(xintercept = lod[1], color = "black", size = 0.7) +
      	geom_vline(xintercept = lod[2], color = "black", linetype = "dashed", size = 0.7) +
        ggtitle(paste0("Target: ", target)) +
      	xlab(expression(standard~concentrations~(copies~reaction^{-1}))) +
      	ylab("Cq - value") +
      	theme_calib()
    } else {
      stop(paste0(lod, " contains more than two values."))
    }
  } else if(!is.null(lod) && is.data.frame(lod)){
      if("loq" %in% colnames(lod)){
        lod.tmp <- as.numeric(lod[lod$Target == target, "lod"])
        loq.tmp <- as.numeric(lod[lod$Target == target, "loq"])
        ggplot(sub.data, aes(x = .data$SQ, y = .data$Cq)) +
          geom_point(alpha = 0.65, size = 2) +
          scale_x_continuous(trans = 'log10') +
          geom_smooth(method = "lm", se = FALSE, colour = "black") +
          geom_vline(xintercept = lod.tmp, color = "black", size = 0.7) +
          geom_vline(xintercept = loq.tmp, color = "black", linetype = "dashed", size = 0.7) +
          ggtitle(paste0("Target: ", target)) +
          xlab(expression(standard~concentrations~(copies~reaction^{-1}))) +
          ylab("Cq - value") +
          theme_calib()
      } else if(!"loq" %in% colnames(lod)) {
        lod.tmp <- as.numeric(lod[lod$Target == target, "lod"])
        ggplot(sub.data, aes(x = .data$SQ, y = .data$Cq)) +
          geom_point(alpha = 0.65, size = 2) +
          scale_x_continuous(trans = 'log10') +
          geom_smooth(method = "lm", se = FALSE, colour = "black") +
          geom_vline(xintercept = lod.tmp, color = "black", size = 0.7) +
          ggtitle(paste0("Target: ", target)) +
          xlab(expression(standard~concentrations~(copies~reaction^{-1}))) +
          ylab("Cq - value") +
          theme_calib()
      }
  } else {
    ggplot(sub.data, aes(x = .data$SQ, y = .data$Cq)) +
      geom_point(alpha = 0.65, size = 2) +
      scale_x_continuous(trans = 'log10') +
      geom_smooth(method = "lm", se = FALSE, colour = "black") +
      xlab(expression(standard~concentrations~(copies~reaction^{-1}))) +
      ylab("Cq - value") +
      ggtitle(paste0("Target: ", target)) +
      theme_calib()
  }
}
