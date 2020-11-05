#' @title LOD Plot
#'
#' @description Creates a LOD plot for a specified 'Target'. The LOD plot shows the relative detection rates for
#'     each standard as well as the LOD model curve. This code has been modified and functionalised from the original
#'     by Merkes et al. (see references).
#'
#' @param lod_obj A \code{lod} object containing the output from the  \code{\link{calib_lod}} function.
#' @param target A character string containing a unique identifier for the target calibration curve to be plotted.
#' @param legend A logical indicating whether a lengen should be included on the plot, Default: \code{TRUE}.
#' @param ... Placeholder for further arguments that might be needed by future implementations.
#'
#' @details The LOD plot shows the relative detection rates for each standard as well as the LOD model curve.
#'     The 95% LOD is identified, and effective LODs with confidence intervals for analysing samples with
#'     multiple replicates are displayed. The logarithmic function that was used to determine the LODs is
#'     shown in the plot subtitle along with the p‐value for a lack of fit test on the model.
#'
#' @return A LOD plot in the current plotting device and a \code{\link{data.frame}} containing the
#'     effective LOD estimates, standard errors, upper and lower confidence intervals for each replicate.
#'
#' @references Merkes CM, Klymus KE, Allison MJ, Goldberg C, Helbing CC, Hunter ME, Jackson CA, Lance RF,
#'     Mangan AM, Monroe EM, Piaggio AJ, Stokdyk JP, Wilson CC, Richter C. (2019) Generic qPCR Limit of
#'     Detection (LOD) / Limit of Quantification (LOQ) calculator. R Script.
#'     Available at: https://github.com/cmerkes/qPCR_LOD_Calc. DOI: https://doi.org/10.5066/P9GT00GB.
#'
#'     Klymus, Katy E., Christopher M. Merkes, Michael J. Allison, Caren S. Goldberg, Caren C. Helbing,
#'     Margaret E. Hunter, Craig A. Jackson, et al. ‘Reporting the Limits of Detection and Quantification
#'     for Environmental DNA Assays’. Environmental DNA 2, no. 3 (July 2020): 271–82.
#'     https://doi.org/10.1002/edn3.29.
#'
#' @export
#'
#' @importFrom drc ED modelFit
#' @importFrom grDevices rgb
#' @importFrom graphics legend lines mtext plot points
#'
#' @examples
#' \dontrun{
#' qpcr_lod <- calib_lod(data = calib_data, threshold = 0.35,
#'               lod.fit = "best", loq.fit = "best")
#' lod_CI <- calib_plot_lod(lod_obj = qpcr_lod, target = "10720201")
#' lod_CI
#'}
#'
calib_plot_lod <- function(lod_obj, target, legend = TRUE, ...) {

	stopifnot(!missing(lod_obj))
	stopifnot(!missing(target))

	calib_obj_name <- deparse(substitute(lod_obj))
	target_name <- deparse(substitute(target))

	## make sure fit is a lod object
	if(!inherits(lod_obj, "lod")){
		stop(paste(calib_obj_name, "is not an lod object"))
	}

  if(class(target) != "character"){
  	stop(paste0(target_name, " must be class 'character'."))
  }

	data_obj <- lod_obj

	if(sum(grepl(paste0("\\b", target, "\\b"), names(data_obj$LODlist))) == 0) {
		stop(paste0("no target name: ", target_name, " found in ", calib_obj_name))
	}

  mod_obj <- data_obj$LODlist[grepl(paste0("\\b", target, "\\b"), names(data_obj$LODlist))]

  if (!is.na(names(mod_obj))) {
    DAT4 <- rbind(
      drc::ED(mod_obj[[1]], 0.95, interval = "delta", type = "absolute"),
      drc::ED(mod_obj[[1]], 1 - sqrt(0.05), interval = "delta", type = "absolute"),
      drc::ED(mod_obj[[1]], 1 - 0.05^(1 / 3), interval = "delta", type = "absolute"),
      drc::ED(mod_obj[[1]], 1 - 0.05^0.25, interval = "delta", type = "absolute"),
      drc::ED(mod_obj[[1]], 1 - 0.05^0.2, interval = "delta", type = "absolute"),
      drc::ED(mod_obj[[1]], 1 - 0.05^0.125, interval = "delta", type = "absolute")
    )
    if (substr(mod_obj[[1]]$fct$name, 1, 4) == "LL2") {
      DAT4 <- exp(DAT4)
    }
    DAT4 <- data.frame(DAT4,
      LoD = c(
        "1rep.LOD", "2rep.LOD", "3rep.LOD", "4rep.LOD",
        "5rep.LOD", "8rep.LOD"
      ),
      Assay = rep(target, nrow(DAT4))
    )
    DAT4$Assay <- as.character(DAT4$Assay)
    LOD.CI <- DAT4
    if (sum(!is.na(DAT4[, 3]) & DAT4[, 3] <= 0) > 0) { # Unable to plot negative lower limits, converting any lower limit values to 0.0001
      DAT4[!is.na(DAT4[, 3]) & DAT4[, 3] <= 0, 3] <- 0.0001
    }
    plot(mod_obj[[1]],
      main = paste0("LOD Plot for: ", target),
      ylab = "Detection Probability", xlab = "Standard concentrations (Copies / Reaction)",
      xlim = c(min(DAT4[, 1:4], na.rm = TRUE), max(mod_obj[[1]]$origData$SQ, na.rm = TRUE))
    )
    LODS <- sum(!is.na(DAT4[, 1]))
    COLS <- c(
      rgb(0.8, 0.47, 0.65), rgb(0, 0.45, 0.7), rgb(0.94, 0.89, 0.26),
      rgb(0.84, 0.37, 0), rgb(0, 0.62, 0.45), rgb(0.90, 0.62, 0)
    )
    PNTS <- c(15, 16, 17, 18, 25, 3)
    YS <- c(0.95, 1 - sqrt(0.05), 1 - 0.05^(1 / 3), 1 - 0.05^0.25, 1 - 0.05^0.2, 1 - 0.05^0.125)
    LODS2 <- c(
      "Limit of Detection", "2 Replicates LoD", "3 Replicates LoD",
      "4 Replicates LoD", "5 Replicates LoD", "8 Replicates LoD"
    )
    if (LODS < 6) {
      LODS2[(LODS + 1):6] <- gsub("licates LoD", "s: Insufficient Data", LODS2[(LODS + 1):6])
    }
    points(x = DAT4[1:LODS, 1], y = YS[1:LODS], pch = PNTS, col = COLS, cex = 1.2)
    for (j in 1:LODS) {
      lines(x = DAT4[j, 3:4], y = rep(YS[j], 2), col = COLS[j], lwd = 2)
      lines(x = rep(DAT4[j, 3], 2), y = c(YS[j] - 0.02, YS[j] + 0.02), lwd = 2, col = COLS[j])
      lines(x = rep(DAT4[j, 4], 2), y = c(YS[j] - 0.02, YS[j] + 0.02), lwd = 2, col = COLS[j])
    }
    if(legend == TRUE){
    legend("bottomright", legend = LODS2, pch = PNTS, col = COLS, text.col = COLS)
    }
    Pval <- drc::modelFit(mod_obj[[1]])[[5]][2]
    mtext(paste0("FCT used: ", mod_obj[[1]]$fct$name, "    Lack of fit test: p = ", round(Pval, digits = 4)), side = 3)
  }
  if (is.na(names(mod_obj))) {
    plot(data_obj$standardsSum$Rate[data_obj$standardsSum$Target == target] ~ log10(data_obj$standardsSum$Standards[data_obj$standardsSum$Target == target]),
      ylim = c(0, 1), ylab = "Detection Probability",
      xlab = expression("Log of standard concentrations (Log"[10] * "Copies / Reaction)"),
      main = paste0("LoD for: ", target, " unsolvable")
    )
  }
  invisible(LOD.CI)
}
