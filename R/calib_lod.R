#' @title Calculate Model Based LOD and LOQ for qPCR eDNA Assays
#'
#' @description Calculate threshold and model based LOD and LOQ for qPCR eDNA assays based on one or
#'     more standard curves. This code has been modified and functionalised from the original
#'     by Merkes et al. (see references).
#'
#' @param data A \code{\link{data.frame}} or \code{\link{tibble}} object containing calibration curve data.
#' @param threshold Threshold value for LOQ calculation using coefficient of variation (CV) of Cq (Ct)
#'     values. Default value = 0.35.
#' @param lod.fit Character string to specify method to estimate LOD. Default: \code{"best"}. See details for other
#'     possible values.
#' @param loq.fit Character string to specify method to estimate LOQ based on CV curves. Default: \code{'best'}.
#'     Other option include \code{'decay'}, \code{'linear'} and \code{'Pn'} (see details).
#' @param robust A logical value indicating whether the fitted model should exclude standards with less than 50\%
#'     detections. Default = \code{FALSE}.
#' @param ... Placeholder for further arguments that might be needed by future implementations.
#'
#' @details The \code{data} object contains data for at least one qPCR calibration curve
#'    usually presented as Cq (Ct) value and corresponding copy number from from a series of serial dilutions.
#'    The \code{data.frame} must contain the headers \code{Target}, \code{Cq} and \code{SQ}. The \code{Target}
#'    column must contain unique identifiers for each calibration curve. The \code{Cq} column contains the
#'    Cq (Ct) values and \code{SQ} contains the copy number data. Additional columns will be ignored. Non-detections
#'    in \code{data} should be represented as \code{NA}.
#'
#'    The \code{lod.fit} argument is used to select the method to estimate LOD. Selecting "best" will automatically
#'    fit all models available from the package \code{drc} and perform model selection based on log
#'    likelihood values, Akaike's information criterion and residual variance using the \code{\link[drc]{mselect}}
#'    function. Information on possible models can be displayed with the \code{\link[drc]{getMeanFunctions}}
#'    function. For example, to fit a 4 parameter Weibull type II model use \code{lod.fit = "W2.4"}.
#'
#'    The \code{loq.fit} argument is used to select the method to estimate LOQ. Selecting "best" will automatically
#'    select the model with lowest residual standard error. 'decay' uses the exponential decay model, 'linear'
#'    uses a linear model and 'Pn' fits an nth-order polynomial model where n is numerical. For example,
#'    \code{loq.fit = "P2"} will fit a 2nd-order polynomial model and \code{loq.fit = "P3"} will use a 3rd-order.
#'    Selecting "best" will test polynomial models up to 6th-order.
#'
#' @return A \code{\link{list}} object containing the following objects:
#'
#' \code{dataSum}: A summary \code{data.frame} of the original qPCR calibration data including copy number estimate (\code{Copy.Estimate}).
#'
#' \code{standardsSum}: A summary \code{data.frame} of the standards containing the following variables: number of replicates (\code{Rep}), number of detections (\code{Detects})
#'    and detection rate (\code{Rate}) for each standard (\code{Standards}) and Target (\code{Target}) combination.
#'    The mean (\code{Cq.mean}), standard deviation (\code{Cq.sd}) and coefficient of variation (\code{Cq.CV}) of
#'    Cq values and also the coefficient of variation of copy number (\code{Copy.CV}) is given.
#'
#' \code{assaySum}: A summary \code{data.frame} of the qPCR assay containing the following variables:
#' \describe{
#' \item{R.squared:}{The R-squared value of linear regression of all standards Cq-values vs log10 of the starting quantities}
#' \item{Slope:}{The slope of the linear regression}
#' \item{Intercept:}{The intercept of the linear regression}
#' \item{Low.95:}{The lowest standard with at least 95 percent positive detection}
#' \item{LOD:}{The 95 percent limit of detection as determined by dose-response modelling}
#' \item{LOQ:}{The limit of quantification as determined by decay modelling, using the user-selected CV threshold}
#' \item{rep2.LOD:}{The effective limit of detection if analysing in 2 replicates}
#' \item{rep3.LOD:}{The effective limit of detection if analysing in 3 replicates}
#' \item{rep4.LOD:}{The effective limit of detection if analysing in 4 replicates}
#' \item{rep5.LOD:}{The effective limit of detection if analysing in 5 replicates}
#' \item{rep8.LOD:}{The effective limit of detection if analysing in 8 replicates}
#' \item{LOD.model}{The optimal model when using \code{lod.fit = "best"} model selection, or the user specified model}
#' }
#'
#' \code{LOQlist}: A \code{list} object containing LOQ model components for each Target.
#'
#' \code{LODlist}: A \code{list} object containing LOD model components for each Target.
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
#' @examples
#' \dontrun{
#' #
#' qpcr_lod <- calib_lod(data = calib_data, threshold = 0.35,
#'               lod.fit = "best", loq.fit = "best")
#' qpcr_lod
#'
#' # decay modelling estimate for LOQ
#' qpcr_lod_decay <- calib_lod(data = calib_data, threshold = 0.35,
#'                     lod.fit = "best", loq.fit = "decay")
#' qpcr_lod_decay
#'
#' # 4 parameter Weibull type II model fitted to estimate LOD
#'
#' qpcr_lod_W2 <- calib_lod(data = calib_data, threshold = 0.35,
#'                  lod.fit = "W2.4", loq.fit = "best")
#' qpcr_lod_W2
#'  }
#'
#'
#' @export
#'
#' @importFrom stats lm nls predict coef median sd
#' @importFrom drc drm mselect getMeanFunctions ED AR.2 AR.3 LL.2 LL.3 LL.3u LL.4 LL.5 MM.2 MM.3 W1.2 W1.3 W1.4 W2.2 W2.3 W2.4
#'
calib_lod <- function(data, threshold = 0.35, lod.fit = "best", loq.fit = "best", robust = FALSE, ...) {

	stopifnot(!missing(data) || class(data == "data.frame"))
	stopifnot(class(threshold ) == "numeric")
	if (threshold > 1 | threshold < 0){
		stop("'threshold' must be between 0 and 1")
	}
	stopifnot(class(loq.fit) == "character")
	stopifnot(class(lod.fit) == "character")
	DAT <- data

	# check columns names supplied
	if(!any(colnames(DAT) == "Target")) {
		stop(paste0(data, " does not contain Targets column."))
	}
	if(!any(colnames(DAT) == "Cq")) {
		stop(paste0(data, " does not contain Cq column."))
	}
	if(!any(colnames(DAT) == "SQ")) {
		stop(paste0(data, " does not contain SQ column."))
	}

  LOQ.Threshold <- threshold
  if(lod.fit == "best") {
  	LOD.FCT <- lod.fit
  } else {
  	LOD.FCT <- drc::getMeanFunctions(fname = lod.fit)[[1]]
  }
  LOQ.FCT <- loq.fit

  ## Ensure data is in the proper format
  DAT$Target <- as.factor(DAT$Target)

  ## Check for wild outliers that the user should go back and review:
  Targets <- unique(DAT$Target)
  ## Get matchups of all standards and markers used:
  for(i in seq_along(Targets)) {
  	if(i == 1) {
  		Standards <- unique(DAT$SQ[DAT$Target == Targets[i] & !is.na(DAT$SQ)])
  		Target <- rep(as.character(Targets[i]), length(Standards))
  	}
  	else {
  		Standards <- c(Standards, unique(DAT$SQ[DAT$Target == Targets[i] & !is.na(DAT$SQ)]))
  		Target <- c(Target, rep(as.character(Targets[i]),
  													 length(unique(DAT$SQ[DAT$Target == Targets[i] & !is.na(DAT$SQ)]))))
  	}
  }
  OUTS <- data.frame(Target = Target, Standard = Standards,Outliers = NA)

  ## Identify any wells where the Cq value is more than 10%  away from the median for
  ##   that standard.
 for (i in 1:nrow(OUTS)) {
  MED <- median(DAT$Cq[DAT$SQ == OUTS$Standard[i] & DAT$Target == OUTS$Target[i]], na.rm = TRUE)
  A <- which(DAT$SQ == OUTS$Standard[i] & DAT$Target == OUTS$Target[i] & DAT$Cq < 0.9 * MED & !is.na(DAT$Cq))
  B <- which(DAT$SQ == OUTS$Standard[i] & DAT$Target == OUTS$Target[i] & DAT$Cq > 1.1 * MED & !is.na(DAT$Cq))
  if (length(c(A, B)) > 0) {
    OUTS$Outliers[i] <- paste(c(A, B), collapse = ",")
    }
  }
  ## If any outliers are detected
  if (sum(!is.na(OUTS$Outliers)) > 0) {
  OUT.ROW <- paste(OUTS$Outliers[!is.na(OUTS$Outliers)], collapse = ",")
  OUT.ROW2 <- unlist(strsplit(OUT.ROW, split = ","))
  message("warning: Potential outliers detected: ")
  print(DAT[OUT.ROW2,])
  }

  curve.list <- rep(NA, length(Targets))
  DAT$Copy.Estimate <- rep(NA, nrow(DAT))
  DAT$Mod <- rep(0, nrow(DAT))
    for (i in seq_along(Targets)) {
    STDS <- data.frame(S = unique(DAT$SQ[DAT$Target == Targets[i]]), R = NA)
    ## Calculate detection rates for each standard:
    for (j in 1:nrow(STDS)) {
      STDS$R[j] <- sum(!is.na(DAT$Cq) & DAT$SQ == STDS$S[j] & DAT$Target == Targets[i], na.rm = TRUE) / sum(DAT$SQ == STDS$S[j] & DAT$Target == Targets[i], na.rm = TRUE)
    }
    if(robust == TRUE){
    	if (sum(STDS$R >= 0.5, na.rm = TRUE) > 2) {
    		# Only use standards with 50% or greater detection rates for linear regression
    		STDS2 <- STDS$S[STDS$R >= 0.5 & !is.na(STDS$R) & !is.na(STDS$S)]
    	}
    	# If there isn't at least 3 standards with 50% or greater detection, use the top 3
    	if (sum(STDS$R >= 0.5, na.rm = TRUE) < 3) {
    		STDS2 <- STDS$S[order(STDS$R, decreasing = TRUE)][1:3]
    	}
    } else {
    	STDS2 <- STDS$S[!is.na(STDS$R) & !is.na(STDS$S)]
    }

    # fit lm's
    for (j in seq_along(STDS2)) {
      ind.cq <- DAT$Cq[DAT$Target == Targets[i] & DAT$SQ == STDS2[j]]

      DAT$Mod[DAT$Target == Targets[i] & DAT$SQ == STDS2[j] & !is.na(DAT$Cq) & !is.na(DAT$SQ)] <- 1
    }
    assign(paste0("curve", i), stats::lm(Cq ~ log10(SQ), data = DAT[DAT$Target == Targets[i] & DAT$Mod == 1, ]))
    curve.list[i] <- paste0("curve", i)
    Intercept <- coef(get(curve.list[i]))[1]
    Slope <- coef(get(curve.list[i]))[2]
    DAT$Copy.Estimate[DAT$Target == Targets[i]] <- 10^((DAT$Cq[DAT$Target == Targets[i]] - Intercept) / Slope)
  }

  ## Summarise the data:
  DAT2 <- data.frame(
  Standards = Standards, Target = Target, Reps = NA, Detects = NA, Cq.mean = NA,
  Cq.sd = NA, Copy.CV = NA, Cq.CV = NA)

  ## Fill in replicate counts, positive detect counts, mean Cq values, standard
  ##   deviations of Cq values, and coefficient of variation of copy estimates for
  ##   each standard and marker combination:
  for(i in 1:nrow(DAT2)) {
  	DAT2$Reps[i] <- sum(DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i],na.rm=TRUE)
  	DAT2$Detects[i] <- sum(!is.na(DAT$Cq)&DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i],na.rm=TRUE)
  	DAT2$Cq.mean[i] <- mean(DAT$Cq[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]],na.rm=TRUE)
  	DAT2$Cq.sd[i] <- sd(DAT$Cq[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]],na.rm=TRUE)
  	DAT2$Copy.CV[i] <- sd(DAT$Copy.Estimate[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]],na.rm=TRUE)/mean(DAT$Copy.Estimate[DAT$SQ==DAT2$Standards[i]&DAT$Target==DAT2$Target[i]],na.rm=TRUE)
  	DAT2$Cq.CV[i] <- sqrt(2^(DAT2$Cq.sd[i]^2*log(2))-1) # enhance to use specific E for each plate
  }

  ## Calculate positive detection rate for each standard and marker combination:
  DAT2$Rate <- DAT2$Detects/DAT2$Reps
  ## Determine LoD and LoQ by modeling, and summarize each assay:
  DAT$Detect <- as.numeric(!is.na(DAT$Cq))
  LOD.list2 <- vector(mode = "list", length = length(Targets))
  LOD.list3 <- vector(mode = "list", length = length(Targets))
  LOQ.list <- vector(mode = "list", length = length(Targets))
  LOD.out <- vector(mode = "list", length = length(Targets))
  LOQ.out <- vector(mode = "list", length = length(Targets))
  DAT3 <- data.frame(Assay = Targets, R.squared = NA, Slope = NA, Intercept = NA, Low.95 = NA,
                     LOD = NA, LOQ = NA, rep2.LOD = NA, rep3.LOD = NA, rep4.LOD = NA, rep5.LOD = NA,
  									 rep8.LOD = NA, Efficiency = NA, LOD.model = NA)
  LOD.FCTS <- list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), W1.2(), W1.3(), W1.4(), W2.2(), W2.3(),
                   W2.4(), AR.2(), AR.3(), MM.2(), MM.3())

  for(i in seq_along(Targets)) {
  	## Check input suitability for dose-response modelling:
  	if(sum(DAT2$Rate[DAT2$Target == Targets[i]] != 1 & DAT2$Rate[DAT2$Target == Targets[i]] != 0) == 0) {
  		ToWrite <- paste0("WARNING: For ", Targets[i],", all standards detected fully or failed fully.  Therefore, the LoD model will not converge.")
  		print(ToWrite)
  	}
  	if(sum(DAT2$Rate[DAT2$Target == Targets[i]] != 1 & DAT2$Rate[DAT2$Target == Targets[i]] != 0) == 1) {
  		ToWrite <- paste0("WARNING: For ", Targets[i],", only 1 standard detected in the informative range (not 0% and not 100%).  Therefore, the LoD model results will be less reliable.")
  		print(ToWrite)
  	}

  	## Define LOQ model using lowest residual standard error selection:
  	if(LOQ.FCT == "best") {
  		## Remove previous marker LOQ models from environment if they exist:
  		suppressWarnings(rm(LOQ1, LOQ2, LOQ3, LOQ4, LOQ5, LOQ6, LOQ7))

  		tryCatch({ #skip if model cannot be determined.
  			LOQ1 <- stats::nls(Cq.CV ~ SSasymp(log10(Standards), Asym, R0, lrc),
                 data = DAT2[DAT2$Target == Targets[i], ])
  		}, error = function(e) {
  			e
  			cat("warning: decay LOQ model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  		})
  		tryCatch({ #skip if model cannot be determined.
  			LOQ2 <- stats::lm(Cq.CV ~ log10(Standards), data = DAT2[DAT2$Target == Targets[i], ])
  		}, error = function(e) {
  			e
  			cat("warning: linear LOQ model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  		})
  		tryCatch({ #skip if model cannot be determined.
  			LOQ3 <- stats::lm(Cq.CV ~ poly(log10(Standards), 2), data = DAT2[DAT2$Target == Targets[i], ])
  		}, error = function(e) {
  			e
  			cat("warning: 2nd polynomial LOQ model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  		})
  		tryCatch({ #skip if model cannot be determined.
  			LOQ4 <- stats::lm(Cq.CV ~ poly(log10(Standards), 3), data = DAT2[DAT2$Target == Targets[i], ])
  		}, error = function(e) {
  			e
  			cat("warning: 3rd polynomial LOQ model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  		})
  		tryCatch({ #skip if model cannot be determined.
  			LOQ5 <- stats::lm(Cq.CV ~ poly(log10(Standards), 4), data = DAT2[DAT2$Target == Targets[i], ])
  		}, error = function(e) {
  			e
  			cat("warning: 4th polynomial LOQ model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  		})
  		tryCatch({ #skip if model cannot be determined.
  			LOQ6 <- stats::lm(Cq.CV ~ poly(log10(Standards), 5), data = DAT2[DAT2$Target == Targets[i], ])
  		}, error = function(e) {
  			e
  			cat("warning: 5th polynomial LOQ model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  		})
  		tryCatch({ #skip if model cannot be determined.
  			LOQ7 <- stats::lm(Cq.CV ~ poly(log10(Standards), 6), data = DAT2[DAT2$Target == Targets[i], ])
  		}, error = function(e) {
  			e
  			cat("warning: 6th polynomial LOQ model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  		})

  		## Determine which models were able to be determined:
  		A <- sapply(c("LOQ1", "LOQ2", "LOQ3", "LOQ4", "LOQ5", "LOQ6", "LOQ7"), function (x) { exists(x) })
  		B <- names(A)[A == TRUE]
  		## If at least 1 LOQ model was determined, select the one with the lowest
  		##   residual standard error:
  		if(length(B) > 0) {
  			LOQ.res <- rep(NA, length(B))
  			for(j in seq_along(B)) {
  				LOQ.res[j] <- summary(get(B[j]))$sigma
  			}
  			C <- which(LOQ.res == min(LOQ.res, na.rm = TRUE))
  			assign(paste0("LOQ.mod", i), get(B[C]))
  			LOQ.list[[i]] <- paste0("LOQ.mod", i)
  		}
  	}
  	## Define LOQ model by exponential decay modelling:
  	if(LOQ.FCT == "decay") {
  		tryCatch({ #skip if model cannot be determined.
  			assign(paste0("LOQ.mod", i), stats::nls(Cq.CV ~ SSasymp(log10(Standards), Asym, R0, lrc),
                                      data = DAT2[DAT2$Target == Targets[i], ]))
  			LOQ.list[[i]] <- paste0("LOQ.mod", i)
  		}, error = function(e) {
  			e
  			cat("warning: decay LOQ model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  		})
  	}
  	## Define LOQ model by linear modelling:
  	if(LOQ.FCT == "linear") {
  		tryCatch({ #skip if model cannot be determined.
  			assign(paste0("LOQ.mod", i), stats::lm(Cq.CV ~ log10(Standards),
  																		data = DAT2[DAT2$Target == Targets[i], ]))
  			LOQ.list[[i]] <- paste0("LOQ.mod", i)
  		}, error = function(e) {
  			e
  			cat("warning: linear LOQ model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  		})
  	}
  	## Define LOQ model by polynomial modelling:
  	if(substr(LOQ.FCT, 1, 1) == "P") {
  		Z <- as.numeric(substr(LOQ.FCT, 2, nchar(LOQ.FCT)))
  		tryCatch({ #skip if model cannot be determined.
  			assign(paste0("LOQ.mod", i),stats::lm(Cq.CV ~ poly(log10(Standards), Z),
  																		data = DAT2[DAT2$Target == Targets[i], ]))
  			LOQ.list[[i]] <- paste0("LOQ.mod", i)
  		}, error = function(e) {
  			e
  			cat("warning: ", Z, "-order polynomial LOQ model cannot be defined for ",
  					as.character(Targets[i]), "\n", sep = "")
  		})
  	}
  	## Signal undetermined model with NA:
  	if(length(LOQ.list) < i) {
  		LOQ.list[[i]] <- NA
  	}
  	## Define the logarithmic model for LOD using user-selected function:
  	if(is.list(LOD.FCT) == TRUE) {
  		tryCatch({ #skip if model cannot be determined.
  			assign(paste0("LOD.mod2", i), drc::drm(Detect ~ SQ, data = DAT[DAT$Target == Targets[i], ], fct = LOD.FCT))
  			LOD.list2[[i]] <- paste0("LOD.mod2", i)
  			LOD.list3[[i]] <- LOD.FCT$name
  		}, error = function(e) {
  			e
  			cat("warning: LOD model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  		})
  	}
  	## Define the logarithmic model with function automatically selected:
  		if(LOD.FCT[1] == "best") {
  			tryCatch({ #skip if model cannot be determined.
  				## Pull out data for specific assay:
  				TEMP.DAT <- DAT[DAT$Target == Targets[i], ]
  				## Define a model to start with:

  				LOD.mod <- drc::drm(Detect ~ SQ, data = TEMP.DAT, fct = W2.4())
  				## Test all available models and select the best one:
  				LOD.FCT2 <- suppressWarnings(row.names(drc::mselect(LOD.mod, LOD.FCTS))[1])
  				LOD.FCT3 <- drc::getMeanFunctions(fname = LOD.FCT2)
  				assign(paste0("LOD.mod2", i), drc::drm(Detect ~ SQ, data = DAT[DAT$Target == Targets[i], ], fct = LOD.FCT3[[1]]))
  				LOD.list2[[i]] <- paste0("LOD.mod2", i)
  				LOD.list3[[i]] <- LOD.FCT2
  			}, error = function(e) {
  				e
  				cat("warning: LOD model cannot be defined for ", as.character(Targets[i]), "\n", sep = "")
  			})
  		}

  	## Signal undetermined model with NA:
  	if(length(LOD.list2) < i) {
  		LOD.list2[[i]] <- NA
  		LOD.list3[[i]] <- NA
  	}

  	# LOD.list2 output
  	  LOD.out[[i]] <- get(LOD.list2[[i]])
  	  names(LOD.out)[i] <- LOD.list2[[i]]

  	# LOQ.list output
  	  LOQ.out[[i]] <- get(LOQ.list[[i]])
  		names(LOQ.out)[i] <- LOQ.list[[i]]

  	## Populate summary data:
  	DAT3$R.squared[i] <- summary(get(curve.list[i]))$r.squared
  	DAT3$Slope[i] <- coef(get(curve.list[i]))[2]
  	DAT3$Intercept[i] <- coef(get(curve.list[i]))[1]
  	DAT3$Efficiency[i] <- 10^((-1/as.numeric(coef(get(curve.list[i]))[2])))-1
  	DAT3$Low.95[i] <- min(DAT2$Standards[DAT2$Rate >= 0.95 & DAT2$Target == Targets[i]])
  	## Only get LOD values if the LOD model is defined:
  	if(!is.na(LOD.list2[i])) {
  		DAT3$LOD.model[i] <- LOD.list3[[i]]
  		DAT3$LOD[i] <- suppressWarnings(drc::ED(get(LOD.list2[[i]]), 0.95, type = "absolute", display = FALSE)[1])
  		DAT3$rep2.LOD[i] <- suppressWarnings(drc::ED(get(LOD.list2[[i]]), 1 - sqrt(0.05), type = "absolute", display = FALSE)[1])
  		DAT3$rep3.LOD[i] <- suppressWarnings(drc::ED(get(LOD.list2[[i]]), 1 - 0.05^(1/3), type = "absolute", display = FALSE)[1])
  		DAT3$rep4.LOD[i] <- suppressWarnings(drc::ED(get(LOD.list2[[i]]), 1 - 0.05^0.25, type = "absolute", display = FALSE)[1])
  		DAT3$rep5.LOD[i] <- suppressWarnings(drc::ED(get(LOD.list2[[i]]), 1 - 0.05^0.2, type = "absolute", display = FALSE)[1])
  		DAT3$rep8.LOD[i] <- suppressWarnings(drc::ED(get(LOD.list2[[i]]), 1 - 0.05^0.125, type = "absolute", display = FALSE)[1])
  	}
  	## Generate prediction data for LOQ:
  	## Only get LOQ if LOQ model is determined:
  	if(!is.na(LOQ.list[i])) {
  		newData <- data.frame(Standards = seq(1, 10000))
  		newData$Cq.CV <- predict(get(LOQ.list[[i]]), newData)
  		## Determine what type of LOQ model is used and calculate LOQ accordingly:
  		## For exponential decay:
  		if(as.character(get(LOQ.list[[i]])$call)[1] == "stats::nls") {
  			## Look up lowest modelled standard below the CV threshold:
  			DAT3$LOQ[i] <- min(newData$Standards[newData$Cq.CV <= LOQ.Threshold])
  			## Unless... If the background variation exceeds the CV threshold, adjust threshold:
  			## Determine the highest standard used:
  			A <- max(DAT2$Standards[DAT2$Target == Targets[i]])
  			if(min(newData$Cq.CV[newData$Standards <= A]) > LOQ.Threshold) {
  				## Set the adjusted threshold to 1.5x the lowest simulated Cq.CV
  				##   within the range of data tested:
  				B <- min(newData$Cq.CV[newData$Standards <= A])
  				DAT3$LOQ[i] <- min(newData$Standards[newData$Cq.CV <= B * 1.5])
  				## Make a note of the adjusted threshold in the analysis log:
  				ToWrite <- paste0("Note: All standards tested for ",Targets[i],
  													" yielded higher Cq.CV values than the user-defined CV threshold of ",
  													LOQ.Threshold,". The CV threshold has been adjusted to ",
  													B*1.5," for the LOQ of this marker.")
  				print(ToWrite)
  			}
  		}
  		if(as.character(get(LOQ.list[[i]])$call)[1] == "stats::lm") {
  			## For polynomial:
  			if(grepl("poly", as.character(get(LOQ.list[[i]])$call)[2]) == TRUE) {
  				## Determine the highest standard used:
  				A <- max(DAT2$Standards[DAT2$Target == Targets[i]])
  				## Adjust if the tested range does not cross below the CV threshold:
  				if(min(DAT2$Cq.CV[DAT2$Target == Targets[i]], na.rm = TRUE) > LOQ.Threshold) {
  					B <- min(DAT2$Cq.CV[DAT2$Target == Targets[i]], na.rm = TRUE) * 1.5
  					## Make a note of the adjusted threshold in the analysis log:
  					ToWrite <- paste0("Note: All standards tested for ",Targets[i],
  														" yielded higher Cq.CV values than the user-defined CV threshold of ",
  														LOQ.Threshold,". The CV threshold has been adjusted to ",
  														B," for the LOQ of this marker.")
  					print(ToWrite)
  				}
  				else {
  					B <- LOQ.Threshold
  				}
  				## Look up highest modelled standard below the CV threshold:
  				C <- max(newData$Standards[newData$Cq.CV <= B & newData$Standards <= A])
  				## Look up the highest modelled standard above the CV threshold...
  				##   and also below the highest standard below the CV threshold.
  				##   This captures the farthest right crossing point on a downward slope.
  				D <- max(newData$Standards[newData$Cq.CV > B & newData$Standards < C])
  				## LOQ is D + 1 to get back less than or equal to the CV threshold.
  				DAT3$LOQ[i] <- D + 1
  			}
  			# For linear:
  			else {
  				## Look up lowest modelled standard below the CV threshold:
  				DAT3$LOQ[i] <- min(newData$Standards[newData$Cq.CV <= LOQ.Threshold])
  			}
  		}
  		## If modelled LOQ is calculated to be below the 95% LOD, set LOQ as LOD:
  		if(!is.na(DAT3$LOD[i])) {
  			if(DAT3$LOQ[i] < DAT3$LOD[i]) {
  				DAT3$LOQ[i] <- DAT3$LOD[i]
  			}
  		}
  		## If modelled LOQ is calculated to be below the lowest standard tested,
  		##   set the lowest standard as the LOQ:
  		if(DAT3$LOQ[i] < min(DAT2$Standards[DAT2$Target == Targets[i]])) {
  			DAT3$LOQ[i] <- min(DAT2$Standards[DAT2$Target == Targets[i]])
  		}
  	}
  }
  out.list <- list(dataSum = DAT, standardsSum = DAT2, assaySum = DAT3,
  								 LOQlist = LOQ.out, LODlist = LOD.out)
  return(out.list)
}
