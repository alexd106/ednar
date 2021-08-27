#' @title Predict Copy Number From Ct Values Using qPCR Standard Curves
#'
#' @description Predicts copy number (CN) from Ct (Cq) values from one
#'     or more standard curves using inverse prediction.
#'
#' @param calib_df A \code{\link{data.frame}} or \code{\link{tibble}} object containing calibration curve data.
#' @param ct_df A \code{\link{data.frame}} or \code{\link{tibble}} object containing Ct values from environmental samples.
#' @param ... Placeholder for further arguments that might be needed by future implementations.
#'
#' @details The \code{calib_df} object contains data for at least one qPCR calibration curve
#'     usually presented as Cq (Ct) value and corresponding copy number from from a series of serial dilutions.
#'     The \code{data.frame} must contain the headers \code{Target}, \code{Cq} and \code{SQ}. The \code{Target}
#'     column must contain unique identifiers for each calibration curve. The \code{Cq} column contains the
#'     Cq (Ct) values and \code{SQ} contains the copy number data. Additional columns will be ignored.
#'
#'     The \code{ct_df} data.frame object usually contains qPCR Cq (Ct) values from sample data. The data.frame must
#'     contain the headers \code{calib.curve} and \code{Ct.value}. The \code{calib.curve} column must contain
#'     at least one unique calibration curve identifier corresponding to those in \code{calib_df}. Any identifiers
#'     found in \code{ct_df} but not in \code{calib_df} will be ignored and a \code{warning} displayed. The
#'     \code{Ct.value} column should contain the Cq (Ct) values from environmental samples. Any additional columns
#'     in the data.frame will remain unchanged and be included in the returned object.
#'
#'     Non-detections in either the \code{calib_df} or \code{ct_df} \code{data.frames} should be represented as \code{NA}.
#'
#'     Copy number is predicted from Cq (Ct) values by fitting a linear model (\code{\link{lm}}) to data from each
#'     calibration curve and then using the \code{\link{inverse.predict}} function from the \code{chemCal} package.
#'     package.
#'
#' @note Copy number values in \code{calib_df} should not be log transformed prior to using this function.
#'
#' @references Massart, L.M, Vandenginste, B.G.M., Buydens, L.M.C., De Jong, S., Lewi, P.J., Smeyers-Verbeke, J. (1997)
#'      Handbook of Chemometrics and Qualimetrics: Part A, p. 200.
#'
#' @return A \code{\link{tibble}} object containing original data in \code{ct_df} with back transformed
#'     copy number predictions (\code{CN.back}) and associated standard errors (\code{CN.back.se}).
#'
#' @export
#' @importFrom chemCal inverse.predict
#' @importFrom plyr ldply
#' @importFrom stats lm
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' my_pred <- calib_predict(calib_df = calib_data, ct_df = field_data)
#' }
calib_predict <- function(calib_df, ct_df, ...){
	# check correct arguments supplied
	stopifnot(!missing(calib_df) || class(calib_df) == "data.frame")
	stopifnot(!missing(ct_df) || class(ct_df) == "data.frame")

	calib_df_name <- deparse(substitute(calib_df))
	ct_df_name <- deparse(substitute(ct_df))

	calib <- calib_df
	dataf <- ct_df

	# check columns names supplied
	if(!any(colnames(calib) == "Target")) {
		stop(paste0("can't find column `Target` in ", calib_df_name, "."))
	}
	if(!any(colnames(calib) == "Cq")) {
		stop(paste0("can't find variable `Cq` in ", calib_df_name, "."))
	}
	if(!any(colnames(calib) == "SQ")) {
		stop(paste0("can't find variable `SQ` in ", calib_df_name, "."))
	}
	if(!any(colnames(dataf) == "calib.curve")) {
		stop(paste0("can't find variable `calib.curve` in ", ct_df_name, "."))
	}
	if(!any(colnames(dataf) == "Ct.value")) {
		stop(paste0("can't find variable `Ct.value` in ", ct_df_name, "."))
	}

	calib$Target <- as.character(calib$Target)
	dataf$calib.curve <- as.character(dataf$calib.curve)

	# check if ct_df calib.curve column contains at
	# least one identifier in calib_df Targets
	if(sum(unique(dataf$calib.curve) %in% unique(calib$Target)) == 0) {
		stop(paste0("no `calib.curve` identifiers in ", ct_df_name, " found in `Targets` ", calib_df_name, "."))
	}
  calib.targets <- unique(calib$Target)
  dataf.targets <- unique(dataf$calib.curve)
  # check to see what calib.curve identifiers in ct_df
  # are in calib_df Target
  in.idx <- dataf.targets[dataf.targets %in% calib.targets]
  out.idx <- dataf.targets[!(dataf.targets %in% calib.targets)]
  if(length(out.idx) > 0){
  	print(paste0("The following calibration curve identifiers were found in ", ct_df_name, " but not in ", calib_df_name, " and will be ignored: "))
  	print(out.idx)
  }
	mod.list <- list(0)
	new.data <- data.frame()
  for(i in seq_along(in.idx)){
  	# fit calibration curve linear models
    mod.list[[i]] <- lm(Cq ~ log10(SQ), data = calib, subset = calib$Target == in.idx[i])
    names(mod.list)[[i]] <- paste0("mod_", in.idx[i])
    # create predictions of CN from Ct using model
    dataf.tmp <- dataf[dataf$calib.curve == in.idx[i],]
    mod.name <- paste0("mod_", in.idx[i])
    xpred <- t(sapply(dataf.tmp$Ct.value, function(y) chemCal::inverse.predict(mod.list[[mod.name]], y)[1:2]))
    pred.tmp <- data.frame(plyr::ldply(xpred[,1], data.frame),
    											 abs(plyr::ldply(xpred[,2], data.frame)))
    names(pred.tmp) <- c("CN.pred", "CN.se")
    # backtransform predicted CN values
    pred.tmp$CN.back <- 10^pred.tmp$CN.pred
    pred.tmp$CN.back.se <- 10^pred.tmp$CN.se
    tmp.df <- cbind(dataf.tmp, CN.back = pred.tmp$CN.back, CN.back.se = pred.tmp$CN.back.se)
    new.data <- rbind(new.data, tmp.df)
  }
  return(tibble(new.data))
}


