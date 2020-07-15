#' @title Predict copy number from Ct value using one or more standard curves.
#'
#' @description Converts Ct value (or Cq value) to copy number (CN) using one
#'     or more standard curves using inverse prediction. Predictions are back transformed
#'     along with associated standard errors.
#'
#' @param calib_df A \code{\link{data.frame}} object containing calibration curve data.
#' @param ct_df A data.frame object containing Ct value data.
#'
#' @details The \code{calib_df} data.frame object contains data from at least one qPCR calibration curve
#'     usually presented as Cq (Ct) value and corresponding copy number from from a series of serial dilutions.
#'     The data.frame must contain the headers \code{Target}, \code{Cq} and \code{SQ}. The \code{Target}
#'     column must contain unique identifiers for each calibration curve. The \code{Cq} column contains the
#'     Cq (Ct) values and \code{SQ} contains the copy number data. Additional columns will be ignored.
#'
#'     The \code{ct_df} data.frame object usually contains qPCR Cq (Ct) values from sample data. The data.frame must
#'     contain the headers \code{calib_curve} and \code{Ct.value}. The \code{calib_curve} column must contain
#'     at least one unique calibration curve identifier corresponding to those in \code{calib_df}. Any identifiers
#'     found in \code{ct_df} but not in \code{calib_df} will be ignored. The \code{Ct.value} should contain the Cq (Ct)
#'     values. Any additional columns in the data.frame will remain unchanged and be included in the returned object.
#'
#'     Non-detections in either the \code{calib_df} or \code{ct_df} data.frames should be represented as \code{NA}.
#'
#' @return A \code{\link{data.frame}} object containing original data in \code{ct_df} with back transformed
#'     CN predictions and associated standard errors.
#'
#' @export
#' @importFrom chemCal inverse.predict
#' @importFrom plyr ldply
#' @importFrom stats lm
#'
#' @examples
#' \dontrun{
#' predictCN(calib_df = curve_data, ct_df = field_data)
#' }
predictCN <- function(calib_df, ct_df){
	# check correct arguments supplied
	stopifnot(!missing(calib_df) || class(calib_df == "data.frame"))
	stopifnot(!missing(ct_df) || class(ct_df == "data.frame"))

	calib <- calib_df
	dataf <- ct_df
  calib$Target <- as.character(calib$Target)
  dataf$calib_curve <- as.character(dataf$calib_curve)

	# check columns names supplied
	if(sum(colnames(calib)=="Target")!=1) {
		stop(paste0(calib_df, " does not contain Targets column."))
	}
	if(sum(colnames(calib)=="Cq")!=1) {
		stop(paste0(calib_df, " does not contain Cq column."))
	}
	if(sum(colnames(calib)=="SQ")!=1) {
		stop(paste0(calib_df, " does not contain SQ column."))
	}
	if(sum(colnames(dataf)=="calib_curve")!=1) {
		stop(paste0(ct_df, " does not contain calib_curve column."))
	}
	if(sum(colnames(dataf)=="Ct.value")!=1) {
		stop(paste0(ct_df, " does not contain Ct.value column."))
	}

	# check if ct_df calib_curve column contains at
	# least one identifier in calib_df Targets
	if(sum(unique(dataf$calib_curve) %in% unique(calib$Target)) == 0) {
		stop(paste0("no Target identifiers in ", ct_df, " found in ", calib_df))
	}
  calib.targets <- unique(calib$Target)
  dataf.targets <- unique(dataf$calib_curve)
  # check to see what calib_curve identifiers in data_file
  # are in calib_file Target
  in.idx <- dataf.targets[dataf.targets %in% calib.targets]
  out.idx <- dataf.targets[!(dataf.targets %in% calib.targets)]
  if(length(out.idx) != 0){
  	print(paste0("The following calibration curve identifiers were found in ", ct_df, " but not in ", calib_df, " and will be ignored: "))
  	print(out.idx)
  }
	mod.list <- list(0)
	new.data <- data.frame()
  for(i in seq_along(in.idx)){
  	# fit calibration curve linear models
    mod.list[[i]] <- lm(Cq ~ log10(SQ), data = calib, subset = calib$Target == in.idx[i])
    names(mod.list)[[i]] <- paste0("mod_", in.idx[i])
    # create predictions of CN from Ct using model
    dataf.tmp <- dataf[dataf$calib_curve == in.idx[i],]
    mod.name <- paste0("mod_", in.idx[i])
    xpred <- t(sapply(dataf.tmp$Ct.value,function(y) chemCal::inverse.predict(mod.list[[mod.name]], y)[1:2]))
    pred.tmp <- data.frame(plyr::ldply(xpred[,1], data.frame),
    											 abs(plyr::ldply(xpred[,2], data.frame)))
    names(pred.tmp) <- c("CN_pred", "CN_se")
    # backtransform predicted CN values
    pred.tmp$CN_back <- 10^pred.tmp$CN_pred
    pred.tmp$CN_back_se <- 10^pred.tmp$CN_se
    tmp.df <- cbind(dataf.tmp, CN_back = pred.tmp$CN_back, CN_back_se = pred.tmp$CN_back_se)
    new.data <- rbind(new.data, tmp.df)
  }
  return(new.data)
}


