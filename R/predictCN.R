library(dplyr)
library(chemCal)

#' @title Predict copy number from Ct value using one or more standard curves.
#'
#' @description This function converts Ct value (or Cq value) to copy number (CN) using one 
#'     or more standard curves using inverse prediction.
#'
#' @param calib_df A \code{\link{data.frame}} object containing calibration curve data. 
#'     The data.frame must contain column names 'Target', 'Cq' and 'SQ'. Other columns will remain unchanged.
#' @param ct_df A \code{\link{data.frame}} object containing Ct value data. The data.frame must 
#'     contain column names 'calib_curve' and 'Ct.value'.Other columns will remain unchanged.
#'
#' @return A \code{\link{data.frame}} object containing CN predictions
#' @export
#'
#' @examples
predictCN <- function(calib_df, ct_df){
	# check correct arguments supplied
	stopifnot(!missing(calib_df) || class(calib_df == "character"))
	stopifnot(!missing(ct_df) || class(ct_df == "character"))
	
	# load data
	calib <- read.table(calib_df, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	dataf <- read.table(ct_df, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
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
	
	# check if data_file calib_curve column contains at 
	# least one identifier in calib_file Targets 
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
  	# fit calibration curve models
    mod.list[[i]] <- lm(Cq ~ log10(SQ), data = calib, subset = Target == in.idx[i])
    names(mod.list)[[i]] <- paste0("mod_", in.idx[i])
    # create predictions of CN from Ct using model
    dataf.tmp <- dataf[dataf$calib_curve == in.idx[i],]
    mod.name <- paste0("mod_", in.idx[i])
    xpred <- t(sapply(dataf.tmp$Ct.value,function(y) inverse.predict(mod.list[[mod.name]], y)[1:2]))
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

# test function
my_test <- conv.CN('data/creran_calib_curve.txt', 'data/creran_field_data.txt')
write.table(my_test, "output/creran_field_data_CN_test.csv", sep = ",", row.names = FALSE, col.names = TRUE)

my_neg_test <- conv.CN('data/creran_calib_curve.txt', 'data/test_data.txt')
write.table(my_neg_test, "output/creran_field_data_CN_neg_test.csv", sep = ",", row.names = FALSE, col.names = TRUE)
