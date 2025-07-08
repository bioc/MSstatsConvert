#' Clean raw Spectronaut output.
#' @param msstats_object an object of class `MSstatsSpectronautFiles`.
#' @param intensity chr, specifies which column will be used for Intensity.
#' @return `data.table`
#' @keywords internal
.cleanRawSpectronaut = function(msstats_object, intensity,
                                calculateAnomalyScores, 
                                anomalyModelFeatures) {
  FFrgLossType = FExcludedFromQuantification = NULL
  
  spec_input = getInputFile(msstats_object, "input")
  .validateSpectronautInput(spec_input)
  spec_input = spec_input[FFrgLossType == "noloss", ]

  f_charge_col = .findAvailable(c("FCharge", "FFrgZ"), colnames(spec_input))
  pg_qval_col = .findAvailable(c("PGQvalue"), colnames(spec_input))
  interference_col = .findAvailable(c("FPossibleInterference"), 
                                    colnames(spec_input))
  exclude_col = .findAvailable(c("FExcludedFromQuantification"), 
                               colnames(spec_input))
  cols = c("PGProteinGroups", "EGModifiedSequence", "FGCharge", "FFrgIon", 
           f_charge_col, "RFileName", "RCondition", "RReplicate", 
           "EGQvalue", pg_qval_col, interference_col, exclude_col,
           paste0("F", intensity))
  if (calculateAnomalyScores){
    cols = c(cols, anomalyModelFeatures)
  }
  cols = intersect(cols, colnames(spec_input))
  spec_input = spec_input[, cols, with = FALSE]
  data.table::setnames(
    spec_input, 
    c("PGProteinGroups", "EGModifiedSequence", "FGCharge", "FFrgIon",
      f_charge_col, "RFileName", paste0("F", intensity), 
      "RCondition", "RReplicate"),
    c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
      "ProductCharge", "Run", "Intensity", "Condition", "BioReplicate"), 
    skip_absent = TRUE)
  .logSuccess("Spectronaut", "clean")
  spec_input
}

#' Helper method to validate input has necessary columns
#' @param spec_input dataframe input
#' @noRd
.validateSpectronautInput = function(spec_input) {
    required_columns = c(
        "FFrgLossType", "FExcludedFromQuantification", "PGProteinGroups",
        "FGCharge")
    missing_columns = setdiff(required_columns, colnames(spec_input))
    if (length(missing_columns) > 0) {
        msg = paste("The following columns are missing from the input data:", 
                    paste(missing_columns, sep = ", ", collapse = ", "))
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
}
