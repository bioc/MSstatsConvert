#' Check validity of parameters to the `MSstatsImport` function.
#' @inheritParams MSstatsImport
#' @return none, throws an error if any of the assertions fail
#' @keywords internal
.checkMSstatsParams = function(input, annotation, 
                               feature_columns,
                               remove_shared_peptides,
                               remove_single_feature_proteins,
                               feature_cleaning) {
    checkmate::assertDataTable(input)
    checkmate::assertTRUE(is.character(annotation) | 
                              inherits(annotation, "data.frame") | 
                              is.null(annotation))
    checkmate::assertCharacter(feature_columns, min.len = 1)
    checkmate::assertLogical(remove_shared_peptides)
    checkmate::assertLogical(
        feature_cleaning[["remove_features_with_few_measurements"]], 
        )
    checkmate::assertLogical(remove_single_feature_proteins)
    checkmate::assertLogical(feature_cleaning[["remove_psms_with_all_missing"]],
                             null.ok = TRUE)
    checkmate::assertFunction(feature_cleaning[["summarize_multiple_psms"]])
}


#' Select an available options from a set of possibilities
#' @param possibilities possible legal values of a variable
#' @param option_set set of values that includes one of the `possibilities`
#' @param fall_back if there is none of the `possibilities` in `option_set`,
#' or there are multiple hits, default to `fall_back`
#' @return same as option_set, usually character
#' @keywords internal
.findAvailable = function(possibilities, option_set, fall_back = NULL) {
    chosen = option_set[option_set %in% possibilities]
    if (length(chosen) != 1L) {
        if (is.null(fall_back)) {
            NULL
        } else {
            if (fall_back %in% possibilities) {
                fall_back
            } else {
                NULL
            }
        }
    } else {
        chosen
    }
}


#' Generic parameter validation for all MSstats converters
#' @param input input data - required for all converters
#' @param annotation annotation data - required for all converters
#' @param intensity intensity type - common parameter
#' @param excludedFromQuantificationFilter filter setting - converter specific
#' @param filter_with_Qvalue Q-value filter setting - common parameter
#' @param qvalue_cutoff Q-value cutoff - common parameter
#' @param useUniquePeptide unique peptide setting - common parameter
#' @param removeFewMeasurements remove few measurements setting - common parameter
#' @param removeProtein_with1Feature remove single feature proteins setting - common parameter
#' @param summaryforMultipleRows aggregation function - common parameter
#' @param calculateAnomalyScores anomaly detection setting - converter specific
#' @param anomalyModelFeatures anomaly model features - converter specific
#' @param anomalyModelFeatureTemporal temporal features - converter specific
#' @param removeMissingFeatures missing feature threshold - converter specific
#' @param anomalyModelFeatureCount feature count for anomaly model - converter specific
#' @param runOrder run order data - converter specific
#' @param n_trees number of trees - converter specific
#' @param max_depth max tree depth - converter specific
#' @param numberOfCores number of cores - common parameter
#' @param use_log_file logging setting - common parameter
#' @param append append setting - common parameter
#' @param verbose verbose setting - common parameter
#' @param log_file_path log file path - common parameter
#' @param converter_name name of the converter for error messages
#' @param ... additional converter-specific parameters
#' @return NULL (throws error if validation fails)
#' @keywords internal
.validateMSstatsConverterParameters <- function(
        input, 
        annotation = NULL,
        intensity = NULL,
        excludedFromQuantificationFilter = NULL,
        filter_with_Qvalue = FALSE,
        qvalue_cutoff = 0.01,
        useUniquePeptide = TRUE,
        removeFewMeasurements = TRUE,
        removeProtein_with1Feature = FALSE,
        summaryforMultipleRows = max,
        calculateAnomalyScores = FALSE,
        anomalyModelFeatures = c(),
        anomalyModelFeatureTemporal = c(),
        removeMissingFeatures = 0.5,
        anomalyModelFeatureCount = 100,
        runOrder = NULL,
        n_trees = 100,
        max_depth = "auto",
        numberOfCores = 1,
        use_log_file = TRUE,
        append = FALSE,
        verbose = TRUE,
        log_file_path = NULL,
        converter_name = "MSstats converter",
        ...
) {
    
    # Input data validation (fail immediately if data is invalid)
    if (is.null(input)) {
        stop("Input data cannot be NULL")
    }
    
    if (is.character(input)) {
        if (!file.exists(input)) {
            stop("Input file does not exist: ", input)
        }
        # Quick file size check for very large files
        file_size_mb <- file.size(input) / (1024^2)
        if (file_size_mb > 1000) {  # Warn for files > 1GB
            warning("Large input file detected (", round(file_size_mb, 1), " MB). ",
                    "Consider validating parameters on a subset first.")
        }
    } else if (is.data.frame(input) || data.table::is.data.table(input)) {
        # Quick structural validation
        if (nrow(input) == 0) {
            stop("Input data is empty (0 rows)")
        }
        if (ncol(input) == 0) {
            stop("Input data has no columns")
        }
    } else {
        stop("Input must be a file path, data.frame, or data.table")
    }
    
    # Annotation validation
    if (!is.null(annotation)) {
        if (is.character(annotation)) {
            if (!file.exists(annotation)) {
                stop("Annotation file does not exist: ", annotation)
            }
        } else if (!is.data.frame(annotation) && !data.table::is.data.table(annotation)) {
            stop("Annotation must be NULL, a file path, data.frame, or data.table")
        }
    }
    
    # Intensity validation (if provided)
    if (!is.null(intensity)) {
        checkmate::assertString(intensity)
    }
    
    # Q-value filtering parameters
    checkmate::assertLogical(filter_with_Qvalue, len = 1)
    checkmate::assertNumber(qvalue_cutoff, lower = 0, upper = 1)
    
    # Common processing parameters
    checkmate::assertLogical(useUniquePeptide, len = 1)
    checkmate::assertLogical(removeFewMeasurements, len = 1)
    checkmate::assertLogical(removeProtein_with1Feature, len = 1)
    
    # Aggregation function validation
    if (!is.null(summaryforMultipleRows)) {
        checkmate::assertFunction(summaryforMultipleRows)
    }
    
    # Core system parameters
    checkmate::assertInt(numberOfCores, lower = 1)
    checkmate::assertLogical(use_log_file, len = 1)
    checkmate::assertLogical(append, len = 1)
    checkmate::assertLogical(verbose, len = 1)
    
    # Converter-specific boolean parameters (if provided)
    if (!is.null(excludedFromQuantificationFilter)) {
        checkmate::assertLogical(excludedFromQuantificationFilter, len = 1)
    }
    
    checkmate::assertLogical(calculateAnomalyScores, len = 1)
    
    # Anomaly detection parameter validation (converter-specific)
    if (calculateAnomalyScores) {
        # These validations only matter if anomaly detection is enabled
        if (length(anomalyModelFeatures) == 0) {
            stop("anomalyModelFeatures cannot be empty when calculateAnomalyScores=TRUE")
        }
        checkmate::assertCharacter(anomalyModelFeatures, min.len = 1)
        
        if (length(anomalyModelFeatureTemporal) > 0) {
            if (length(anomalyModelFeatureTemporal) != length(anomalyModelFeatures)) {
                stop("anomalyModelFeatureTemporal must have same length as anomalyModelFeatures or be empty")
            }
            valid_temporal <- c("mean_decrease", "mean_increase", "dispersion_increase")
            invalid_temporal <- anomalyModelFeatureTemporal[
                !is.null(anomalyModelFeatureTemporal) & 
                    !anomalyModelFeatureTemporal %in% valid_temporal
            ]
            if (length(invalid_temporal) > 0) {
                stop("Invalid temporal directions: ", paste(invalid_temporal, collapse = ", "),
                     ". Must be one of: ", paste(valid_temporal, collapse = ", "), " or NULL")
            }
        }
        
        checkmate::assertInt(n_trees, lower = 1)
        if (is.character(max_depth)) {
            checkmate::assertChoice(max_depth, choices = "auto")
        } else {
            checkmate::assertInt(max_depth, lower = 1)
        }
        checkmate::assertInt(anomalyModelFeatureCount, lower = 1)
        checkmate::assertNumber(removeMissingFeatures, lower = 0, upper = 1)
        
        if (!is.null(runOrder)) {
            if (!is.data.frame(runOrder) && !data.table::is.data.table(runOrder)) {
                stop("runOrder must be a data.frame or data.table")
            }
            required_cols <- c("Run", "Order")
            missing_cols <- setdiff(required_cols, colnames(runOrder))
            if (length(missing_cols) > 0) {
                stop("runOrder is missing required columns: ", paste(missing_cols, collapse = ", "))
            }
            if (!is.numeric(runOrder$Order)) {
                stop("runOrder$Order must be numeric")
            }
        }
    } else {
        # When anomaly detection is disabled, these parameters are ignored but warn user
        if (length(anomalyModelFeatures) > 0) {
            warning("anomalyModelFeatures provided but calculateAnomalyScores=FALSE, ignoring")
        }
        if (!is.null(runOrder)) {
            warning("runOrder provided but calculateAnomalyScores=FALSE, ignoring")
        }
    }
    
    # Log file validation
    if (!is.null(log_file_path)) {
        checkmate::assertString(log_file_path)
        log_dir <- dirname(log_file_path)
        if (!dir.exists(log_dir)) {
            stop("Log file directory does not exist: ", log_dir)
        }
        if (!file.access(log_dir, mode = 2) == 0) {
            stop("No write permission for log file directory: ", log_dir)
        }
    }
}