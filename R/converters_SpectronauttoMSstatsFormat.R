#' Import Spectronaut files
#' 
#' @param input name of Spectronaut output, which is long-format. ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity, F.ExcludedFromQuantification are required. Rows with F.ExcludedFromQuantification=True will be removed.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Spectronaut, use annotation=NULL (default). It will use the annotation information from input.
#' @param intensity 'PeakArea'(default) uses not normalized peak area. 'NormalizedPeakArea' uses peak area normalized by Spectronaut.
#' @param excludedFromQuantificationFilter Remove rows with F.ExcludedFromQuantification=TRUE Default is TRUE.
#' @param possibleInterferenceFilter Remove rows with F.PossibleInterference=TRUE Default is TRUE.
#' @param filter_with_Qvalue FALSE(default) will not perform any filtering. TRUE will filter out the intensities that have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .sharedParametersAmongConverters
#' 
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
#' @examples 
#' spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
#'                               package = "MSstatsConvert")
#' spectronaut_raw = data.table::fread(spectronaut_raw)
#' spectronaut_imported = SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE)
#' head(spectronaut_imported)
#' 
SpectronauttoMSstatsFormat = function(
        input, annotation = NULL, intensity = 'PeakArea', 
        excludedFromQuantificationFilter = TRUE,
        possibleInterferenceFilter = TRUE,
        filter_with_Qvalue = FALSE, qvalue_cutoff = 0.01, 
        useUniquePeptide = TRUE, removeFewMeasurements=TRUE,
        removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
        calculateAnomalyScores=FALSE, anomalyModelFeatures=c(),
        anomalyModelFeatureTemporal=c(),
        runOrder=NULL, n_trees=100, max_depth="auto", numberOfCores=1, 
        use_log_file = TRUE, append = FALSE, verbose = TRUE, 
        log_file_path = NULL, ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "Spectronaut", ...)

    input = MSstatsConvert::MSstatsClean(input, intensity = intensity,
                                         calculateAnomalyScores, 
                                         anomalyModelFeatures)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    pq_filter = list(score_column = "PGQvalue", 
                     score_threshold = 0.01, 
                     direction = "smaller", 
                     behavior = "fill", 
                     handle_na = "keep", 
                     fill_value = NA_real_,
                     filter = filter_with_Qvalue, 
                     drop_column = TRUE)
    qval_filter = list(score_column = "EGQvalue", 
                       score_threshold = qvalue_cutoff, 
                       direction = "smaller", 
                       behavior = "fill", 
                       handle_na = "keep", 
                       fill_value = NA_real_, 
                       filter = filter_with_Qvalue, 
                       drop_column = TRUE)
    # excluded_interference_filter = list(
    #     col_name = "FPossibleInterference", 
    #     filter_symbols = TRUE, 
    #     behavior = "remove", 
    #     fill_value = NA_real_, 
    #     filter = TRUE, 
    #     drop_column = TRUE)
    excluded_quant_filter = list(
        col_name = "FExcludedFromQuantification", 
        filter_symbols = TRUE, 
        behavior = "fill", 
        fill_value = NA_real_, 
        filter = TRUE, 
        drop_column = TRUE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", 
                        "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(remove_features_with_few_measurements = removeFewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        score_filtering = list(pgq = pq_filter, 
                               psm_q = qval_filter),
        exact_filtering = list(#interference = excluded_interference_filter,
                               excluded_quant = excluded_quant_filter),
        columns_to_fill = list("IsotopeLabelType" = "L"),
        anomaly_metrics = anomalyModelFeatures)
    input[, Intensity := ifelse(Intensity == 0, NA, Intensity)]
    
    input = MSstatsConvert::MSstatsBalancedDesign(
        input, feature_columns, 
        remove_few = removeFewMeasurements,
        anomaly_metrics = anomalyModelFeatures)
    
    if (calculateAnomalyScores){
        input = MSstatsConvert::MSstatsAnomalyScores(
            input, anomalyModelFeatures, anomalyModelFeatureTemporal,
            runOrder, n_trees, max_depth, numberOfCores)
    }
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}