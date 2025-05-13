.checkMissing = function(input){
    
    missing_percent = sum(is.na(input$Intensity)) / nrow(input)
    msg = paste0(as.character(round(missing_percent,4)*100),
                "% of the values are missing across all intensities")
    getOption("MSstatsMsg")("INFO", msg)
    
    return(missing_percent)
}

.checkIntensityDistribution = function(input){
    log_intensities = log2(input$Intensity)
    log_intensities = log_intensities[!is.na(log_intensities)]
    mean_int = mean(log_intensities, na.rm=TRUE) #quantile
    sd_int = sd(log_intensities, na.rm=TRUE)
    
    if (mean_int - (2*sd_int) < 0){
        below_zero_warning = TRUE
        
        msg = paste0("Intensity distribution indicates zero truncated data.", 
                     " More missing values than normal should be expected.")
        getOption("MSstatsMsg")("INFO", msg)
    } else {
        below_zero_warning = FALSE
        
        msg = paste0("Intensity distribution standard truncated data.")
        getOption("MSstatsMsg")("INFO", msg)
    }
    
    return(below_zero_warning)
    
}

.checkFeatureSD = function(input){
    sd_data = input[, .(sd_Intensity = sd(log2(Intensity), na.rm = TRUE),
                        mean_Intensity = mean(log2(Intensity), na.rm = TRUE),
                        ratio = sd(Intensity, na.rm = TRUE) / mean(
                            Intensity, na.rm = TRUE)), 
                    by = Feature]
    
    high_ratio = nrow(sd_data[sd_data$ratio > .5]) / nrow(sd_data) # by condition
    msg = paste0(as.character(round(high_ratio,4)*100),
                  "% of features have a ratio of standard deviation to mean > 0.5.",
                 " A high value here could indicate problematic measurements.")
    getOption("MSstatsMsg")("INFO", msg)
    
    return(sd_data)
}

.checkFeatureOutliers = function(input, feature_data){
    
    full_dt = merge(input, feature_data, by = "Feature", all.x = TRUE)
    full_dt[, is_outlier := abs(
        log2(Intensity) - mean_Intensity) > 2 * sd_Intensity] # switch to quantile
    percent_outlier = sum(full_dt$is_outlier == TRUE, na.rm=TRUE
        ) / nrow(full_dt[!is.na(full_dt$Intensity)])
    
    outlier_info = full_dt[, .(outliers = sum(is_outlier, na.rm=TRUE)), 
                           by = Feature]
    feature_data = merge(feature_data, outlier_info, 
                         by = "Feature", all.x = TRUE)
    
    msg = paste0(as.character(round(percent_outlier,4)*100),
                 "% of all measured intensities fall within the outlier range")
    getOption("MSstatsMsg")("INFO", msg)
    
    return(list(feature_data, percent_outlier))
}

.checkFeatureCoverage = function(input, feature_data){
    n_run = length(unique(input$Run))
    missing_check = input[, .(percent_missing = sum(is.na(Intensity)) / n_run), 
                          by=Feature]
    feature_data = merge(feature_data, missing_check, 
                         by = "Feature", all.x = TRUE)
    
    high_missing = nrow(feature_data[feature_data$percent_missing > .5]
                        ) / nrow(feature_data)
    msg = paste0(as.character(round(high_missing, 4)*100),
                 "% of features are missing in more than half the runs.",
                 " These features may need to be removed.")
    getOption("MSstatsMsg")("INFO", msg)
    
    return(feature_data)
}

pearson_skewness = function(x) {
    n = length(x)
    mean_x = mean(x)
    sd_x = sd(x)  # sample standard deviation
    skewness = sum((x - mean_x)^3) / n / (sd_x^3)
    return(skewness)
}

.checkAnomalySkew = function(input){
    
    input$PSM = paste(input$PeptideSequence,
                      input$PrecursorCharge, sep="_")

    skew_results = input[, .(skew = pearson_skewness(AnomalyScores)), by = PSM]
    
    return(skew_results)
    
}