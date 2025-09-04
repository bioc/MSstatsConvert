#' Prepare data for anomaly model, including feature engineering
#' @noRd
.prepareSpectronautAnomalyInput = function(input,
                                           quality_metrics,
                                           run_order,
                                           n_feat=100,
                                           missing_run_count=.5){
    input = as.data.table(input)
    
    input$Fragment = paste(input$PeptideSequence,
                           input$PrecursorCharge, 
                           input$FragmentIon,
                           input$ProductCharge, sep="_")
    input$PSM = paste(input$PeptideSequence,
                           input$PrecursorCharge, sep="_")

    # Remove fragments with more than 50% missing
    runs = length(unique(input$Run))
    remove = input[is.na(Intensity), .(count = .N / runs), by = Fragment
                   ][count > missing_run_count, Fragment]
    input = input[!(Fragment %in% remove)]
    
    # Select top N features
    feature_counts = input[
        , .(MeanAbundance = mean(as.numeric(Intensity), na.rm = TRUE)),
        by = .(ProteinName, Fragment)]
    feature_counts[, feature_rank := frank(-MeanAbundance, ties.method = "min"),
                   by = ProteinName]
    feature_counts = feature_counts[feature_rank <= n_feat,
                                    .(ProteinName, Fragment)]
    
    input = input[feature_counts, on = .(ProteinName, Fragment)]
    
    # Model at PSM level
    input = merge(input, run_order, by="Run", 
                  all.x=TRUE, all.y=FALSE)
    cols=c("ProteinName", "PSM", "Run", "Order", quality_metrics)
    feature_data = unique(
        input[, ..cols])
    setorder(feature_data, PSM, Order)
    
    # Scale quality metrics
    feature_data[, (quality_metrics) := lapply(
        .SD, function(x) as.numeric(scale(x))), 
        by = PSM, .SDcols = quality_metrics]
    
    # Apply add_features to each metric (avoids for loop)
    feature_data = feature_data[, {
        feature_list = lapply(.SD, .add_features)  # Apply add_features to each column
        feature_data = do.call(cbind, feature_list)  # Combine results
        feature_data[, Order := Order] # Ensure order is returned for join
        feature_data
    }, by = PSM, .SDcols = quality_metrics]
    
    # Add metrics back to df
    input = merge(input, feature_data,
                       all.x=TRUE, all.y=FALSE, 
                       by=c("PSM", "Order"))
    
    return(input)
}

#' Add moving window features given a quality metric vector
#' @noRd
.add_features = function(quality_vector){
    
    mean_increase = .add_mean_increase(quality_vector)
    mean_decrease = .add_mean_decrease(quality_vector)
    dispersion_increase = .add_dispersion_increase(quality_vector)

    
    feature_df = data.table(
        "mean_increase"=mean_increase,
        "mean_decrease"=mean_decrease,
        "dispersion_increase"=dispersion_increase
    )
    return(feature_df)
}

#' Calculate mean increase
#' @noRd
.add_mean_increase = function(quality_vector){
    
    mean_increase = numeric(length(quality_vector))
    mean_increase[1] = 0
    d = 0.5
    
    for(k in 2:length(quality_vector)) {
        # 5 is reference (3 sigma)
        if (mean_increase[k] > 5){
            mean_increase[k] = max(0,(quality_vector[k] - d), na.rm = TRUE)
        } else {
            mean_increase[k] = max(0,
                                   (quality_vector[k] - d + mean_increase[k-1]),
                                   na.rm = TRUE) # positive CuSum
        }
    }
    return(mean_increase)
}

#' Calculate mean decrease
#' @noRd
.add_mean_decrease = function(quality_vector){
    
    mean_decrease = numeric(length(quality_vector))
    mean_decrease[1] = 0
    d = -0.5
    
    for(k in 2:length(quality_vector)) {
        # 5 is reference (3 sigma)
        if (mean_decrease[k]>5){
            mean_decrease[k] <- max(0,(d - quality_vector[k] + 0), na.rm = TRUE)
        } else {
            mean_decrease[k] <- max(0,
                                    (d - quality_vector[k] + mean_decrease[k-1]),
                                    na.rm = TRUE) # negative CuSum
        }
    }
    return(mean_decrease)
}

#' Calculate dispersion increase
#' @noRd
.add_dispersion_increase = function(quality_vector){
    dispersion_increase = numeric(length(quality_vector))
    v = numeric(length(quality_vector))
    v[1] = (sqrt(abs(quality_vector[1]))-0.822)/0.349
    d = 0.5
    for(k in 2:length(quality_vector)) {
        
        v[k] = (sqrt(abs(quality_vector[k]))-0.822)/0.349 
        
        if (dispersion_increase[k] > 5){
            dispersion_increase[k] = max(0,(v[k] - d),
                                         na.rm = TRUE)
        } else {
            dispersion_increase[k] = max(0, 
                                         (v[k] - d + dispersion_increase[k-1]),
                                         na.rm = TRUE) # CuSum variance
        }
    }
    return(dispersion_increase)
}

#' Train isolation forest model in parallel
#' @import parallel
#' @import Rcpp
#' 
#' @noRd
.runAnomalyModel = function(input_data, 
                            n_trees, 
                            max_depth, 
                            cores,
                            split_column,
                            quality_metrics){
    
    function_environment = environment()
    cl = parallel::makeCluster(cores)
    
    parallel::clusterExport(cl, c(
        "calculate_anomaly_score",
        "as.data.table",
        "max_depth",
        "n_trees",
        "split_column",
        "quality_metrics"), 
        envir = function_environment)
    
    psm_list = unique(input_data[, get(split_column)])
    num_psm = length(psm_list)
    
    cat(paste0("Number of PSMs to process: ", num_psm), 
        sep = "\n", file = "MSstats_anomaly_model_progress.log")
    
    model_results = parallel::parLapply(
        cl, seq_len(num_psm), 
        function(i){
            if (i %% 100 == 0) {
                cat("Finished processing an additional 100 PSMs", 
                    sep = "\n", file = "MSstats_anomaly_model_progress.log", append = TRUE)
            }
            single_psm = input_data[get(split_column) == psm_list[[i]], 
                                    ..quality_metrics]
            
            if (max_depth == "auto"){
                max_depth = round(log2(nrow(single_psm)))
            }
            
            forest = calculate_anomaly_score(
                single_psm, n_trees, max_depth)
            forest$anomaly_score
        }
    )
    
    model_results = unlist(model_results)
    # Clip anomaly scores to stop them from exploding
    input_data$AnomalyScores = pmax(model_results, .001)
    
    return(input_data)
}
