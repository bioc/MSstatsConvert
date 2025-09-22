# Test data setup
create_test_data = function() {
    data.table::data.table(
        ProteinName = c(rep("Q9UFW8", 10), rep("Q96S19", 15)),
        PeptideSequence = c(rep("AEFEEQNVR", 5), rep("TALYVTPLDR", 5),
                            rep("AFPLAEWQPSDVDQR", 5), rep("ASGLLLER", 5),
                            rep("LowAbundancePeptide", 5)),
        PrecursorCharge = rep(2, 25),
        FragmentIon = rep("y3", 25),
        ProductCharge = rep(1, 25),
        IsotopeLabelType = rep("L", 25),
        Condition = rep(c("A", "A", "A", "B", "B"), 5),
        BioReplicate = rep(seq(1:5), 5),
        Run = rep(paste0("Run", seq(1:5)), 5),
        Intensity = c(1000, 1500, 2000, 2500, 3000,
                      1100, 1600, 2100, 2600, 3100,
                      1200, 1700, 2200, 2000, 1700,
                      1300, 1800, 1200, 2800, 1800,
                      100, 200, 300, 400, 500),
        AnomalyScores = rep(c(0.01,0.01,0.01,0.7,0.5), 5)
    )
}

# Test .checkMissing function
test_data = create_test_data()

# Test with no missing values
result = MSstatsConvert:::.checkMissing(test_data)
expect_equal(result, 0, info = ".checkMissing returns 0 for no missing values")
expect_true(is.numeric(result), info = ".checkMissing returns numeric value")

# Test with some missing values
test_data_missing = test_data
test_data_missing$Intensity[1:5] = NA
result_missing = MSstatsConvert:::.checkMissing(test_data_missing)
expect_equal(result_missing, 5/25, info = ".checkMissing correctly calculates missing percentage")

# Test with all missing values
test_data_all_missing = test_data
test_data_all_missing$Intensity = NA
result_all_missing = MSstatsConvert:::.checkMissing(test_data_all_missing)
expect_equal(result_all_missing, 1, info = ".checkMissing returns 1 for all missing values")

# Test .checkIntensityDistribution function
test_data = create_test_data()

# Test normal distribution
result = MSstatsConvert:::.checkIntensityDistribution(test_data)
expect_true(is.logical(result), info = ".checkIntensityDistribution returns logical value")

# Test with low intensities that would trigger zero truncation warning
test_data_low = test_data
test_data_low$Intensity = c(rep(1, 10), rep(2, 15))
result_low = MSstatsConvert:::.checkIntensityDistribution(test_data_low)
expect_true(is.logical(result_low), info = ".checkIntensityDistribution handles low intensity data")

# Test .checkFeatureSD function
test_data = create_test_data()
test_data$Feature = paste(test_data$PeptideSequence,
                           test_data$PrecursorCharge,
                           test_data$FragmentIon,
                           test_data$ProductCharge, sep="_")

result = MSstatsConvert:::.checkFeatureSD(test_data)

# Check structure
expect_true(inherits(result, "data.table"), info = ".checkFeatureSD returns data.table")
expect_true("Feature" %in% names(result), info = ".checkFeatureSD includes Feature column")
expect_true("sd_Intensity" %in% names(result), info = ".checkFeatureSD includes sd_Intensity column")
expect_true("mean_Intensity" %in% names(result), info = ".checkFeatureSD includes mean_Intensity column")
expect_true("ratio" %in% names(result), info = ".checkFeatureSD includes ratio column")

# Check that we get one row per unique feature
expected_features = unique(test_data$Feature)
expect_equal(nrow(result), length(expected_features), info = ".checkFeatureSD returns one row per feature")

# Check that calculations make sense
expect_true(all(result$mean_Intensity > 0, na.rm = TRUE), info = ".checkFeatureSD mean intensities are positive")
expect_true(all(result$sd_Intensity >= 0, na.rm = TRUE), info = ".checkFeatureSD standard deviations are non-negative")
expect_true(all(result$ratio >= 0, na.rm = TRUE), info = ".checkFeatureSD ratios are non-negative")

# Test .checkFeatureOutliers function
test_data = create_test_data()
test_data$Feature = paste(test_data$PeptideSequence,
                           test_data$PrecursorCharge,
                           test_data$FragmentIon,
                           test_data$ProductCharge, sep="_")

feature_data = MSstatsConvert:::.checkFeatureSD(test_data)
result = MSstatsConvert:::.checkFeatureOutliers(test_data, feature_data)

# Check structure
expect_true(is.list(result), info = ".checkFeatureOutliers returns list")
expect_equal(length(result), 2, info = ".checkFeatureOutliers returns list of length 2")

enhanced_feature_data = result[[1]]
outlier_percent = result[[2]]

expect_true("outliers" %in% names(enhanced_feature_data), info = ".checkFeatureOutliers adds outliers column")
expect_true(inherits(enhanced_feature_data, "data.table"), info = ".checkFeatureOutliers returns data.table as first element")
expect_true(is.numeric(outlier_percent), info = ".checkFeatureOutliers returns numeric outlier percentage")
expect_true(outlier_percent >= 0 && outlier_percent <= 1, info = ".checkFeatureOutliers outlier percentage is between 0 and 1")

# Test .checkFeatureCoverage function
test_data = create_test_data()
test_data$Feature = paste(test_data$PeptideSequence,
                           test_data$PrecursorCharge,
                           test_data$FragmentIon,
                           test_data$ProductCharge, sep="_")

feature_data = MSstatsConvert:::.checkFeatureSD(test_data)
result = MSstatsConvert:::.checkFeatureCoverage(test_data, feature_data)

# Check structure
expect_true(inherits(result, "data.table"), info = ".checkFeatureCoverage returns data.table")
expect_true("percent_missing" %in% names(result), info = ".checkFeatureCoverage includes percent_missing column")

# With no missing values, percent_missing should be 0
expect_true(all(result$percent_missing == 0), info = ".checkFeatureCoverage shows 0% missing for complete data")

# Test with missing values
test_data_missing = test_data
test_data_missing$Intensity[1:3] = NA
result_missing = MSstatsConvert:::.checkFeatureCoverage(test_data_missing, feature_data)

# First feature should have some missing percentage
first_feature = result_missing$Feature[1]
first_feature_missing = result_missing[Feature == first_feature]$percent_missing
expect_true(first_feature_missing > 0, info = ".checkFeatureCoverage detects missing values")

# Test pearson_skewness function
# Test symmetric distribution (should be near 0)
symmetric_data = c(1, 2, 3, 4, 4, 3, 2, 1, 5)
skew_symmetric = MSstatsConvert:::pearson_skewness(symmetric_data)
expect_true(abs(skew_symmetric) < 0.5, info = "pearson_skewness returns near-zero for symmetric data")

# Test right-skewed distribution (should be positive)
right_skewed = c(1, 1, 1, 1, 2, 2, 3, 4, 10)
skew_right = MSstatsConvert:::pearson_skewness(right_skewed)
expect_true(skew_right > 0, info = "pearson_skewness returns positive for right-skewed data")

# Test left-skewed distribution (should be negative)
left_skewed = c(10, 10, 10, 10, 9, 9, 8, 7, 1)
skew_left = MSstatsConvert:::pearson_skewness(left_skewed)
expect_true(skew_left < 0, info = "pearson_skewness returns negative for left-skewed data")

# Test .checkAnomalySkew function
test_data = create_test_data()
result = MSstatsConvert:::.checkAnomalySkew(test_data)

# Check structure
expect_true(inherits(result, "data.table"), info = ".checkAnomalySkew returns data.table")
expect_true("PSM" %in% names(result), info = ".checkAnomalySkew includes PSM column")
expect_true("skew" %in% names(result), info = ".checkAnomalySkew includes skew column")

# Check that we get one row per unique PSM
expected_psms = unique(paste(test_data$PeptideSequence, 
                              test_data$PrecursorCharge, sep="_"))
expect_equal(nrow(result), length(expected_psms), info = ".checkAnomalySkew returns one row per PSM")

# Check that skewness values are numeric and finite (FAIL)
expect_true(all(is.numeric(result$skew)), info = ".checkAnomalySkew returns numeric skewness values")
expect_true(all(is.finite(result$skew)), info = ".checkAnomalySkew returns finite skewness values")

# Verify PSM creation
expected_psm_1 = paste("AEFEEQNVR", "2", sep="_")
expect_true(expected_psm_1 %in% result$PSM, info = ".checkAnomalySkew creates correct PSM identifiers")

# Test CheckDataHealth main function
test_data = create_test_data()
result = CheckDataHealth(test_data)

# Check structure
expect_true(is.list(result), info = "CheckDataHealth returns list")
expect_equal(length(result), 2, info = "CheckDataHealth returns list of length 2")

feature_data = result[[1]]
skew_results = result[[2]]

# Check feature_data structure
expect_true(inherits(feature_data, "data.table"), info = "CheckDataHealth returns feature_data as data.table")
expected_cols = c("Feature", "sd_Intensity", "mean_Intensity", 
                   "ratio", "outliers", "percent_missing")
expect_true(all(expected_cols %in% names(feature_data)), info = "CheckDataHealth feature_data has all expected columns")

# Check skew_results structure
expect_true(inherits(skew_results, "data.table"), info = "CheckDataHealth returns skew_results as data.table")
expect_true(all(c("PSM", "skew") %in% names(skew_results)), info = "CheckDataHealth skew_results has expected columns")

# Verify Feature creation in main function
expected_feature = paste("AEFEEQNVR", "2", "y3", "1", sep="_")
expect_true(expected_feature %in% feature_data$Feature, info = "CheckDataHealth creates correct Feature identifiers")

# Test with all same intensities (zero variance)
test_data_same = create_test_data()
test_data_same$Intensity = 1000
result_same = CheckDataHealth(test_data_same)
expect_true(is.list(result_same), info = "CheckDataHealth handles zero variance data")

# Check that SD is 0 for features with same intensities
feature_data_same = result_same[[1]]
expect_true(all(feature_data_same$sd_Intensity == 0), info = "CheckDataHealth correctly calculates zero SD for identical intensities")

# Test CheckDataHealth with different anomaly score patterns
test_data = create_test_data()

# Test with highly skewed anomaly scores
test_data_skewed = test_data
test_data_skewed$AnomalyScores = c(rep(0.01, 20), rep(0.99, 5))
result_skewed = CheckDataHealth(test_data_skewed)

skew_results = result_skewed[[2]]
expect_true(inherits(skew_results, "data.table"), info = "CheckDataHealth handles skewed anomaly scores")
expect_true(all(c("PSM", "skew") %in% names(skew_results)), info = "CheckDataHealth skew results maintain structure with skewed data")

# Test with uniform anomaly scores
test_data_uniform = test_data
test_data_uniform$AnomalyScores = 0.5
result_uniform = CheckDataHealth(test_data_uniform)

skew_uniform = result_uniform[[2]]
# Uniform distribution should have skewness near 0 (but will be NaN for truly uniform)
expect_true(all(is.finite(skew_uniform$skew) | is.nan(skew_uniform$skew)), info = "CheckDataHealth handles uniform anomaly scores appropriately")

# Test additional edge case: different anomaly score distributions per PSM
test_data_mixed = create_test_data()
# Create more varied anomaly scores to test skewness calculation
test_data_mixed$AnomalyScores[1:5] = c(0.1, 0.2, 0.3, 0.8, 0.9)  # Right skewed for first PSM
test_data_mixed$AnomalyScores[6:10] = c(0.9, 0.8, 0.3, 0.2, 0.1)  # Left skewed for second PSM
result_mixed = CheckDataHealth(test_data_mixed)

skew_mixed = result_mixed[[2]]
expect_equal(nrow(skew_mixed), 5, info = "CheckDataHealth calculates skewness for all unique PSMs")
expect_true(any(skew_mixed$skew > 0) || any(skew_mixed$skew < 0), info = "CheckDataHealth detects skewness in anomaly scores")

