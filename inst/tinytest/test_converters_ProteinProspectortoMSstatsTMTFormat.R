# Test ProteinProspectortoMSstatsTMTFormat ---------------------------
input = system.file("tinytest/raw_data/ProteinProspector/Prospector_TotalTMT.txt",
    package = "MSstatsConvert")
input = data.table::fread(input)
annot = system.file("tinytest/raw_data/ProteinProspector/Annotation.csv",
                                package = "MSstatsConvert")
annot = data.table::fread(annot)
output = ProteinProspectortoMSstatsTMTFormat(input, annot)
expect_equal(ncol(output), 11)
expect_equal(nrow(output), 528)
expect_true("Run" %in% colnames(output))
expect_true("ProteinName" %in% colnames(output))
expect_true("PeptideSequence" %in% colnames(output))
expect_true("Charge" %in% colnames(output))
expect_true("Intensity" %in% colnames(output))
expect_true("TechRepMixture" %in% colnames(output))
expect_true("PSM" %in% colnames(output))
expect_true("Mixture" %in% colnames(output))
expect_true("Condition" %in% colnames(output))
expect_true("BioReplicate" %in% colnames(output))
expect_true("Channel" %in% colnames(output))

# Test ProteinProspectortoMSstatsTMTFormat with missing value ------------
zero_value_entry = 
    output[output$PeptideSequence == "DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR" 
       & output$BioReplicate == "S1",]$Intensity
expect_true(is.na(zero_value_entry))