# MSstatsConvert

<!-- badges: start -->
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/MSstatsConvert.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MSstatsConvert/)
[![Codecov test coverage](https://codecov.io/gh/Vitek-Lab/MSstatsConvert/branch/devel/graph/badge.svg)](https://codecov.io/gh/Vitek-Lab/MSstatsConvert/branch/devel)
[![Bioconductor Downloads Rank](https://bioconductor.org/shields/downloads/release/MSstatsConvert.svg)](https://bioconductor.org/packages/stats/bioc/MSstatsConvert/)
[![Years in Bioconductor](https://bioconductor.org/shields/years-in-bioc/MSstatsConvert.svg)](https://bioconductor.org/packages/release/bioc/html/MSstatsConvert.html#since)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->

MSstatsConvert is an R/Bioconductor package that imports reports from a wide
range of mass spectrometry signal-processing tools into the standardized format
used by [MSstats](https://github.com/Vitek-Lab/MSstats) and
[MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT) for statistical analysis.
It provides the shared converter functions — one row per feature, run, and
condition — used across the MSstats ecosystem, along with a common three-step
import pipeline (`MSstatsImport()`, `MSstatsClean()`, `MSstatsPreprocess()`)
that reads, cleans, and reformats each tool's output. Because MSstats and its
companion packages call these converters internally, most users interact with
MSstatsConvert indirectly, but the functions can also be used directly for
custom or scripted workflows.

MSstatsConvert is part of the [MSstats](https://github.com/Vitek-Lab/MSstats)
family of packages, developed and maintained by the
[Vitek Lab](https://olga-vitek-lab.khoury.northeastern.edu/) at Northeastern
University. The package and its documentation are also available at
[msstats.org](http://msstats.org).

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsConvert")
```

The development version can be installed directly from this repository:

```r
remotes::install_github("Vitek-Lab/MSstatsConvert")
```

## Quick Start

```r
library(MSstatsConvert)

# Read your search tool's report, then convert it to MSstats format.
# Example: a Skyline report exported as CSV
raw <- data.table::fread("Skyline_report.csv")
msstats_input <- SkylinetoMSstatsFormat(raw)

head(msstats_input)
```

The converted data frame is then passed to `MSstats::dataProcess()` (or
`MSstatsTMT::proteinSummarization()` for the `*toMSstatsTMTFormat()` converters).
See the data-format vignette for the columns each converter produces and the
options it accepts.

## Supported Converters

**Label-free / label-based (for use with [MSstats](https://github.com/Vitek-Lab/MSstats)):**

| Search tool / format | Converter function |
| --- | --- |
| Skyline | `SkylinetoMSstatsFormat()` |
| MaxQuant | `MaxQtoMSstatsFormat()` |
| Progenesis | `ProgenesistoMSstatsFormat()` |
| Spectronaut | `SpectronauttoMSstatsFormat()` |
| Proteome Discoverer | `PDtoMSstatsFormat()` |
| DIA-NN | `DIANNtoMSstatsFormat()` |
| DIA-Umpire | `DIAUmpiretoMSstatsFormat()` |
| FragPipe | `FragPipetoMSstatsFormat()` |
| OpenMS | `OpenMStoMSstatsFormat()` |
| OpenSWATH | `OpenSWATHtoMSstatsFormat()` |
| Metamorpheus | `MetamorpheusToMSstatsFormat()` |
| MZmine | `MZMinetoMSstatsFormat()` |

**Isobaric labeling (for use with [MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT)):**

| Search tool / format | Converter function |
| --- | --- |
| MaxQuant | `MaxQtoMSstatsTMTFormat()` |
| OpenMS | `OpenMStoMSstatsTMTFormat()` |
| Proteome Discoverer | `PDtoMSstatsTMTFormat()` |
| Philosopher / FragPipe | `PhilosophertoMSstatsTMTFormat()` |
| Protein Prospector | `ProteinProspectortoMSstatsTMTFormat()` |
| SpectroMine | `SpectroMinetoMSstatsTMTFormat()` |

## Documentation

- [MSstats data format vignette](vignettes/msstats_data_format.Rmd) — the required MSstats input format and how converters produce it
- [Official website: msstats.org](https://msstats.org)
- [Bioconductor package page and reference manual](https://bioconductor.org/packages/MSstatsConvert)

## Getting Help / Reporting Bugs

- **Questions about usage, statistical methods, or troubleshooting:** please
  post to the [MSstats Google Group](https://groups.google.com/forum/#!forum/msstats).
  This is monitored by the development team and searchable, so it's the fastest
  way to get help and to see if your question has already been answered.
- **Bug reports and feature requests for this repository:** please open a
  [GitHub issue](https://github.com/Vitek-Lab/MSstatsConvert/issues).

## References

MSstatsConvert provides the converters used by the MSstats workflow. If you use
it, please cite the MSstats publication(s):

1. Kohler D, Staniak M, Tsai TH, Huang T, Shulman N, Bernhardt OM, MacLean BX,
   Nesvizhskii AI, Reiter L, Sabido E, Choi M, Vitek O. **MSstats Version 4.0:
   Statistical Analyses of Quantitative Mass Spectrometry-Based Proteomic
   Experiments with Chromatography-Based Quantification at Scale.**
   *J Proteome Res*. 2023;22(5):1466-1482.
   [DOI: 10.1021/acs.jproteome.2c00834](https://doi.org/10.1021/acs.jproteome.2c00834)
2. Choi M, Chang CY, Clough T, Broudy D, Killeen T, MacLean B, Vitek O.
   **MSstats: an R package for statistical analysis of quantitative mass
   spectrometry-based proteomic experiments.** *Bioinformatics*.
   2014;30(17):2524-2526.
   [DOI: 10.1093/bioinformatics/btu305](https://doi.org/10.1093/bioinformatics/btu305)

## Funding

MSstats development has been supported by the Chan Zuckerberg Initiative's
[Essential Open Source Software for Science](https://chanzuckerberg.com/eoss/proposals/).

## License

MSstatsConvert is released under the [Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
