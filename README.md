# RangeShiftR <img src="RangeShiftR/man/figures/RSRlogo.png" align="right" height = 150/>

The RangeShiftR package implements the RangeShifter simulation platform for R.

[RangeShifter](https://rangeshifter.github.io/)
is a state-of-the-art eco-evolutionary modelling platform that is becoming 
increasingly used worldwide for both theoretical and applied purposes
[(Bocedi et al. 2014)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12162).

RangeShifter is a spatially-explicit, individual-based simulation platform that 
allows modelling species’ range dynamics, such as expansion and shifting, and 
patch connectivity by linking complex local population dynamics and dispersal 
behaviour, while also taking into account inter-individual variability and 
evolutionary processes. RangeShifter is highly flexible in terms of the spatial 
resolution and extent, and regarding the complexity of the considered ecological 
processes. Due to its modular structure, the level of detail in demographic and 
dispersal processes can be easily adapted to different research questions and 
available data.


## Installation

RangeShiftR is only available from this github repository.
(It may move to CRAN in the future.)

RangeShiftR has to be built from source and requires the package `Rcpp` as
well as a functional C++ compiler toolchain.

```r
# Install RangeShiftR from GitHub:
devtools::install_github("RangeShifter/RangeShiftR-pkg", ref="main", subdir="RangeShiftR")
```

## Usage and help

Please refer to our [website](https://rangeshifter.github.io/) for more information about RangeShifter simulation 
platform. RangeShifter is accompanied by extensive documentation. 

For getting acquainted with the software, we recommend to first read the [manual](https://raw.githubusercontent.com/RangeShifter/RangeShifter-software-and-documentation/master/RangeShifter_v2.0_UserManual.pdf) to understand the conceptual underpinnings of RangeShifter.

Analogous to the RangeShifter GUI, we provide [tutorials](https://rangeshifter.github.io/RangeshiftR-tutorials/) to learn the different features of RangeshiftR using example applications from Bocedi et al. (2014, 2021) and Malchow et al. (2021). These cover some of the main features of RangeShifter, and help becoming familiar with the software.

If you have any further question related to the general concepts and usage of RangeShifter, please browse earlier topics in the [forum pages](https://github.com/RangeShifter/RangeshiftR-tutorials/discussions) or add a new one. Often it is also helpful to review [published studies](https://rangeshifter.github.io/site/references/) using the RangeShifter modelling platform.

For technical questions related to the RangeShiftR package interface and which cannot be answered with the documentation provided above, please browse the [issues section](https://github.com/RangeShifter/RangeShiftR-package-dev/issues) of this repository and open a new issues if required. We also offer *technical support* for the RangeShiftR package via mail (rangeshiftr@uni-potsdam.de) if you follow the guidelines of how to ask for help, e.g. guidelines given by [StackOverflow](https://stackoverflow.com/help/how-to-ask):

## Contributing

See [Contributing guidelines](https://github.com/RangeShifter/RangeShiftR-package-dev/CONTRIBUTING.md)

## See also

-    [Compiled software and documentation](https://github.com/RangeShifter/RangeShifter-software-and-documentation)
-    [RScore](https://github.com/RangeShifter/RScore), source for RangeShifter's core code
-    [RangeShifter batch mode](https://github.com/RangeShifter/RangeShifter_batch_dev), source for the batch mode interface


## References

 - Bocedi G, Palmer SCF, Pe’er G, Heikkinen RK, Matsinos YG, Watts K, Travis JMJ (2014). 
 *RangeShifter: A Platform for Modelling Spatial Eco-Evolutionary Dynamics and 
 Species’ Responses to Environmental Changes.* Methods in Ecology and Evolution 5: 388–96. 

 - Bocedi G, Palmer SCF, Malchow AK, Zurell D, Watts K, Travis JMJ (2021) RangeShifter 2.0: An extended and enhanced platform for modelling spatial eco-evolutionary dynamics and species’ responses to environmental changes. Ecography 44:1453-1462.
 
 - Malchow AK, Bocedi G, Palmer SCF, Travis JMJ, Zurell D (2021) RangeShiftR: an R package for individual-based simulation of spatial eco-evolutionary dynamics and species’ responses to environmental change. Ecography 4: 1443-1452.
