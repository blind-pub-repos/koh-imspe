# koh-imspe
Blind repo for the article "Active Learning for Simulator Calibration"

## Example Code

Selected example code with a reasonable computation time can be found in the `KOH-IMSPE-example-code.Rmd` file with a description.  This code sources `.R` files from both the `code` and `exp` (shorthand for experiments) folder.

## Monte Carlo Experiments

Code for each of the Monte Carlo (MC) experiments run can be found in the `exp` folder.  Code and data specific to each of the three examples can be found in the `exp/1d-sinusoid`, `exp/surrogates-ch8p2`, and `exp/sx` folders.  To start the MC experiment, source the `.R` file containing `exp` in the name.  For some examples, this file contains the code for each design method tested.  For other examples, the `.R` file containing the name `append` must be run subsequently to the `exp` file in order to collect data for all the design types, after setting the `previous.file` and `previous.file.short` variables to contain the name for the previously obtained `.RData` file.

## Publication

In order to compile `active-learning-for-simulator-calibration.Rmd` the `.RData` files are necessary and available upon request.  They must be placed in the sub directory under `exp` for the relevant experiment, eg. the sinusoid data set must be placed in `exp/1d-sinusoid`.