# BCM3 - Bayesian inference for Computational Models
BCM3 is a C++/R software package for Bayesian inference with computational models using sampling.

The original BCM implemented a variety of sampling algorithms, but as we found one sampler to generally work the best, BCM3 includes only a single sampling algorithm (parallel-tempered MCMC).

The focus is on inference for models that require extensive computation; likelihood function therefore need to be implemented in C++.

BCM3 includes several likelihoods:
- cellpop; a model describing populations of heterogeneous biological cells
- fISA (feedback-Inference of signaling activity); a model of steady state signaling activity in biological cells (discontinued)
- PopPK; a population pharmacokinetic model
- several theoretical test/example likelihoods: a banana distribution and a bimodal circular distribution

## Installing BCM3

#### Platforms
BCM3 is actively used on both Linux and Windows; and should work with GCC and Visual Studio (community edition or other). It may also work on MacOS.

#### Dependencies:
The following dependencies are requires:
- cmake
- Boost C++ libraries
- NetCDF
- hdf5
- libsbml (It should be straightforward to remove this dependency if you don't need it, by disabling the modules 'sbml', 'cellpop' and 'fISA')
- (Eigen and SUNDIALS are used, but versions of these libraries are included as part of BCM3 and do not need to be installed separately)

#### Installation steps
1) Copy "external_dependency_locations_template.txt" to "external_dependency_locations.txt", and modify the paths as necessary. (I find this easier to work with than cmake variables)
2) Create an environment variable BCM3_ROOT that points to the root folder of BCM3
3) Create a folder "build"
4) Run `cmake ..` within the build folder (from command line or with the CMake GUI)
5) Build the generated makefile/project file

## Using BCM3

#### Running an inference
To start an inference run, use 
```
$BCM3_ROOT\bin\bcminf
```
(or `$BCM3_ROOT\bin\Release\bcminf` when built with Visual Studio)

The working directory should contain three files:
prior.xml -- specifying the variables which are to be sampled
likelihood.xml -- specifying the likelihood
config.txt -- specifying configuration of the sampling algorith, output directory, etc

Several examples are provided in the examples folder.

Running bcminf will generate a netCDF file containing the sampling output, which can be loaded into R for further analysis.

#### Loading sampling results in R
Sampling results can be loaded into R through

```
source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
bcm3res <- bcm3.load.results(..)
```

and various plotting functions are provided in 
```
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
```

The C++-defined likelihood can also be accessed from R, for example for generating simulations with modified parameter values easily from within R. This requires an R<->C++ bridge to be implemented for the likelihood function; these bridges are provided for cellpop, PopPK and fISA.
The R<->C++ bridges are somewhat crude and makeshift, and have little argument checking. They may crash your R session if not used exactly right; please use with caution.

#### Creating a custom likelihood
To use BCM3 with your own model/likelihood:
- Make a class that derives from bcm3::Likelihood and implements the function EvaluateLogProbability.
- Add the class to LikelihoodFactory::CreateLikelihood.
- In the likelihood.xml file, you can specify the name of your likelihood, such that the factory will create a class instance. You can add configuration, etc in the likelihood.xml as needed (see the existing likelihood classes for examples).

## Author information
BCM3 was developed by Bram Thijssen; contact: b.thijssen -at- nki -dot- nl
If you find useful BCM3 for anything, I'd love to hear about it. If you use it in a scientific publication, a citation would be appreciated:
Bram Thijssen, Tjeerd M. H. Dijkstra, Tom Heskes & Lodewyk F. A. Wessels. BCM: toolkit for Bayesian analysis of Computational Models using samplers. BMC Systems Biology, 2016. https://doi.org/10.1186/s12918-016-0339-3
