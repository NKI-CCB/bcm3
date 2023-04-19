# BCM3 - Bayesian inference with Computational Models
BCM3 is a C++/R software package for Bayesian inference with computational models using MCMC sampling.

The focus of BCM3 is on inference for models that require extensive computation to evaluate the likelihood, such as simulations or ODE integrations. Likelihood function are therefore implemented in C++. BCM3 is fully multithreading; supporting parallelization over tempered sampling chains as well as parallelization inside each likelihood evaluation (if applicable).

BCM3 uses an advanced sampling algorithm that can efficiently sample from multimodal distributions or distributions containing strong ridges or parameter correlations. The original BCM implemented a variety of sampling algorithms; in this new version, BCM3 includes only parallel-tempered MCMC with adaptive proposals; as we found that one to generally work most efficiently and has been continuously optimized.

BCM3 includes several likelihoods:
- cellpop; a model describing populations of heterogeneous biological cells
- fISA (feedback-Inference of signaling activity); a model of steady state signaling activity in biological cells (discontinued)
- PopPK; a population pharmacokinetic model
- Several theoretical test/example likelihoods: a banana distribution and a bimodal circular distribution
- You can add your own likelihood as a C++ class

## Installing BCM3

#### Platforms
BCM3 should work on Windows, Linux and Mac; with Visual Studio, GCC or clang. For Visual Studio the community edition is sufficient.

#### Dependencies:
The following dependencies are required:
- [cmake](https://cmake.org/)
- [Boost C++ libraries](https://www.boost.org/)
- [NetCDF](https://www.unidata.ucar.edu/software/netcdf/)
- [hdf5](https://www.hdfgroup.org/solutions/hdf5/)
- [libsbml](https://synonym.caltech.edu/software/libsbml/) (You can remove this dependency by disabling the modules 'sbml', 'cellpop' and 'fISA')
- ([Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and [SUNDIALS](https://computing.llnl.gov/projects/sundials) are used, but versions of these libraries are included as part of BCM3 and do not need to be installed separately)

On Linux, these dependencies can be installed through the package manager; you will need the -dev packages for the C++ library (libboost-dev, libsbml-dev, libnetcdf-dev, libhdf5-dev).
On Windows, pre-compiled versions may work, or you may need to build the dependencies from source.
On Mac, it may be possible to install the dependencies with Homebrew.

For using the R interface and plotting functions, the following R packages are needed:
- XML, hdf5r, extraDistr, crch, pracma

#### Installation steps
1) Copy "external_dependency_locations_template.txt" to "external_dependency_locations.txt", and modify the paths as necessary. (I find this easier to work with than using cmake variables or cmake's FindPackage functionality. If you used Homebrew to install the dependencies on Mac, you may also need to add the relevant locations to your LIBRARY_PATH.)
2) Create an environment variable BCM3_ROOT that points to the root folder of BCM3. (For Mac users, the environment variable might not be active inside R; in that case you might need to add a call to Sys.setenv() to set the BCM3_ROOT variable within R before calling BCM3's R scripts.)
3) Create a folder "build"
4) Run `cmake ..` within the build folder (from command line or with the CMake GUI)
5) Build the generated makefile/project file (with `make` on Linux/Mac, or by building the solution in Visual Studio. In Visual Studio, be sure to use the Release build for production runs.)

## Using BCM3

#### Running an inference
To start an inference run, use 
```
$BCM3_ROOT\bin\bcminf
```
(or `$BCM3_ROOT\bin\Release\bcminf` when built with Visual Studio)

The working directory should contain three files:
- prior.xml -- specifying the variables which are to be sampled and their prior distributions
- likelihood.xml -- specifying the likelihood
- config.txt -- specifying configuration of the sampling algorith, output directory, etc

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

The C++-defined likelihood can also be accessed from R, for example for generating simulations with modified parameter values easily from within R. This requires an R<->C++ bridge which has to be implemented for the likelihood. Such R<->C++ bridges are provided for cellpop, PopPK and fISA. The R<->C++ bridges are somewhat crude and makeshift, and have little argument checking. They may crash your R session if not used exactly right; please use with caution.

#### Creating a custom likelihood
To use BCM3 with your own computational model/likelihood:
- Make a class that derives from bcm3::Likelihood and implements the function EvaluateLogProbability.
- Add the class to LikelihoodFactory::CreateLikelihood.
- In the likelihood.xml file, you can specify the name of your likelihood, such that the factory will create a class instance. You can add configuration/hyperparameters in the likelihood.xml as needed (see the existing likelihood classes for examples).

## Author information
BCM3 was developed by Bram Thijssen; contact: b.thijssen -at- nki -dot- nl

If you find BCM3 useful, I'd love to hear about it. If you use BCM3 in a scientific publication, a citation would be appreciated:
> Bram Thijssen, Tjeerd M. H. Dijkstra, Tom Heskes & Lodewyk F. A. Wessels. BCM: toolkit for Bayesian analysis of Computational Models using samplers. BMC Systems Biology, 2016. https://doi.org/10.1186/s12918-016-0339-3
