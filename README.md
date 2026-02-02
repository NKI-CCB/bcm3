# BCM3 - Bayesian inference with Computational Models
BCM3 is a C++/R software package for Bayesian inference with computational models using MCMC sampling.

The focus of BCM3 is on inference for models that require extensive computation to evaluate the likelihood, such as simulations or ODE integrations. Likelihood function are therefore implemented in C++. BCM3 is fully multithreaded; supporting parallelization over tempered sampling chains as well as parallelization inside each likelihood evaluation.

BCM3 uses an advanced sampling algorithm that can efficiently sample from multimodal distributions or distributions containing strong ridges or parameter correlations. The original BCM implemented a variety of sampling algorithms; in this new version, BCM3 includes only parallel-tempered MCMC with adaptive proposals, as we found that sampling algorithm to be the most generally efficient. Several proposal schemes are available.

You can important your own likelihood through C++ evaluation code. BCM3 also includes several likelihoods:
- cellpop; a model describing populations of heterogeneous biological cells
- fISA (feedback-Inference of signaling activity); a model of steady state signaling activity in biological cells (discontinued)
- PopPK; a population pharmacokinetic model
- Several theoretical test/example likelihoods: a banana distribution and a bimodal circular distribution
- You can add your own likelihood either as a C++ class or through a shared library/DLL

## Installing BCM3

BCM3 should work on Windows, Linux and Mac; with Visual Studio, GCC or clang. For Visual Studio the community edition is sufficient.

#### 1. Install dependencies:
The following dependencies are required:
- [cmake](https://cmake.org/)
- [Boost C++ libraries](https://www.boost.org/)
- [NetCDF](https://www.unidata.ucar.edu/software/netcdf/)
- [hdf5](https://www.hdfgroup.org/solutions/hdf5/)
- Optional: [libsbml](https://synonym.caltech.edu/software/libsbml/) - a version of libsbml is included to simplify installation, but a system-wide installation will be used if available. You can remove this dependency by disabling the modules 'sbml', 'cellpop' and 'fISA'
- Note: BCM3 uses [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and [SUNDIALS](https://computing.llnl.gov/projects/sundials), but versions of these libraries are included as part of BCM3 and do not need to be installed separately.

On <strong>Linux</strong>, these dependencies can be installed through the package manager; e.g.:
`apt install libboost-all-dev libnetcdf-dev libhdf5-dev libsbml-dev`.  
On <strong>Windows</strong>, vcpkg is useful to install these dependencies. This can be done with a manifest if you wish to install the dependencies specifically for this project, but personally I use vcpkg in classic mode to have a reusable installation of these packages. This can be done by installing vcpkg separately (the one bundled with Visual Studio can only work in manifest mode), and then installing the packages with the stand-alone vcpkg ( "vcpkg.exe install [boost/netcdf-c/libsbml]").  
On <strong>Mac</strong>, it may be easiest to install these dependencies with Homebrew.
For <strong>Conda</strong> users, you can install the dependencies by `conda install conda-forge::libnetcdf`, `conda install conda-forge::hdf5`, `conda-forge::boost`

For using the <strong>R interface and plotting functions</strong>, the following R packages are needed: ```install.packages("XML", "hdf5r", "extraDistr", "crch", "pracma", "sm")```.

#### 2. Installation of BCM3
Do the following steps:
1) Create an environment variable BCM3_ROOT that points to the root folder of BCM3.  
   For <strong>Windows</strong> users, press Windows+break; click advanced system settings; go to environment variables, and add a new entry BCM3_ROOT pointing to the BCM3 folder.  
   For <strong>Linux</strong> users, ```export BCM3_ROOT="/path/to/bcm3"```.  
   For <strong>Mac</strong> users, the environment variable might not be active inside R; in that case you might need to add a call to Sys.setenv() to set the BCM3_ROOT variable within R before calling BCM3's R scripts.

2) Create a folder "build" within the BCM3 folder

3) Generate solution/makefiles with cmake.  
   For <strong>Windows</strong> when using vcpkg, run     `cmake .. -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake`  within the build folder (change the path to vcpkg).  
   For <strong>Linux/mac/other systems</strong>, run     `cmake ..`        within the build folder.

   If you want to include support for cellpop, fISA or SBML, set the cmake options for these modules to ON.  
   On <strong>Windows</strong> you can do this with the CMake GUI
   On <strong>Linux/mac/other systems</strong>, you can do this by adding "-DINCLUDE_CELLPOP=ON", "-DINCLUDE_fISA=ON", "-DINCLUDE_SBML=ON" to the call to cmake.


4) Build the generated makefile/project file.  
   For <strong>Windows</strong>: build the solution in Visual Studio. Be sure to use the Release build for production runs.  
   For <strong>Linux/Mac</strong>: run
   `make` 
   within the build folder.

## Using BCM3

#### Creating a likelihood

Several likelihoods are available in the package that you can use. If you want to use BCM for your own likelihood function, you can incorporate it in one of two ways: either by creating a DLL that exports a likelihood function, or by adding a class to BCM3.

- The class approach:
  - Make a class that derives from bcm3::Likelihood
  - Implement the function EvaluateLogProbability.
  - Let the sampler know whether EvaluateLogProbability is re-entrant or not by overloading the IsReentrant function.
  - Add the class to LikelihoodFactory::CreateLikelihood
  - In the likelihood.xml file, you can specify the name of your likelihood, such that the factory will create a class instance. You can add configuration/hyperparameters in the likelihood.xml as needed (see the existing likelihood classes for examples).

- The DLL approach:
  - Copy the code in examples/dll_likelihood to a new folder
  - Add initialization or data loading code to the "initialize_likelihood" function; this will be called once at the start of sampling.
  - Add likelihood evaluation code to the "evaluate_log_probability" function. The current paramter values are provided, and the resulting log likelihood should be stored into the log_p variable. 
    <strong>Important</strong>: this function should be re-entrant; i.e. this function will be called from multiple threads simultaneously with different parameter values.

#### Running an inference
To start an inference run, use  
<strong>Windows/Visual studio</strong>: `%BCM3_ROOT%\bin\Release\bcminf`  
<strong>Linux/Mac</strong>: `$BCM3_ROOT/bin/bcminf`  

The working directory should contain three files:
- prior.xml -- specifying the variables which are to be sampled and their prior distributions
- likelihood.xml -- specifying the likelihood
- config.txt -- specifying the configuration of the sampling algorith, output directory, etc

Several examples of these files are provided in the examples folder.

Running bcminf will generate a netCDF file containing the sampling output, which can be loaded into R for further analysis.

Run ```bcminf --help``` to see all the configuration options. These options can be specified on the command line and in the config.txt file.

#### Loading sampling results in R
Sampling results can be loaded into R through

```
source(paste(Sys.getenv("BCM3_ROOT"), "/R/load.r", sep=""))
bcm3res <- bcm3.load.results(..)
```

The samples are available in bcm3res$posterior$samples, and other information is in the bcm3res structure. See str(bcm3res).

Various plotting functions are provided in: 
```
source(paste(Sys.getenv("BCM3_ROOT"), "/R/plots_functions.r", sep=""))
```

The C++-defined likelihood can also be accessed from R, for example for generating simulations with modified parameter values easily from within R. This requires an R-to-C++ bridge which has to be implemented for the likelihood. Such R-to-C++ bridges are provided for cellpop, PopPK and fISA. The R-to-C++ bridges are somewhat crude and makeshift, and have little argument checking. They may crash your R session if not used exactly right; please use with caution.

## Author information
BCM3 was developed by Bram Thijssen; contact: b.thijssen -at- nki -dot- nl

If you find BCM3 useful, I'd love to hear about it. If you use BCM3 in a scientific publication, a citation would be appreciated:
> Bram Thijssen, Tjeerd M. H. Dijkstra, Tom Heskes & Lodewyk F. A. Wessels. BCM: toolkit for Bayesian analysis of Computational Models using samplers. BMC Systems Biology, 2016. https://doi.org/10.1186/s12918-016-0339-3
