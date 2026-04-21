#include "SolverCodeGenerator.h"
#include "Experiment.h"
#include "Utils.h"
#include "SBMLModel.h"
#include "SBMLSpecies.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <fstream>
#include <cstdlib>

#if !PLATFORM_WINDOWS
#include <dlfcn.h>
#endif

SolverCodeGenerator::SolverCodeGenerator()
    : derivative_dll(NULL)
{
}

SolverCodeGenerator::~SolverCodeGenerator()
{
    if (derivative_dll) {
#if PLATFORM_WINDOWS
        FreeLibrary(derivative_dll);
#else
        dlclose(derivative_dll);
#endif
    }
}

bool SolverCodeGenerator::GenerateAndCompileSolverCode(Experiment* experiment, SBMLModel& cell_model, const std::string& codegen_name)
{
    std::string hash = (boost::format("%x") % std::hash<std::string>{}(codegen_name)).str();
    std::string output_folder = std::string("codegen_") + hash + "/";

    std::string codename_file = output_folder + "codegen_name";
    if (boost::filesystem::exists(codename_file)) {
        std::ifstream infile(codename_file);
        std::string modelfn;
        infile >> modelfn;
        if (modelfn != codegen_name) {
            // Hash mismatch
            LOGERROR("Code generation hashing collision; delete \"%s\" if you wish to regenerate the simulation code", output_folder.c_str());
            return false;
        }
    }

    // Test if the dll is already there.
#if PLATFORM_WINDOWS
#if BUILD_DEBUG
    std::string derivative_dll_fn = output_folder + "Debug/generated_derivatives.dll";
#else
    std::string derivative_dll_fn = output_folder + "Release/generated_derivatives.dll";
#endif
    derivative_dll = LoadLibrary(derivative_dll_fn.c_str());
    if (derivative_dll) {
        derivative = (derivative_fn)GetProcAddress(derivative_dll, "generated_derivative");
        jacobian = (jacobian_fn)GetProcAddress(derivative_dll, "generated_jacobian");
#else
    std::string derivative_dll_fn = output_folder + "libgenerated_derivatives.so";
    derivative_dll = dlopen(derivative_dll_fn.c_str(), RTLD_NOW);
    if (derivative_dll) {
        derivative = (derivative_fn)dlsym(derivative_dll, "generated_derivative");
        jacobian = (jacobian_fn)dlsym(derivative_dll, "generated_jacobian");
#endif
        if (!derivative || !jacobian) {
            LOG("Unable to find generated derivative or jacobian in the dll, recompiling derivative");
            printf("Unable to find generated derivative or jacobian in the dll, recompiling...\n");
        } else {
            LOG("Found DLL with derivative function, reusing it.");
            return true;
        }
    } else {
        LOG("Can't find dll with derivative functions for experiment \"%s\", generating and compiling derivative code.", experiment->GetName().c_str());
        printf("Can't find dll with derivative functions for experiment \"%s\", generating and compiling derivative code...\n", experiment->GetName().c_str());
    }

    std::string bcm_path;
    char* bcm_root_env = getenv("BCM3_ROOT");
    if (!bcm_root_env) {
        LOGERROR("BCM3_ROOT environment variable has not been specified!");
        return false;
    } else {
        bcm_path = bcm_root_env;
        bcm3::fix_path(bcm_path);
    }

    std::string code = cell_model.GenerateCode();

    if (!boost::filesystem::is_directory(output_folder)) {
        ASSERT(*output_folder.rbegin() == '/');
        std::string create_path = output_folder.substr(0, output_folder.size() - 1);
        if (!boost::filesystem::create_directories(create_path)) {
            LOGERROR("Unable to make directory for compiling generated code");
            return false;
        }
    }

    std::ofstream f;
    f.open(output_folder + "code.cpp", std::fstream::out);
    if (f.is_open()) {
#if PLATFORM_WINDOWS
        f << "#define WIN32_LEAN_AND_MEAN\n";
        f << "#define NOMINMAX\n";
        f << "#define EIGEN_NO_IO\n";
        f << "#include <windows.h>\n";
#endif
        f << "#include <cmath>\n";
        f << "#include \"sundials/sundials_matrix.h\"\n";
        f << "#include <Eigen/Dense>\n";
        f << "typedef Eigen::VectorXd VectorReal;\n";
        f << "typedef Eigen::MatrixXd MatrixReal;\n";
        f << "#include \"LinearAlgebraSelector.h\"\n";
        f << std::endl;
        for (size_t i = 0; i < cell_model.GetNumODEIntegratedSpecies(); i++) {
            f << "// species[" + std::to_string(i) + "] -- " + cell_model.GetODEIntegratedSpeciesName(i) + " -- " + cell_model.GetODEIntegratedSpecies(i)->GetFullName() + "\n";
        }
        f << std::endl;
        f << "#define EXPORT_PREFIX  extern \"C\"\n";
        f << std::endl;
        f << "inline OdeReal square(OdeReal x) { return x * x; }\n";
        f << std::endl;
        f << code;
        f.close();
    } else {
        LOGERROR("Unable to open output file for generated code");
        return false;
    }

    // Dynamic library definitions file
    f.open(output_folder + "definitions.def", std::fstream::out);
    if (f.is_open()) {
        f << "EXPORTS\n";
        f << "\tgenerated_derivative\n";
        f << "\tgenerated_jacobian\n";
        f.close();
    } else {
        LOGERROR("Unable to open output file for generated code definitions");
        return false;
    }

    f.open(output_folder + "CMakeLists.txt", std::fstream::out);
    if (f.is_open()) {
        f << "cmake_minimum_required(VERSION 3.16)\n";
        f << "project(generated_derivatives CXX)\n";
        f << "\n";
        f << "if(CMAKE_HOST_WIN32)\n";
        f << "\tset(CMAKE_CXX_FLAGS \\\"/arch:AVX2\\\")\n";
        f << "endif(CMAKE_HOST_WIN32)\n";
        f << "if(CMAKE_HOST_UNIX)\n";
        f << "\tset(CMAKE_BUILD_TYPE Release)\n";
        f << "\tset(CMAKE_CXX_FLAGS \\\"-O3 -std=c++11 -march=native\\\")\n";
        f << "endif(CMAKE_HOST_UNIX)\n";
        f << "\n";
        f << "include_directories(" << bcm_path << "dependencies/eigen-3.4-rc1/ " << bcm_path << "dependencies/cvode-5.3.0/include/ " << bcm_path << "src/odecommon " << bcm_path << "src/utils)\n";
        f << "add_library(generated_derivatives MODULE code.cpp definitions.def)\n";
        f.close();
    } else {
        LOGERROR("Unable to open output file for generated CMake file");
        return false;
    }

    // Compile
    LOG("Compiling generated code...");
    try {
#if PLATFORM_WINDOWS
        std::string vcvarsfn = "C:\\Program Files (x86)\\Microsoft Visual Studio\\2017\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat";
        FILE *file = fopen(vcvarsfn.c_str(), "r");
        if (file) {
            fclose(file);
#if BUILD_DEBUG
            int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 15 2017 Win64\" . & msbuild generated_derivatives.sln /p:Configuration=Debug").c_str());
#else
            int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 15 2017 Win64\" . & msbuild generated_derivatives.sln /p:Configuration=Release").c_str());
#endif
        } else {
            vcvarsfn = "C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat";
            FILE *file = fopen(vcvarsfn.c_str(), "r");
            if (file) {
                fclose(file);
#if BUILD_DEBUG
                int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 16 2019\" . & msbuild generated_derivatives.sln /p:Configuration=Debug").c_str());
#else
                int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 16 2019\" . & msbuild generated_derivatives.sln /p:Configuration=Release").c_str());
#endif
            } else {
                vcvarsfn = "C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat";
                FILE* file = fopen(vcvarsfn.c_str(), "r");
                if (file) {
                    fclose(file);
#if BUILD_DEBUG
                    int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 17 2022\" . & msbuild generated_derivatives.sln /p:Configuration=Debug").c_str());
#else
                    int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 17 2022\" . & msbuild generated_derivatives.sln /p:Configuration=Release").c_str());
#endif
                } else {
                    vcvarsfn = "C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat";
                    FILE* file = fopen(vcvarsfn.c_str(), "r");
                    if (file) {
                        fclose(file);
#if BUILD_DEBUG
                        int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 17 2022\" . & msbuild generated_derivatives.sln /p:Configuration=Debug").c_str());
#else
                        int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 17 2022\" . & msbuild generated_derivatives.sln /p:Configuration=Release").c_str());
#endif
                    } else {
                        LOGERROR("Not sure where to find vcvars64.bat for compiling the code...");
                        return false;
                    }
                }
            }
        }
#else
        int result = system(("cd " + output_folder + " ; cmake . ; make").c_str());
#endif
    } catch (const std::exception &e) {
        LOGERROR("Process error: %s", e.what());
        return false;
    }

    // Load the compiled dynamic library
#if PLATFORM_WINDOWS
    derivative_dll = LoadLibrary(derivative_dll_fn.c_str());
    if (!derivative_dll) {
        LOGERROR("LoadLibrary failed: %u", GetLastError());
        return false;
    }
    derivative = (derivative_fn)GetProcAddress(derivative_dll, "generated_derivative");
    jacobian = (jacobian_fn)GetProcAddress(derivative_dll, "generated_jacobian");
#else
    derivative_dll = dlopen(derivative_dll_fn.c_str(), RTLD_NOW);
    if (!derivative_dll) {
        LOGERROR("dlopen failed");
        return false;
    }
    derivative = (derivative_fn)dlsym(derivative_dll, "generated_derivative");
    jacobian = (jacobian_fn)dlsym(derivative_dll, "generated_jacobian");
#endif

    if (!derivative || !jacobian) {
        LOGERROR("Unable to find generated derivative or jacobian in the dll");
        return false;
    }

    f.open(codename_file, std::fstream::out);
    if (f.is_open()) {
        f << codegen_name;
        f.close();
    } else {
        LOGERROR("Unable to open codegen name file");
        return false;
    }

    return true;
}
