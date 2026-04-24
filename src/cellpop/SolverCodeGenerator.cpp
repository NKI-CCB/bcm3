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
		// All of the numerical functions assume that parameters k and n are chosen such that real_epsilon < k^n < real_max/2
		// Then x^n + k^n == 0 is due to x^n being 0, and x^n + k^n == inf is due to x^n being too large
		// Also assume that KM's are such that KM^2 >= minimum
		f << "inline OdeReal hill_function(OdeReal x, OdeReal k, OdeReal n)\n";
		f << "{\n";
		f << "\tif (x <= 0.0) return 0.0;\n";
		f << "\tOdeReal xn = pow(x, n);\n";
		f << "\tOdeReal kn = pow(k, n);\n";
		f << "\tOdeReal xnpkn = xn + kn;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\tif (xnpkn > 3e38f) return 1.0;\n";
		f << "\treturn xn / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn2(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x <= 0.0) return 0.0;\n";
		f << "\tOdeReal x2 = x * x;\n";
		f << "\tOdeReal k2 = k * k;\n";
		f << "\tOdeReal xnpkn = x2 + k2;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\tif (xnpkn > 3e38f) return 10.0;\n";
		f << "\treturn x2 / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn4(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x <= 0.0) return 0.0;\n";
		f << "\tOdeReal x2 = x * x;\n";
		f << "\tOdeReal x4 = x2 * x2;\n";
		f << "\tOdeReal k2 = k * k;\n";
		f << "\tOdeReal k4 = k2 * k2;\n";
		f << "\tOdeReal xnpkn = x4 + k4;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\tif (xnpkn > 3e38f) return 1.0;\n";
		f << "\treturn x4 / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn10(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x <= 0.0) return 0.0;\n";
		f << "\tOdeReal x2 = x * x;\n";
		f << "\tOdeReal x4 = x2 * x2;\n";
		f << "\tOdeReal x8 = x4 * x4;\n";
		f << "\tOdeReal x10 = x2 * x8;\n";
		f << "\tOdeReal k2 = k * k;\n";
		f << "\tOdeReal k4 = k2 * k2;\n";
		f << "\tOdeReal k8 = k4 * k4;\n";
		f << "\tOdeReal k10 = k2 * k8;\n";
		f << "\tOdeReal xnpkn = x10 + k10;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\tif (xnpkn > 3e38f) return 1.0;\n";
		f << "\treturn x10 / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn16(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x <= 0.0) return 0.0;\n";
		f << "\tOdeReal x2 = x * x;\n";
		f << "\tOdeReal x4 = x2 * x2;\n";
		f << "\tOdeReal x8 = x4 * x4;\n";
		f << "\tOdeReal x16 = x8 * x8;\n";
		f << "\tOdeReal k2 = k * k;\n";
		f << "\tOdeReal k4 = k2 * k2;\n";
		f << "\tOdeReal k8 = k4 * k4;\n";
		f << "\tOdeReal k16 = k8 * k8;\n";
		f << "\tOdeReal xnpkn = x16 + k16;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\tif (xnpkn > 3e38f) 1.0;\n";
		f << "\treturn x16 / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_fixedn100(OdeReal x, OdeReal k)\n";
		f << "{\n";
		f << "\tif (x <= 0.0) return 0.0;\n";
		f << "\tOdeReal x2 = x * x;\n";
		f << "\tOdeReal x4 = x2 * x2;\n";
		f << "\tOdeReal x8 = x4 * x4;\n";
		f << "\tOdeReal x16 = x8 * x8;\n";
		f << "\tOdeReal x32 = x16 * x16;\n";
		f << "\tOdeReal x64 = x32 * x32;\n";
		f << "\tOdeReal x100 = x64 * x32 * x4;\n";
		f << "\tOdeReal k2 = k * k;\n";
		f << "\tOdeReal k4 = k2 * k2;\n";
		f << "\tOdeReal k8 = k4 * k4;\n";
		f << "\tOdeReal k16 = k8 * k8;\n";
		f << "\tOdeReal k32 = k16 * k16;\n";
		f << "\tOdeReal k64 = k32 * k32;\n";
		f << "\tOdeReal k100 = k64 * k32 * k4;\n";
		f << "\tOdeReal xnpkn = x100 + k100;\n";
		f << "\tif (xnpkn < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\tif (xnpkn > 3e38f) return 1.0;\n";
		f << "\treturn x100 / xnpkn;\n";
		f << "}\n";
		f << "inline OdeReal hill_function_derivative(OdeReal x, OdeReal k, OdeReal n)\n";
		f << "{\n";
		f << "\tif (x <= 0.0) return 0.0;\n";
		f << "\tOdeReal xn = pow(x, n);\n";
		f << "\tOdeReal kn = pow(k, n);\n";
		f << "\tOdeReal denom = (square(xn + kn));\n";
		f << "\tif (denom < std::numeric_limits<OdeReal>::min()) return 0.0;\n";
		f << "\tif (denom > 3e38f) return 0.0;\n";
		f << "\tOdeReal xnm1 = pow(x, n - 1);\n";
		f << "\treturn kn * n * xnm1 / denom;\n";
		f << "}\n";
		f << "inline OdeReal michaelis_menten_function(OdeReal kcat, OdeReal KM, OdeReal e, OdeReal s)\n";
		f << "{\n";
		f << "\tif (e <= 0) return 0.0;\n";
		f << "\tif (s + KM < 0.1 * KM) {\n";
		f << "\t\tOdeReal bound = -KM + 0.1 * KM;\n";
		f << "\t\tOdeReal offset = (e * kcat * bound / (0.01 * KM) - e * kcat * bound / (KM + bound));\n";
		f << "\t\treturn e * kcat * s / (0.01 * KM) - offset;\n";
		f << "\t}\n";
		f << "\treturn kcat * e * s / (KM + s);\n";
		f << "}\n";
		f << "inline OdeReal michaelis_menten_derivative_enzyme(OdeReal kcat, OdeReal KM, OdeReal e, OdeReal s)\n";
		f << "{\n";
		f << "\tif (e <= 0) return 0.0;\n";
		f << "\tif (s + KM < 0.1 * KM) {\n";
		f << "\t\tOdeReal bound = -KM + 0.1 * KM;\n";
		f << "\t\tOdeReal offset = (kcat * bound / (0.01 * KM) - kcat * bound / (KM + bound));\n";
		f << "\t\treturn kcat * s / (0.01 * KM) - offset;\n";
		f << "\t}\n";
		f << "\treturn kcat * s / (KM + s);\n";
		f << "}\n";
		f << "inline OdeReal michaelis_menten_derivative_substrate(OdeReal kcat, OdeReal KM, OdeReal e, OdeReal s)\n";
		f << "{\n";
		f << "\tif (e <= 0) return 0.0;\n";
		f << "\tif (s + KM <= 0.1 * KM) return e * kcat / (0.01 * KM);\n";
		f << "\treturn e * kcat * KM / (square(KM + s));\n";
		f << "}\n";
		f << "inline OdeReal safepow(OdeReal x, OdeReal n)\n";
		f << "{\n";
		f << "\tif (x <= 0) {\n";
		f << "\t\treturn 0.0;\n";
		f << "\t} else {\n";
		f << "\t\treturn pow(x, n);\n";
		f << "\t}\n";
		f << "}\n";
		f << "inline OdeReal synthcap(OdeReal x)\n";
		f << "{\n";
		f << "\tif (x <= 0) {\n";
		f << "\t\treturn 1.0;\n";
		f << "\t} else {\n";
		f << "\t\tOdeReal x2 = x * x;\n";
		f << "\t\tOdeReal x4 = x2 * x2;\n";
		f << "\t\tOdeReal x8 = x4 * x4;\n";
		f << "\t\treturn 1.0 - x8 * x2;\n";
		f << "\t}\n";
		f << "}\n";
		f << "inline OdeReal synthcap_derivative(OdeReal x, OdeReal dx)\n";
		f << "{\n";
		f << "\tif (x <= 0) {\n";
		f << "\t\treturn 0.0;\n";
		f << "\t} else {\n";
		f << "\t\tOdeReal x2 = x * x;\n";
		f << "\t\tOdeReal x4 = x2 * x2;\n";
		f << "\t\tOdeReal x8 = x4 * x4;\n";
		f << "\t\treturn -10.0 * x8 * x * dx;\n";
		f << "\t}\n";
		f << "}\n";
		f << "inline OdeReal tQSSA(OdeReal k, OdeReal km, OdeReal e, OdeReal s)\n";
		f << "{\n";
		f << "\tOdeReal ekms = e + km + s;\n";
		f << "\treturn 0.5 * k * (ekms - sqrt(ekms * ekms - 4 * e * s));\n";
		f << "}\n";
		f << "inline OdeReal tQSSA_derivative_enzyme(OdeReal k, OdeReal km, OdeReal e, OdeReal s)\n";
		f << "{\n";
		f << "\tOdeReal ekms = e + km + s;\n";
		f << "\treturn k * (0.5 - 0.5 * (km - s + e) / (sqrt(ekms * ekms - 4 * e * s)));\n";
		f << "}\n";
		f << "inline OdeReal tQSSA_derivative_substrate(OdeReal k, OdeReal km, OdeReal e, OdeReal s)\n";
		f << "{\n";
		f << "\tOdeReal ekms = e + km + s;\n";
		f << "\treturn k * (0.5 - 0.5 * (km + s - e) / (sqrt(ekms * ekms - 4 * e * s)));\n";
		f << "}\n";
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
		f << "	generated_derivative\n";
		f << "	generated_jacobian\n";
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
		f << "	set(CMAKE_CXX_FLAGS \"/arch:AVX2\")\n";
		f << "endif(CMAKE_HOST_WIN32)\n";
		f << "if(CMAKE_HOST_UNIX)\n";
		f << "	set(CMAKE_BUILD_TYPE Release)\n";
		f << "	set(CMAKE_CXX_FLAGS \"-O3 -std=c++11 -march=native\")\n";
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
		// HACKY - assume only one of these is installed and that it's the same one used to compile BCM...
		// Look for vcvars64.bat in the common places
		std::string vcvarsfn = "C:\\Program Files (x86)\\Microsoft Visual Studio\\2017\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat";
		FILE* file = fopen(vcvarsfn.c_str(), "r");
		if (file) {
			fclose(file);
#if BUILD_DEBUG
			int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 15 2017 Win64\" . & msbuild generated_derivatives.sln /p:Configuration=Debug").c_str());
#else
			int result = system(("cd " + output_folder + " & \"" + vcvarsfn + "\" & cmake -G \"Visual Studio 15 2017 Win64\" . & msbuild generated_derivatives.sln /p:Configuration=Release").c_str());
#endif
		} else {
			vcvarsfn = "C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat";
			FILE* file = fopen(vcvarsfn.c_str(), "r");
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
	} catch (const std::exception& e) {
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
