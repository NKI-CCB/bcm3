#pragma once

#if 1
	typedef double Real;
	typedef Eigen::VectorXd VectorReal;
	typedef Eigen::MatrixXd MatrixReal;
#else
	typedef float Real;
	typedef Eigen::VectorXf VectorReal;
	typedef Eigen::MatrixXf MatrixReal;
#endif

#if defined(_MSC_VER) || defined(__BORLANDC__)
	typedef unsigned __int64 uint64;
	typedef signed __int64 int64;
#else
	typedef unsigned long long uint64;
	typedef signed long long int64;
#endif
