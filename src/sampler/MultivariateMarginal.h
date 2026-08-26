#pragma once

namespace bcm3 {

	class RNG;

	class MultivariateMarginal
	{
	public:
		MultivariateMarginal();
		~MultivariateMarginal();

		bool CreateDirichlet();
		bool CreateLKJCholesky(Real eta, size_t D);
		bool AddVariable(Real param, size_t& ix);
		bool Initialize();

		bool Sample(Real* values, RNG* rng) const;
		bool EvaluateLogPDF(Real& logp, const Real* values) const;
		bool EvaluateMarginalMean(size_t i, Real& mean) const;
		bool EvaluateMarginalVariance(size_t i, Real& var) const;
		Real GetLowerBound(size_t i) const;
		Real GetUpperBound(size_t i) const;

		inline size_t GetSize() const { return N; }

	private:
		enum EDistribution
		{
			D_Invalid,
			D_Dirichlet,
			D_LKJCholesky,
		};
	
		EDistribution Type;

		size_t N;
		VectorReal alpha;	// Dirichlet
		Real eta;			// LKJCholesky
		size_t D;			// LKJCholesky
		Real log_normalization_constant;
	};

}
