#pragma once

#include "Prior.h"
#include "MultivariateMarginal.h"
#include "UnivariateMarginal.h"

namespace bcm3 {

	class PriorIndependence : public Prior
	{
	public:
		PriorIndependence(std::shared_ptr<VariableSet> varset);
		virtual ~PriorIndependence();
	
		virtual bool LoadFromXML(const boost::property_tree::ptree& xml_node);
		bool SetMarginal(const std::string& name, std::unique_ptr<UnivariateMarginal>& marginal);

		virtual bool EvaluateLogPDF(size_t threadix, const VectorReal& values, Real& logp) const;
		virtual bool EvaluateMarginalMean(size_t i, Real& mean) const;
		virtual bool EvaluateMarginalVariance(size_t i, Real& var) const;
		virtual bool Sample(VectorReal& values, RNG* rng) const;

		virtual size_t GetNumDirichletDistributions() const;
		virtual void GetDirichletMembership(ptrdiff_t ix, ptrdiff_t& residual_ix, std::vector<ptrdiff_t>& var_ixs) const;

	private:
		std::vector< std::unique_ptr<UnivariateMarginal> > Variables;
		std::vector< std::unique_ptr<MultivariateMarginal> > MultivariateDistributions;
		std::vector<size_t> MultivariateMembership;
		std::vector<size_t> MultivariateVarIx;
	};

}
