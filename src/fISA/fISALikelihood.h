#pragma once

#include "Likelihood.h"
#include "TaskManager.h"
#include "VariableSet.h"
#include <boost/multi_array.hpp>

class SignalingNetwork;
class fISAExperiment;

class fISALikelihood : public bcm3::Likelihood
{
public:
	fISALikelihood(size_t numthreads, size_t evaluation_threads);
	~fISALikelihood();
	
	virtual bool Initialize(std::shared_ptr<const bcm3::VariableSet> varset, boost::property_tree::ptree likelihood_node);
	virtual bool IsReentrant() { return false; }
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

	size_t GetNumSignalingMolecules() const;
	const fISAExperiment* GetExperiment(const std::string& experiment) const;

#if 0
	void SetBootstrap(unsigned long rngseed);
	void SetSafeBayesIx(size_t ix);
	void AddSafeBayesCvIx(size_t ix);
#endif

private:
	bool CalculateModel(VectorReal& activities, const Real* values) const;
	
	inline Real GetValue(const Real* values, size_t varix) const {
		ASSERT(varix < VarSet->GetNumVariables());
		return values[varix];
	}

	size_t numthreads;
	size_t evaluation_threads;

	std::shared_ptr<const bcm3::VariableSet> VarSet;
	std::unique_ptr<SignalingNetwork> network;
	std::vector< std::unique_ptr<fISAExperiment> > experiments;
	bcm3::TaskManager cell_line_evaluation_task_manager;
};
