#pragma once

namespace bcm3 {

class RNG;
class VariableSet;

class Prior
{
public:
	Prior(std::shared_ptr<VariableSet> varset);
	virtual ~Prior();

	static std::unique_ptr<Prior> Create(const std::string& filename, std::shared_ptr<VariableSet> varset, size_t numthreads);
		
	virtual bool LoadFromXML(const boost::property_tree::ptree& xml_node) { return false; }
	virtual bool EvaluateLogPDF(size_t threadix, const VectorReal& values, Real& logp) const = 0;
	virtual bool EvaluateMarginalMean(size_t i, Real& mean) const { return false; }
	virtual bool EvaluateMarginalVariance(size_t i, Real& var) const { return false; }
	virtual bool Sample(VectorReal& values, RNG* rng) const = 0;

	inline Real GetLowerBound(size_t i) const { return LowerBounds(i); }
	inline Real GetUpperBound(size_t i) const { return UpperBounds(i); }

protected:
	std::shared_ptr<VariableSet> varset;
	VectorReal LowerBounds;
	VectorReal UpperBounds;
};

}
