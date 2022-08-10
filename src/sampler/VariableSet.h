#pragma once

namespace bcm3 {

class VariableSet
{
public:
	enum Transform {
		Transform_None,
		Transform_Log,
		Transform_Log10,
		Transform_Logit
	};

	VariableSet();
	~VariableSet();

	bool LoadFromXML(const std::string& filename);
	void AddVariable(const std::string& name, bool logspace = false, bool logistic = false);
	Real TransformVariable(size_t i, Real x) const;

	size_t GetVariableIndex(const std::string& name, bool log_error = true) const;
	inline size_t GetNumVariables() const { return variables.size(); }
	inline const std::string& GetVariableName(size_t i) const { return variables[i]; }
	inline const std::vector<std::string>& GetAllVariableNames() const { return variables; }
	inline const Transform GetVariableTransform(size_t i) const { return transforms[i]; }
	inline const std::string& GetVarsetName() const { return Name; }

private:
	std::string Name;
	std::vector<std::string> variables;
	std::vector<Transform> transforms;
};

}
