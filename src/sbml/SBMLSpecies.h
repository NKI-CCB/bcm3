#pragma once

#include <sbml/SBMLTypes.h>

class SBMLSpecies
{
public:
	SBMLSpecies();
	~SBMLSpecies();

	enum EType {
		Gene,
		Transcript,
		Protein,
		Complex,
		Drug,
		Phenotype,
		Sink,
		Unknown,

		Invalid
	};

	bool Initialize(const Species* species, const XMLNode* doc_annotation);

	inline const std::string& GetId() const { return id; }
	inline const std::string& GetName() const { return name; }
	inline EType GetType() const { return type; }
	inline Real GetInitialValue() const { return initial_value; }
	
	std::string GetFullName() const;

private:
	std::string id;
	std::string name;
	EType type;
	std::map<std::string, std::string> residues;
	std::map<std::string, std::string> residue_modifications;
	Real initial_value;
};
