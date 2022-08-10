#pragma once

// From libsbml:
class Species;

class SignalingMolecule
{
public:
	SignalingMolecule();
	~SignalingMolecule();

	static std::unique_ptr<SignalingMolecule> Create(const Species* sp);

	const std::string& GetId() const { return id; }
	const std::string& GetName() const { return name; }

private:
	enum EType {
		Protein,
		ActivatingMutation,
		GeneAmplification,
		CompleteLossMutation,
		Drug,
		Phenotype,
		InvalidType,
	};

	std::string id;
	std::string name;
	EType type;

	friend class SignalingModel;
};
