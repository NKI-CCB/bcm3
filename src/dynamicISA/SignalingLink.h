#pragma once

// From libsbml:
class Reaction;
class SignalingModel;

class SignalingLink
{
public:
	SignalingLink();
	~SignalingLink();

	static std::unique_ptr<SignalingLink> Create(const Reaction* re, const SignalingModel* model);
	void NotifySimulatedSpecies(const std::vector<std::string>& simulated_species);

	const std::string& GetFrom() const { return from; }
	const std::string& GetTo() const { return to; }
	const size_t GetFromIx() const { return from_ix; }
	const size_t GetToIx() const { return to_ix; }
	const size_t GetFromSimulatedIx() const { return from_simulated_ix; }
	const size_t GetToSimulatedIx() const { return to_simulated_ix; }
	const std::string GetStrengthName() const { return std::string("strength_") + from + std::string("_") + to; }
	const bool FromIsConstant() const { return from_constant; }
	const bool IsActivating() const { return activating; }

private:
	std::string from;
	std::string to;
	bool activating;

	size_t from_ix;
	size_t to_ix;
	size_t from_simulated_ix;
	size_t to_simulated_ix;
	bool from_constant;
};
