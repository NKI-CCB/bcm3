#include "Utils.h"
#include "SignalingLink.h"
#include "SignalingModel.h"

#include <sbml/SBMLTypes.h>
#include <boost/algorithm/string.hpp>

SignalingLink::SignalingLink()
	: activating(true)
	, from_ix(std::numeric_limits<size_t>::max())
	, to_ix(std::numeric_limits<size_t>::max())
	, from_simulated_ix(std::numeric_limits<size_t>::max())
	, to_simulated_ix(std::numeric_limits<size_t>::max())
	, from_constant(false)
{
}

SignalingLink::~SignalingLink()
{

}

std::unique_ptr<SignalingLink> SignalingLink::Create(const Reaction* re, const SignalingModel* model)
{
	std::unique_ptr<SignalingLink> l = std::make_unique<SignalingLink>();

	const XMLNode* annotation = re->getAnnotation();
	for (unsigned int i = 0; i < annotation->getNumChildren(); i++) {
		const XMLNode& child = annotation->getChild(i);
		if (child.getPrefix() == "celldesigner" && child.getName() == "extension") {
			std::string reaction_type_str = child.getChild("reactionType").getChild(0).getCharacters();
			if (reaction_type_str.empty()) {
				LOGERROR("Can't find reaction type for reaction %s - SBML file should be a CellDesigner file", re->getId().c_str());
				l.reset();
				return l;
			}

			if (reaction_type_str == "POSITIVE_INFLUENCE") {
				l->activating = true;
			} else if (reaction_type_str == "NEGATIVE_INFLUENCE") {
				l->activating = false;
			} else {
				LOGERROR("Unrecognized reaction type %s for reaction %s - SBML file should use reduced notation", reaction_type_str.c_str(), re->getId().c_str());
				l.reset();
				return l;
			}
		}
	}

	if (re->getNumReactants() != 1) {
		LOGERROR("Can only handle reactions with 1 reactant");
		l.reset();
		return l;
	}
	if (re->getNumProducts() != 1) {
		LOGERROR("Can only handle reactions with 1 product");
		l.reset();
		return l;
	}

	l->from = model->GetMoleculeNameById(re->getReactant(0)->getSpecies());
	l->to = model->GetMoleculeNameById(re->getProduct(0)->getSpecies());
	l->from_ix = model->GetMoleculeIxByName(l->from);
	l->to_ix = model->GetMoleculeIxByName(l->to);

	return l;
}

void SignalingLink::NotifySimulatedSpecies(const std::vector<std::string>& simulated_species)
{
	auto it = std::find(simulated_species.begin(), simulated_species.end(), from);
	if (it == simulated_species.end()) {
		from_constant = true;
	} else {
		from_constant = false;
		from_simulated_ix = it - simulated_species.begin();
	}

	it = std::find(simulated_species.begin(), simulated_species.end(), to);
	ASSERT(it != simulated_species.end());
	to_simulated_ix = it - simulated_species.begin();
}
