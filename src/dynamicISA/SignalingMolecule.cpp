#include "Utils.h"
#include "SignalingMolecule.h"

#include <sbml/SBMLTypes.h>
#include <boost/algorithm/string.hpp>

SignalingMolecule::SignalingMolecule()
	: type(InvalidType)
{
}

SignalingMolecule::~SignalingMolecule()
{
}

std::unique_ptr<SignalingMolecule> SignalingMolecule::Create(const Species* sp)
{
	std::unique_ptr<SignalingMolecule> s = std::make_unique< SignalingMolecule>();

	const XMLNode* annotation = sp->getAnnotation();
	if (!annotation) {
		s.reset();
		return s;
	}

	std::string notes;
	const XMLNode* notes_node = sp->getNotes();
	if (notes_node) {
		const XMLNode& child = notes_node->getChild(0);
		const XMLNode& child1 = child.getChild(1);
		const XMLNode& child2 = child1.getChild(0);
		notes = child2.getCharacters();
		boost::trim(notes);
	}

	s->id = sp->getId();
	s->name = sp->getName();

	for (unsigned int ai = 0; ai < annotation->getNumChildren(); ai++) {
		const XMLNode& child = annotation->getChild(ai);
		if (child.getPrefix() == "celldesigner" && child.getName() == "extension") {
			std::string species_class = child.getChild("speciesIdentity").getChild("class").getChild(0).getCharacters();
			if (species_class.empty()) {
				LOGERROR("Can't find species class for species %s [%s] - SBML file should be a CellDesigner file", sp->getId().c_str(), sp->getName().c_str());
				s.reset();
				return s;
			}

			if (species_class == "PROTEIN") {
				s->type = SignalingMolecule::Protein;
			} else if (species_class == "GENE") {
				if (notes == "complete_loss") {
					s->type = SignalingMolecule::CompleteLossMutation;
				} else if (notes == "gene_amplification") {
					s->type = SignalingMolecule::GeneAmplification;
				} else {
					s->type = SignalingMolecule::ActivatingMutation;
				}
			} else if (species_class == "DRUG") {
				s->type = SignalingMolecule::Drug;

#if 0
				if (notes == "inhibit activity") {
					s.drug_inhibition_type = SignalingMolecule::InhibitActivity;
				} else if (notes == "inhibit activity,alter susceptibility") {
					s.drug_inhibition_type = SignalingMolecule::InhibitActivityAlterSusceptibility;
				} else if (notes == "inhibit activation") {
					s.drug_inhibition_type = SignalingMolecule::InhibitActivation;
				} else if (notes == "activate") {
					s.drug_inhibition_type = SignalingMolecule::Activate;
				} else if (notes == "") {
					LOGERROR("The drug \"%s\" does not have a note to specify whether the drug inhibits the activity or the activation of its targets. Please add a \"note\" to the species in the SBML file specifying either \"inhibit activity\", \"inhibit activation\" or \"activate\".", s.name.c_str());
					delete document;
					return false;
				} else {
					LOGERROR("The drug \"%s\" has notes, but the notes do not appear to specify whether the drug inhibits the activity or the activation of its targets. Please add a \"note\" to the species in the SBML file specifying either \"inhibit activity\", \"inhibit activation\" or \"activate\".", s.name.c_str());
					delete document;
					return false;
				}
#endif
			} else if (species_class == "PHENOTYPE") {
				s->type = SignalingMolecule::Phenotype;
			} else {
				LOGERROR("Unrecognized species type %s for species %s [%s]", species_class.c_str(), sp->getId().c_str(), sp->getName().c_str());
				s.reset();
				return s;
			}
		}
	}

	return s;
}
