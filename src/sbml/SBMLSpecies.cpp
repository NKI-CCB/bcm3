#include "Utils.h"
#include "SBMLSpecies.h"

SBMLSpecies::SBMLSpecies()
	: type(Invalid)
	, initial_value(std::numeric_limits<Real>::quiet_NaN())
{
}

SBMLSpecies::~SBMLSpecies()
{
}

bool SBMLSpecies::Initialize(const Species* species, const XMLNode* model_annotation)
{
	id = species->getId();
	name = species->getName();
	initial_value = species->getInitialAmount();

	const XMLNode* annotation = species->getAnnotation();
	if (annotation) {
		for (unsigned int i = 0; i < annotation->getNumChildren(); i++) {
			const XMLNode& child = annotation->getChild(i);
			if (child.getPrefix() == "celldesigner" && child.getName() == "extension") {
				std::string spclassstr = child.getChild("speciesIdentity").getChild("class").getChild(0).getCharacters();
				if (spclassstr.empty()) {
					LOGERROR("Can't find species class for species %s [%s] - is the SBML file a CellDesigner file?", id.c_str(), name.c_str());
					return false;
				}

				if (spclassstr == "GENE") {
					type = Gene;
				} else if (spclassstr == "RNA") {
					type = Transcript;
				} else if (spclassstr == "PROTEIN") {
					type = Protein;
				} else if (spclassstr == "COMPLEX") {
					type = Complex;
				} else if (spclassstr == "DEGRADED") {
					type = Sink;
				} else if (spclassstr == "DRUG") {
					type = Drug;
				} else if (spclassstr == "PHENOTYPE") {
					type = Phenotype;
				} else if (spclassstr == "UNKNOWN") {
					type = Unknown;
				} else {
					LOGERROR("Unrecognized species type %s for species %s [%s]", spclassstr.c_str(), id.c_str(), name.c_str());
					return false;
				}

				if (type == Transcript) {
					name += "_mRNA";
				}

				if (type == Protein) {
					// Retrieve list of modifications from the model annotation
					std::string protein_reference = child.getChild("speciesIdentity").getChild("proteinReference").getChild(0).getCharacters();
					if (protein_reference.empty()) {
						LOGERROR("Can't find protein reference for species %s [%s] - is the SBML file a CellDesigner file?", id.c_str(), name.c_str());
						return false;
					}

					for (unsigned int j = 0; j < model_annotation->getNumChildren(); j++) {
						const XMLNode& model_child = model_annotation->getChild(j);
						if (model_child.getPrefix() == "celldesigner" && model_child.getName() == "extension") {
							const XMLNode& protein_list = model_child.getChild("listOfProteins");
							for (unsigned int k = 0; k < protein_list.getNumChildren(); k++) {
								if (protein_list.getChild(k).getAttrValue("id") == protein_reference) {
									const XMLNode& residue_list = protein_list.getChild(k).getChild("listOfModificationResidues");
									for (unsigned int l = 0; l < residue_list.getNumChildren(); l++) {
										residues[residue_list.getChild(l).getAttrValue("id")] = residue_list.getChild(l).getAttrValue("name");
									}
								}
							}
						}
					}

					// See if any of the modifications are in fact modified in this species
					const XMLNode& modification_list = child.getChild("speciesIdentity").getChild("state").getChild("listOfModifications");
					for (unsigned int j = 0; j < modification_list.getNumChildren(); j++) {
						const XMLNode& mod = modification_list.getChild(j);
						residue_modifications[mod.getAttrValue("residue")] = mod.getAttrValue("state");
					}
				}
			}
		}
	} else {
		type = Unknown;
	}
	
	return true;
}

std::string SBMLSpecies::GetFullName() const
{
	switch (type) {
	case Gene:
		return name + std::string("_gene");

	case Transcript:
		return name + std::string("_mrna");
		
	case Protein:
		{
			std::string fullname = name + std::string("_protein");
			for (std::map<std::string, std::string>::const_iterator ri = residues.begin(); ri != residues.end(); ++ri) {
				std::map<std::string, std::string>::const_iterator mi = residue_modifications.find(ri->first);
				if (mi == residue_modifications.end()) {
					fullname += std::string("_") + ri->second + std::string("_empty");
				} else {
					fullname += std::string("_") + ri->second + std::string("_") + mi->second;
				}
			}
			return fullname;
		}
	
	case Complex:
	case Drug:
	case Phenotype:
	case Unknown:
		return name;
	
	case Sink:
		return std::string("sink");

	case Invalid:
	default:
		LOGERROR("Invalid species type %d for species %s [%s]", (int)type, id.c_str(), name.c_str());
		return "";
	};
}
