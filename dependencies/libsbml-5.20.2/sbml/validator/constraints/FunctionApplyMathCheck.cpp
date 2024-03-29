/**
 * @cond doxygenLibsbmlInternal
 *
 * @file    FunctionApplyMathCheck.cpp
 * @brief   Ensures <ci> after apply refers to function definition.
 * @author  Sarah Keating
 * 
 * <!--------------------------------------------------------------------------
 * This file is part of libSBML.  Please visit http://sbml.org for more
 * information about SBML, and the latest version of libSBML.
 *
 * Copyright (C) 2020 jointly by the following organizations:
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. University of Heidelberg, Heidelberg, Germany
 *     3. University College London, London, UK
 *
 * Copyright (C) 2019 jointly by the following organizations:
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. University of Heidelberg, Heidelberg, Germany
 *
 * Copyright (C) 2013-2018 jointly by the following organizations:
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. EMBL European Bioinformatics Institute (EMBL-EBI), Hinxton, UK
 *     3. University of Heidelberg, Heidelberg, Germany
 *
 * Copyright (C) 2009-2013 jointly by the following organizations: 
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. EMBL European Bioinformatics Institute (EMBL-EBI), Hinxton, UK
 *  
 * Copyright (C) 2006-2008 by the California Institute of Technology,
 *     Pasadena, CA, USA 
 *  
 * Copyright (C) 2002-2005 jointly by the following organizations: 
 *     1. California Institute of Technology, Pasadena, CA, USA
 *     2. Japan Science and Technology Agency, Japan
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation.  A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution
 * and also available online as http://sbml.org/software/libsbml/license.html
 * ---------------------------------------------------------------------- -->*/

#include <sbml/Model.h>
#include <sbml/Compartment.h>
#include <sbml/Species.h>
#include <sbml/Parameter.h>
#include <sbml/UnitDefinition.h>
#include <sbml/Event.h>
#include <sbml/Reaction.h>
#include <sbml/EventAssignment.h>
#include <sbml/SpeciesReference.h>
#include <sbml/Rule.h>
#include <sbml/math/FormulaFormatter.h>

#include <sbml/units/UnitFormulaFormatter.h>

#include "FunctionApplyMathCheck.h"

/** @cond doxygenIgnored */
using namespace std;
/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN
#ifdef __cplusplus

static const char* PREAMBLE =
    "Outside of a <functionDefinition>, if a 'ci' element is the first "
    "element within a MathML 'apply', then the 'ci''s value can only be "
    "chosen from the set of identifiers of <functionDefinition>s defined in "
    "the SBML model. (References: L2V2 Section 4.3.2.)";


/*
 * Creates a new Constraint with the given @p id.
 */
FunctionApplyMathCheck::FunctionApplyMathCheck (unsigned int id, Validator& v) : MathMLBase(id, v)
{
}


/*
 * Destroys this Constraint.
 */
FunctionApplyMathCheck::~FunctionApplyMathCheck ()
{
}


/*
 * @return the preamble to use when logging constraint violations.
 */
const char*
FunctionApplyMathCheck::getPreamble ()
{
  return PREAMBLE;
}


/*
  * Checks the MathML of the ASTnode 
  * is appropriate for the function being performed
  *
  * If an inconsistency is found, an error message is logged.
  */
void
FunctionApplyMathCheck::checkMath (const Model& m, const ASTNode& node, const SBase & sb)
{

  ASTNodeType_t type = node.getType();

  switch (type) 
  {
    case AST_FUNCTION:

      checkExists(m, node, sb);
      break;

    default:

      checkChildren(m, node, sb);
      break;
  }
}

void
FunctionApplyMathCheck::checkFunctionDefinition(const Model& m, const ASTNode& node,
  const SBase & sb)
{
  // a missing fd is caught by 20301 for l2v3 and lower
  if (m.getLevel() > 2 || (m.getLevel() == 2 && m.getVersion() > 3))
  {
    checkMath(m, node, sb);
  }
}
  
/*
  * Checks that the functionDefinition referred to by a &lt;ci&gt; element exists
  *
  * If <ci> does not refer to functionDefinition id, an error message is logged.
  */
void 
FunctionApplyMathCheck::checkExists (const Model& m, const ASTNode& node, 
                                                const SBase & sb)
{
  std::string name = node.getName();

  if (!m.getFunctionDefinition(name))
    logMathConflict(node, sb);
}


/*
 * @return the error message to use when logging constraint violations.
 * This method is called by logFailure.
 *
 * Returns a message that the given @p id and its corresponding object are
 * in  conflict with an object previously defined.
 */
const string
FunctionApplyMathCheck::getMessage (const ASTNode& node, const SBase& object)
{

  ostringstream oss_msg;

  //oss_msg << getPreamble();
  char * formula = SBML_formulaToString(&node);
  oss_msg << "The formula '" << formula;
  oss_msg << "' in the " << getFieldname() << " element of the <" << object.getElementName();
  oss_msg << "> ";
  switch(object.getTypeCode()) {
  case SBML_INITIAL_ASSIGNMENT:
  case SBML_EVENT_ASSIGNMENT:
  case SBML_ASSIGNMENT_RULE:
  case SBML_RATE_RULE:
    //LS DEBUG:  could use other attribute values, or 'isSetActualId'.
    break;
  default:
    if (object.isSetId()) {
      oss_msg << "with id '" << object.getId() << "' ";
    }
    break;
  }
  oss_msg << "uses '" << node.getName() << "' which is not a function definition id.";
  safe_free(formula);

  return oss_msg.str();
}

#endif /* __cplusplus */

LIBSBML_CPP_NAMESPACE_END
/** @endcond */
