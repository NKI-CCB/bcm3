/**
 * @file    CVTerm.cpp
 * @brief   CVTerm I/O
 * @author  Sarah Keating
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
 * the Free Software Foundation.  A copy of the license agreement is
 * provided in the file named "LICENSE.txt" included with this software
 * distribution.  It is also available online at
 * http://sbml.org/software/libsbml/license.html
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


#include <sbml/common/common.h>
#include <sbml/xml/XMLToken.h>
#include <sbml/xml/XMLNode.h>
#include <sbml/xml/XMLInputStream.h>
#include <sbml/xml/XMLErrorLog.h>

#include <sbml/SBase.h>

#include <sbml/SBMLErrorLog.h>

#include <sbml/util/util.h>
#include <sbml/util/List.h>

#include <sbml/annotation/CVTerm.h>


/** @cond doxygenIgnored */
using namespace std;
/** @endcond */

LIBSBML_CPP_NAMESPACE_BEGIN
#ifdef __cplusplus


/*
 * create a new CVTerm
 */
CVTerm::CVTerm(QualifierType_t type) :
    mHasBeenModified (false)
  , mNestedCVTerms (NULL)
  , mCapturedInStoredAnnotation (false)

{
  mResources = new XMLAttributes();

  mQualifier = UNKNOWN_QUALIFIER;
  mModelQualifier = BQM_UNKNOWN;
  mBiolQualifier = BQB_UNKNOWN;

  setQualifierType(type);

}


/*
 * create a new CVTerm from an XMLNode
 * this assumes that the XMLNode has a prefix 
 * that represents a CV term
 */
CVTerm::CVTerm(const XMLNode node) :
    mHasBeenModified (false)
  , mNestedCVTerms(NULL)
  , mCapturedInStoredAnnotation(false)
{
  const string& name = node.getName();
  const string& prefix = node.getPrefix();
  XMLNode Bag = node.getChild(0);

  mResources = new XMLAttributes();

  mQualifier = UNKNOWN_QUALIFIER;
  mModelQualifier = BQM_UNKNOWN;
  mBiolQualifier = BQB_UNKNOWN;

  if (prefix == "bqbiol")
  {
    setQualifierType(BIOLOGICAL_QUALIFIER);
    setBiologicalQualifierType(name);
  }
  else if (prefix == "bqmodel")
  {
    setQualifierType(MODEL_QUALIFIER);
    setModelQualifierType(name);
  }


  for (unsigned int n = 0; n < Bag.getNumChildren(); n++)
  {
    // need to look to see if we have a nested annotation
    XMLNode child = Bag.getChild(n);

    if (child.getNumChildren() == 0)
    {
      for (int b = 0; b < child.getAttributes().getLength(); b++)
      {
        addResource(child.getAttributes().getValue(b));
      }
    }
    else
    {
      // check we have a correct node
      if (mNestedCVTerms == NULL)
      {
        mNestedCVTerms = new List();
      }
      CVTerm *term = new CVTerm(child);
      mNestedCVTerms->add((void *)(term));
    }
  }

}


/*
 * destructor
 */
CVTerm::~CVTerm()
{
  delete mResources;
  mResources = NULL;

  if (mNestedCVTerms != NULL)
  {
    unsigned int size = mNestedCVTerms->getSize();
    while (size > 0) 
    {
      delete static_cast<CVTerm*>( mNestedCVTerms->remove(0) );
      size--;
    }
    delete mNestedCVTerms;
  }
}

/*
 * Copy constructor; creates a copy of a CVTerm.
 * 
 * @param orig the CVTerm instance to copy.
 */
CVTerm::CVTerm(const CVTerm& orig)
{
  {
    mQualifier      = orig.mQualifier;
    mModelQualifier = orig.mModelQualifier;
    mBiolQualifier  = orig.mBiolQualifier;
    mResources      = new XMLAttributes(*orig.mResources);
    mHasBeenModified = orig.mHasBeenModified;
    mCapturedInStoredAnnotation = orig.mCapturedInStoredAnnotation;

    if(orig.mNestedCVTerms != NULL)
    {
      this->mNestedCVTerms  = new List();
      unsigned int i,iMax = orig.mNestedCVTerms->getSize();
      for(i = 0; i < iMax; ++i)
      {
        this->mNestedCVTerms
          ->add(static_cast<CVTerm*>(orig.mNestedCVTerms->get(i))->clone());
      }
    }
    else
    {
      this->mNestedCVTerms = NULL;
    }
  }
}

/*
 * Assignment operator for CVTerm.
 */
CVTerm& 
CVTerm::operator=(const CVTerm& rhs)
{
  if(&rhs!=this)
  {
    mQualifier       = rhs.mQualifier;
    mModelQualifier  = rhs.mModelQualifier;
    mBiolQualifier   = rhs.mBiolQualifier;

    delete mResources;
    mResources=new XMLAttributes(*rhs.mResources);

    mHasBeenModified = rhs.mHasBeenModified;
    mCapturedInStoredAnnotation = rhs.mCapturedInStoredAnnotation;

    if(this->mNestedCVTerms != NULL)
    {
      unsigned int size = this->mNestedCVTerms->getSize();
      while (size--) delete static_cast<CVTerm*>( this->mNestedCVTerms->remove(0) );
      delete this->mNestedCVTerms;
    }

    if(rhs.mNestedCVTerms != NULL)
    {
      this->mNestedCVTerms  = new List();
      unsigned int i,iMax = rhs.mNestedCVTerms->getSize();
      for(i = 0; i < iMax; ++i)
      {
        this->mNestedCVTerms
          ->add(static_cast<CVTerm*>(rhs.mNestedCVTerms->get(i))->clone());
      }
    }
    else
    {
      this->mNestedCVTerms = NULL;
    }
  }

  return *this;
}





/*
 * clones the CVTerm
 */  
CVTerm* CVTerm::clone() const
{
    CVTerm* term=new CVTerm(*this);
    return term;
}


/*
 * set the qualifier type
 */
int 
CVTerm::setQualifierType(QualifierType_t type)
{
  mQualifier = type;
  
  if (mQualifier == MODEL_QUALIFIER)
  {
    mBiolQualifier = BQB_UNKNOWN;
  }
  else
  {
    mModelQualifier = BQM_UNKNOWN;
  }

  mHasBeenModified = true;
  return LIBSBML_OPERATION_SUCCESS;
}


/*
 * set the model qualifier type
 * this should be consistent with the mQualifier == MODEL_QUALIFIER
 */
int 
CVTerm::setModelQualifierType(ModelQualifierType_t type)
{
  if (mQualifier == MODEL_QUALIFIER)
  {
    mModelQualifier = type;
    mBiolQualifier = BQB_UNKNOWN;
    mHasBeenModified = true;
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    mModelQualifier = BQM_UNKNOWN;
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
}


/*
 * set the biological qualifier type
 * this should be consistent with the mQualifier == BIOLOGICAL_QUALIFIER
 */
int 
CVTerm::setBiologicalQualifierType(BiolQualifierType_t type)
{
  if (mQualifier == BIOLOGICAL_QUALIFIER)
  {
    mBiolQualifier = type;
    mModelQualifier = BQM_UNKNOWN;
    mHasBeenModified = true;
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    mBiolQualifier = BQB_UNKNOWN;
    return LIBSBML_INVALID_ATTRIBUTE_VALUE;
  }
}

/*
 * set the model qualifier type
 * this should be consistent with the mQualifier == MODEL_QUALIFIER
 */
int 
CVTerm::setModelQualifierType(const std::string& qualifier)
{
  ModelQualifierType_t type = ModelQualifierType_fromString(qualifier.c_str());
  
  return setModelQualifierType(type);  
}


/*
 * set the biological qualifier type
 * this should be consistent with the mQualifier == BIOLOGICAL_QUALIFIER
 */
int 
CVTerm::setBiologicalQualifierType(const std::string& qualifier)
{
  BiolQualifierType_t type = BiolQualifierType_fromString(qualifier.c_str());
  return setBiologicalQualifierType(type);
}

/*
 * gets the Qualifier type
 */
QualifierType_t 
CVTerm::getQualifierType()
{
  return mQualifier;
}


QualifierType_t 
CVTerm::getQualifierType() const
{
  return mQualifier;
}


/*
 * gets the Model Qualifier type
 */
ModelQualifierType_t 
CVTerm::getModelQualifierType()
{
  return mModelQualifier;
}


ModelQualifierType_t 
CVTerm::getModelQualifierType() const
{
  return mModelQualifier;
}


/*
 * gets the biological Qualifier type
 */
BiolQualifierType_t 
CVTerm::getBiologicalQualifierType()
{
  return mBiolQualifier;
}


BiolQualifierType_t 
CVTerm::getBiologicalQualifierType() const
{
  return mBiolQualifier;
}


/*
 * gets the resources
 */
XMLAttributes * 
CVTerm::getResources()
{
  return mResources;
}

/*
 * gets the resources
 */
const XMLAttributes * 
CVTerm::getResources() const
{
  return mResources;
}


/*
 * Returns the number of resources for this %CVTerm.
 */
unsigned int 
CVTerm::getNumResources()
{
  return (unsigned int)mResources->getLength();
}

  
unsigned int 
CVTerm::getNumResources() const
{
  return (unsigned int)mResources->getLength();
}

  
/*
 * Returns the value of the nth resource for this %CVTerm.
 */
std::string
CVTerm::getResourceURI(unsigned int n)
{
  return mResources->getValue((int)n);
}

  
/*
 * Returns the value of the nth resource for this %CVTerm.
 */
std::string
CVTerm::getResourceURI(unsigned int n) const
{
  return mResources->getValue((int)n);
}

  
/*
 * adds a resource to the term
 */
int 
CVTerm::addResource(const std::string& resource)
{
  if (resource.empty())
  {
    return LIBSBML_OPERATION_FAILED;
  }
  else
  {
    mHasBeenModified = true;
    return mResources->addResource("rdf:resource", resource);
  }
}


/*
 * removes a resource to the term
 */
int 
CVTerm::removeResource(std::string resource)
{
  int result = LIBSBML_INVALID_ATTRIBUTE_VALUE;
  for (int n = 0; n < mResources->getLength(); n++)
  {
    if (resource == mResources->getValue(n))
    {
      mHasBeenModified = true;
      result = mResources->removeResource(n);
    }
  }

  if (mResources->getLength() == 0)
  {
    if (getQualifierType() == MODEL_QUALIFIER)
    {
      setModelQualifierType(BQM_UNKNOWN);
      setQualifierType(UNKNOWN_QUALIFIER);
    }
    else
    {
      setBiologicalQualifierType(BQB_UNKNOWN);
      setQualifierType(UNKNOWN_QUALIFIER);
    }
  }

  return result;
}

/** @cond doxygenLibsbmlInternal */
bool 
CVTerm::hasRequiredAttributes()
{
  bool valid = true;

  if (getQualifierType() == UNKNOWN_QUALIFIER)
  {
    valid = false;
  }
  else if (getQualifierType() == MODEL_QUALIFIER)
  {
    if (getModelQualifierType() == BQM_UNKNOWN)
    {
      valid = false;
    }
  }
  else
  {
    if (getBiologicalQualifierType() == BQB_UNKNOWN)
    {
      valid = false;
    }
  }

  if (valid)
  {
    if (getResources()->isEmpty())
    {
      valid = false;
    }
  }

  return valid;
}

bool 
CVTerm::hasRequiredAttributes() const
{
  bool valid = true;

  if (getQualifierType() == UNKNOWN_QUALIFIER)
  {
    valid = false;
  }
  else if (getQualifierType() == MODEL_QUALIFIER)
  {
    if (getModelQualifierType() == BQM_UNKNOWN)
    {
      valid = false;
    }
  }
  else
  {
    if (getBiologicalQualifierType() == BQB_UNKNOWN)
    {
      valid = false;
    }
  }

  if (valid)
  {
    if (getResources()->isEmpty())
    {
      valid = false;
    }
  }

  return valid;
}

bool
CVTerm::hasBeenModified()
{
  return mHasBeenModified;
}

void
CVTerm::resetModifiedFlags()
{
  mHasBeenModified = false;
}

void
CVTerm::setHasBeenModifiedFlag()
{
  mHasBeenModified = true;
}

void
CVTerm::setCapturedInStoredAnnotation(bool captured)
{
  mCapturedInStoredAnnotation = captured;
}

bool
CVTerm::getCapturedInStoredAnnotation()
{
   return mCapturedInStoredAnnotation;
}

unsigned int
CVTerm::getNumNestedCVTerms() const
{
  if (mNestedCVTerms != NULL)
  {
    return mNestedCVTerms->getSize();
  }
  else
  {
    return 0;
  }
}

const CVTerm *
CVTerm::getNestedCVTerm(unsigned int n) const
{
  return (mNestedCVTerms) ? static_cast <const CVTerm*> (mNestedCVTerms->get(n)) : NULL;
}


CVTerm *
CVTerm::getNestedCVTerm(unsigned int n)
{
  return (mNestedCVTerms) ? static_cast <CVTerm*> (mNestedCVTerms->get(n)) : NULL;
}


List*
CVTerm::getListNestedCVTerms()
{
  return mNestedCVTerms;
}

const List*
CVTerm::getListNestedCVTerms() const
{
  return const_cast<const List*>(mNestedCVTerms);
}


int
CVTerm::addNestedCVTerm(const CVTerm* term)
{
  if (term == NULL)
  {
    return LIBSBML_OPERATION_FAILED;
  }
  else if (!term->hasRequiredAttributes())
  {
    return LIBSBML_INVALID_OBJECT;
  }

  if (mNestedCVTerms == NULL)
  {
    mNestedCVTerms = new List();
  }

  unsigned int num = mNestedCVTerms->getSize();
  
  mNestedCVTerms->add((void *) term->clone());
  
  if (mNestedCVTerms->getSize() == num+1)
  {
    mHasBeenModified = true;
    return LIBSBML_OPERATION_SUCCESS;
  }
  else
  {
    return LIBSBML_OPERATION_FAILED;
  }
}


CVTerm*
CVTerm::removeNestedCVTerm(unsigned int n)
{
  CVTerm* removed = NULL;

  unsigned int num = mNestedCVTerms->getSize();

  if (n >= num)
  {
    return NULL;
  }

  removed = static_cast<CVTerm*>(mNestedCVTerms->remove(n));
  mHasBeenModified = true;

  return removed;
}
/** @endcond */
#endif /* __cplusplus */


/** @cond doxygenIgnored */
LIBSBML_EXTERN
CVTerm_t*
CVTerm_createWithQualifierType(QualifierType_t type)
{
  return new(nothrow) CVTerm(type);
}

LIBSBML_EXTERN
CVTerm_t*
CVTerm_createFromNode(const XMLNode_t *node)
{
  if (node == NULL) return NULL;
  return new(nothrow) CVTerm(*node);
}


LIBSBML_EXTERN
void
CVTerm_free(CVTerm_t * term)
{
  if (term == NULL ) return;
  delete static_cast<CVTerm*>(term);
}

LIBSBML_EXTERN
CVTerm_t *
CVTerm_clone (const CVTerm_t* c)
{
  if (c == NULL) return NULL;
  return static_cast<CVTerm*>( c->clone() );
}


LIBSBML_EXTERN
QualifierType_t 
CVTerm_getQualifierType(CVTerm_t * term)
{
  if (term == NULL) return UNKNOWN_QUALIFIER;
  return term->getQualifierType();
}

LIBSBML_EXTERN
ModelQualifierType_t 
CVTerm_getModelQualifierType(CVTerm_t * term)
{
  if (term == NULL) return BQM_UNKNOWN;
  return term->getModelQualifierType();
}

LIBSBML_EXTERN
BiolQualifierType_t 
CVTerm_getBiologicalQualifierType(CVTerm_t * term)
{
  if (term == NULL) return BQB_UNKNOWN;
  return term->getBiologicalQualifierType();
}

LIBSBML_EXTERN
XMLAttributes_t * 
CVTerm_getResources(CVTerm_t * term)
{
  if (term == NULL) return NULL;
  return term->getResources();
}

LIBSBML_EXTERN
unsigned int
CVTerm_getNumResources(CVTerm_t* term)
{
  if (term == NULL) return SBML_INT_MAX;
  return term->getNumResources();
}


LIBSBML_EXTERN
char *
CVTerm_getResourceURI(CVTerm_t * cv, unsigned int n)
{
  if (cv == NULL) return NULL;
  return cv->getResourceURI(n).empty() ? NULL : safe_strdup(cv->getResourceURI(n).c_str());
}


LIBSBML_EXTERN
int 
CVTerm_setQualifierType(CVTerm_t * term, QualifierType_t type)
{
  if (term == NULL) return LIBSBML_INVALID_OBJECT;
  return term->setQualifierType(type);
}


LIBSBML_EXTERN
int 
CVTerm_setModelQualifierType(CVTerm_t * term, ModelQualifierType_t type)
{
  if (term == NULL) return LIBSBML_INVALID_OBJECT;
  return term->setModelQualifierType(type);
}


LIBSBML_EXTERN
int 
CVTerm_setBiologicalQualifierType(CVTerm_t * term, BiolQualifierType_t type)
{
  if (term == NULL) return LIBSBML_INVALID_OBJECT;
  return term->setBiologicalQualifierType(type);
}

LIBSBML_EXTERN
int 
CVTerm_setModelQualifierTypeByString(CVTerm_t * term, const char* qualifier)
{
  if (term == NULL) return LIBSBML_INVALID_OBJECT;
  if (qualifier == NULL)
    return term->setModelQualifierType(BQM_UNKNOWN);
  else 
    return term->setModelQualifierType(qualifier);
}


LIBSBML_EXTERN
int 
CVTerm_setBiologicalQualifierTypeByString(CVTerm_t * term, const char* qualifier)
{
  if (term == NULL) return LIBSBML_INVALID_OBJECT;
  if (qualifier == NULL)
    return term->setBiologicalQualifierType(BQB_UNKNOWN);
  else 
    return term->setBiologicalQualifierType(qualifier);
}

LIBSBML_EXTERN
int 
CVTerm_addResource(CVTerm_t * term, const char * resource)
{
  if (term == NULL) return LIBSBML_OPERATION_FAILED;
  return term->addResource(resource);
}

LIBSBML_EXTERN
int 
CVTerm_removeResource(CVTerm_t * term, const char * resource)
{
  if (term == NULL) return LIBSBML_INVALID_OBJECT;
  return term->removeResource(resource);
}


LIBSBML_EXTERN
int
CVTerm_hasRequiredAttributes(CVTerm_t *cvt)
{
  if (cvt == NULL) return (int)false;
  return static_cast<int> (cvt->hasRequiredAttributes());
}


LIBSBML_EXTERN
unsigned int
CVTerm_getNumNestedCVTerms(const CVTerm_t *cvt)
{
  if (cvt == NULL) return SBML_INT_MAX;
  return cvt->getNumNestedCVTerms();
}


LIBSBML_EXTERN
const CVTerm_t*
CVTerm_getNestedCVTerm(const CVTerm_t *cvt, unsigned int n)
{
  if (cvt == NULL) return NULL;
  return cvt->getNestedCVTerm(n);
}



LIBSBML_EXTERN
int
CVTerm_addNestedCVTerm(CVTerm_t *cvt, const CVTerm_t* term)
{
  if (cvt == NULL) return LIBSBML_INVALID_OBJECT;
  return cvt->addNestedCVTerm(term);
}



LIBSBML_EXTERN
CVTerm_t*
CVTerm_removeNestedCVTerm(CVTerm_t *cvt, unsigned int n)
{
  if (cvt == NULL) return NULL;
  return (CVTerm_t *)(cvt->removeNestedCVTerm(n));
}



LIBSBML_EXTERN
const List_t*
CVTerm_getListNestedCVTerms(const CVTerm_t *cvt)
{
  if (cvt == NULL) return NULL;
  return cvt->getListNestedCVTerms();
}



static
const char* MODEL_QUALIFIER_STRINGS[] =
{
    "is"
  , "isDescribedBy"
  , "isDerivedFrom"
  , "isInstanceOf"
  , "hasInstance"
};

static
const char* BIOL_QUALIFIER_STRINGS[] =
{
    "is"
  , "hasPart"
  , "isPartOf"
  , "isVersionOf"
  , "hasVersion"
  , "isHomologTo"
  , "isDescribedBy"
  , "isEncodedBy"
  , "encodes"
  , "occursIn"
  , "hasProperty"
  , "isPropertyOf"
  , "hasTaxon"
};


LIBSBML_EXTERN
const char* 
ModelQualifierType_toString(ModelQualifierType_t type)
{
  int max = BQM_UNKNOWN;

  if (type < BQM_IS || type >= max)
  {
      return NULL;
  }

  return MODEL_QUALIFIER_STRINGS[type];
}

LIBSBML_EXTERN
const char* 
BiolQualifierType_toString(BiolQualifierType_t type)
{
  int max = BQB_UNKNOWN;

  if (type < BQB_IS || type >= max)
  {
      return NULL;
  }

  return BIOL_QUALIFIER_STRINGS[type];
}

LIBSBML_EXTERN
ModelQualifierType_t 
ModelQualifierType_fromString(const char* s)
{
  if (s == NULL) return BQM_UNKNOWN;

  int max = BQM_UNKNOWN;
  for (int i = 0; i < max; i++)
  {
    if (strcmp(MODEL_QUALIFIER_STRINGS[i], s) == 0)
      return (ModelQualifierType_t)i;
  }
  return BQM_UNKNOWN;
}

LIBSBML_EXTERN
BiolQualifierType_t 
BiolQualifierType_fromString(const char* s)
{  
  if (s == NULL) return BQB_UNKNOWN;

  int max = BQB_UNKNOWN;
  for (int i = 0; i < max; i++)
  {
    if (strcmp(BIOL_QUALIFIER_STRINGS[i], s) == 0)
      return (BiolQualifierType_t)i;
  }
  return BQB_UNKNOWN;
}
/** @endcond */

LIBSBML_CPP_NAMESPACE_END

