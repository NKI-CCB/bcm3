###############################################################################
# Description       : CMake build script for libSBML
# Original author(s): Frank Bergmann <fbergman@caltech.edu>
# Organization      : California Institute of Technology
#
# This file is part of libSBML.  Please visit http://sbml.org for more
# information about SBML, and the latest version of libSBML.
#
# Copyright (C) 2013-2018 jointly by the following organizations:
#     1. California Institute of Technology, Pasadena, CA, USA
#     2. EMBL European Bioinformatics Institute (EMBL-EBI), Hinxton, UK
#     3. University of Heidelberg, Heidelberg, Germany
#
# Copyright (C) 2009-2013 jointly by the following organizations:
#     1. California Institute of Technology, Pasadena, CA, USA
#     2. EMBL European Bioinformatics Institute (EMBL-EBI), Hinxton, UK
#
# Copyright (C) 2006-2008 by the California Institute of Technology,
#     Pasadena, CA, USA
#
# Copyright (C) 2002-2005 jointly by the following organizations:
#     1. California Institute of Technology, Pasadena, CA, USA
#     2. Japan Science and Technology Agency, Japan
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation.  A copy of the license agreement is provided
# in the file named "LICENSE.txt" included with this software distribution
# and also available online as http://sbml.org/software/libsbml/license.html
#
###############################################################################

set(LIBSBML_SKIP_SHARED_LIBRARY ON)
set(LIBSBML_SKIP_STATIC_LIBRARY ON)
set(WITH_CHECK OFF)
set(WITH_EXPAT ON)

###############################################################################
#
# Find all sources
#

macro(ADD_FUNCTION directory)

    set(prefix "sbml/")

    if (directory STREQUAL "sbml")
        set(prefix "")
    endif()


    file(GLOB temp
            ${prefix}${directory}/*.cpp
            ${prefix}${directory}/*.cxx
            ${prefix}${directory}/*.c
            ${prefix}${directory}/*.h)

    # Set the constraint files to be 'headers' as they don't have to be compiled
    # since they don't contain symbols but rather are included directly.  Could
    # remove them, but that keeps their content from being found in searches, etc.
    # in IDEs such as VisualStudio.
    foreach(tempFile ${temp})
        if (directory STREQUAL "validator/constraints")
            if (tempFile MATCHES ".*Constraints.cpp")
              set_source_files_properties(
                  ${tempFile}
                  PROPERTIES HEADER_FILE_ONLY true
                  )
#              list(REMOVE_ITEM temp "${tempFile}")
            endif()
        endif()
        if (tempFile MATCHES ".*.cxx$")
          set_source_files_properties(
              ${tempFile}
              PROPERTIES HEADER_FILE_ONLY true
              )
        endif()
    endforeach()

    # create source group for IDEs
    source_group(${directory} FILES ${temp})

    # add to libsbml sources
    set(LIBSBML_SOURCES ${LIBSBML_SOURCES} ${temp})

endmacro(ADD_FUNCTION)

foreach (directory annotation common conversion extension sbml math
                   units util validator validator/constraints)

	ADD_FUNCTION(${directory})
				       
endforeach()

set(XML_SOURCES)
set(XML_SOURCES ${XML_SOURCES}

  sbml/xml/XMLAttributes.cpp
  sbml/xml/XMLBuffer.cpp
  sbml/xml/XMLConstructorException.cpp
  sbml/xml/XMLError.cpp
  sbml/xml/XMLErrorLog.cpp
  sbml/xml/XMLLogOverride.cpp
  sbml/xml/XMLFileBuffer.cpp
  sbml/xml/XMLHandler.cpp
  sbml/xml/XMLInputStream.cpp
  sbml/xml/XMLMemoryBuffer.cpp
  sbml/xml/XMLNamespaces.cpp
  sbml/xml/XMLNode.cpp
  sbml/xml/XMLOutputStream.cpp
  sbml/xml/XMLParser.cpp
  sbml/xml/XMLToken.cpp
  sbml/xml/XMLTokenizer.cpp
  sbml/xml/XMLTriple.cpp
  sbml/xml/XMLAttributes.h
  sbml/xml/XMLBuffer.h
  sbml/xml/XMLConstructorException.h
  sbml/xml/XMLError.h
  sbml/xml/XMLErrorLog.h
  sbml/xml/XMLLogOverride.h
  sbml/xml/XMLFileBuffer.h
  sbml/xml/XMLHandler.h
  sbml/xml/XMLInputStream.h
  sbml/xml/XMLMemoryBuffer.h
  sbml/xml/XMLNamespaces.h
  sbml/xml/XMLNode.h
  sbml/xml/XMLOutputStream.h
  sbml/xml/XMLParser.h
  sbml/xml/XMLToken.h
  sbml/xml/XMLTokenizer.h
  sbml/xml/XMLTriple.h
)

if(WITH_EXPAT)

    set(XML_SOURCES ${XML_SOURCES}
        sbml/xml/ExpatAttributes.cpp
        sbml/xml/ExpatHandler.cpp
        sbml/xml/ExpatParser.cpp
        sbml/xml/ExpatAttributes.h
        sbml/xml/ExpatHandler.h
        sbml/xml/ExpatParser.h
    )
    set(LIBSBML_LIBS ${LIBSBML_LIBS} expat)

endif(WITH_EXPAT)

set (COMPRESS_SOURCES)
set(COMPRESS_SOURCES ${COMPRESS_SOURCES}
    sbml/compress/CompressCommon.h
    sbml/compress/CompressCommon.cpp
    sbml/compress/InputDecompressor.cpp
    sbml/compress/InputDecompressor.h
    sbml/compress/OutputCompressor.cpp
    sbml/compress/OutputCompressor.h
    )

source_group(xml FILES ${XML_SOURCES})
set(LIBSBML_SOURCES ${LIBSBML_SOURCES} ${XML_SOURCES} ${COMPRESS_SOURCES})

add_library (libsbml STATIC ${LIBSBML_SOURCES})
target_include_directories(libsbml PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_SOURCE_DIR}/dependencies/expat-2.6.0/lib)
target_link_libraries(libsbml ${LIBSBML_LIBS})
set_target_properties(libsbml PROPERTIES DEBUG_POSTFIX d)
set_target_properties(libsbml PROPERTIES FOLDER Dependencies)
