macro(_expat_copy_bool_int source_ref dest_ref)
    if(${source_ref})
        set(${dest_ref} 1)
    else()
        set(${dest_ref} 0)
    endif()
endmacro()
_expat_copy_bool_int(EXPAT_ATTR_INFO        XML_ATTR_INFO)
_expat_copy_bool_int(EXPAT_DTD              XML_DTD)
_expat_copy_bool_int(EXPAT_GE               XML_GE)
_expat_copy_bool_int(EXPAT_LARGE_SIZE       XML_LARGE_SIZE)
_expat_copy_bool_int(EXPAT_MIN_SIZE         XML_MIN_SIZE)
_expat_copy_bool_int(EXPAT_NS               XML_NS)
set(XML_CONTEXT_BYTES 1024)
set(_SHARED STATIC)

configure_file(expat_config.h.cmake ${CMAKE_SOURCE_DIR}/dependencies/expat-2.6.0/lib/expat_config.h)

add_library(expat)

target_sources(expat
  PRIVATE
    lib/xmlparse.c
    lib/xmlrole.c
    lib/xmltok.c
)

target_include_directories(expat
  PRIVATE
    lib
)

set_target_properties(expat PROPERTIES DEBUG_POSTFIX d)
set_target_properties(expat PROPERTIES FOLDER Dependencies)
#if(MSVC)
#    add_definitions(-D_CRT_SECURE_NO_WARNINGS -wd4996)
#endif()
