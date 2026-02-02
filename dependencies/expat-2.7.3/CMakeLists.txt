add_library(expat)

target_compile_definitions(expat PRIVATE XML_GE=1)

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
