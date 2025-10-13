CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

project(bib-studies)

include(GNUInstallDirs)

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/tests DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME})

# TESTING
include(CTest)
add_test(NAME ALLEGRO_o1_v03 COMMAND bash ../tests/full_chain_example_for_tests.sh)
