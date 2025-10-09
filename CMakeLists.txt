CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

project(bib-studies)

# TESTING
include(CTest)
add_test(NAME ALLEGRO_o1_v03 COMMAND bash ../tests/full_chain_example_for_tests.sh)
