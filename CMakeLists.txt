# RhoNNO top-level (standalone ROOT build)
# FairROOT-style GENERATE_LIBRARY() is not used here.
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(rhonno LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Release if user did not pass -DCMAKE_BUILD_TYPE=...
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ROOT_GENERATE_DICTIONARY installs pcm/rootmap when CMAKE_LIBRARY_OUTPUT_DIRECTORY
# is set; GNUInstallDirs must be loaded first so DESTINATION is not empty.
include(GNUInstallDirs)

# Keep artifacts under the build tree in bin/ and lib/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Also mirror into the source tree bin/ and lib/ (legacy layout)
set(RHONNO_SOURCE_BIN ${CMAKE_SOURCE_DIR}/bin CACHE PATH "Legacy bin dir")
set(RHONNO_SOURCE_LIB ${CMAKE_SOURCE_DIR}/lib CACHE PATH "Legacy lib dir")

add_subdirectory(RhoNNO)

include(CTest)
if(BUILD_TESTING)
  # strain needs no args; utrain expects a data file
  add_test(NAME Teststrain COMMAND $<TARGET_FILE:strain>)
  add_test(
    NAME Testutrain
    COMMAND $<TARGET_FILE:utrain> ${CMAKE_SOURCE_DIR}/Data/ppe.dat
  )
  # NetworkTrainer steering file generated next to the binary tree
  set(_pid_nno ${CMAKE_BINARY_DIR}/pid.nno)
  file(WRITE ${_pid_nno}
"# Training of PIDSelectors with NNO
xmlp 7 15 10 1
transfer TR_FERMI
momentum 0.2
balance true
plots false
test 10000
start 1
stop 2
datapath ${CMAKE_SOURCE_DIR}/Data
networkpath ${CMAKE_SOURCE_DIR}/Networks
file PidTuple.root
tree PidTuple
input mom:acos(theta):svt:emc:drc:dch:ifr:ifrExp:ifrAdd
autoscale true
output abs(pid)==1
")
  add_test(NAME Testtrainer COMMAND $<TARGET_FILE:NetworkTrainer> ${_pid_nno})
  add_test(
    NAME Testtracker
    COMMAND $<TARGET_FILE:NNOTracker> ${CMAKE_SOURCE_DIR}/Data/event2
  )
  add_test(NAME Testradon COMMAND $<TARGET_FILE:RadonTracker>)
  add_test(NAME Testxmlpint8 COMMAND $<TARGET_FILE:test_xmlpint8>)
  add_test(NAME Testxmlpint16 COMMAND $<TARGET_FILE:test_xmlpint16>)
endif()
