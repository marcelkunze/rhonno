cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# When built via top-level rhonno/CMakeLists.txt this is a subproject;
# when configured alone (cmake -S RhoNNO) it is the main project.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(RhoNNO LANGUAGES CXX)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  endif()
endif()

# ROOT needs CMAKE_INSTALL_LIBDIR when CMAKE_LIBRARY_OUTPUT_DIRECTORY is set,
# otherwise ROOT_GENERATE_DICTIONARY() does install(FILES ...) with empty DESTINATION.
include(GNUInstallDirs)

# Place artifacts under bin/ and lib/ in the build tree (must be set before dict gen
# only after GNUInstallDirs so ROOT's install() rules get a real destination).
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ROOT (dict + runtime)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree TreePlayer Hist Gpad Graf MathCore)
include(${ROOT_USE_FILE})

# Compiler flags (valid GCC/Clang options only; -Wno-#warnings is not a real flag)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -pthread -Wunused-parameter")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow -Wno-unused-result -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -ffast-math -funroll-loops")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -march=native -ffast-math -funroll-loops -g")

set(LIB_SOURCES
  VNeuralNet.cxx
  VNeuralNetPlotter.cxx
  VSupervisedNet.cxx
  VUnsupervisedNet.cxx
  TFD.cxx
  TSGCS.cxx
  TSGNG.cxx
  TPerceptron.cxx
  TMLP.cxx
  TXMLP.cxx
  TXMLPInt8.cxx
  TXMLPInt16.cxx
  TNeuralNetCell.cxx
  TGCS.cxx
  TGNG.cxx
  TLVQ.cxx
  TDataServe.cxx
  TNNK.cxx
  TGNGTracker.cxx
  TRadon.cxx
  Graph.cxx
)

set(LIB_HEADERS
  VNeuralNet.h
  VNeuralNetPlotter.h
  VSupervisedNet.h
  VUnsupervisedNet.h
  TFD.h
  TSGCS.h
  TSGNG.h
  TPerceptron.h
  TMLP.h
  TXMLP.h
  TXMLPInt16.h
  TNeuralNetCell.h
  TGCS.h
  TGNG.h
  TLVQ.h
  TDataServe.h
  TNNK.h
  TGNGTracker.h
  TRadon.h
  Graph.h
)

# ROOT dictionary (produces RhoNNO_Cint.cxx + libRhoNNO_Cint_rdict.pcm + rootmap)
# Do not pass MODULE here: with some ROOT builds it injects install() without DESTINATION.
ROOT_GENERATE_DICTIONARY(RhoNNO_Cint ${LIB_HEADERS} LINKDEF RhoNNO_LinkDef.h)

add_library(RhoNNO SHARED ${LIB_SOURCES} RhoNNO_Cint.cxx)
target_include_directories(RhoNNO
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/RhoNNO>
)
target_link_libraries(RhoNNO PUBLIC ${ROOT_LIBRARIES})
set_target_properties(RhoNNO PROPERTIES
  OUTPUT_NAME RhoNNO
  VERSION 1.0.0
  SOVERSION 1
)

# Executables
foreach(_exe strain utrain NetworkTrainer NNOTracker RadonTracker test_xmlpint8 test_xmlpint16)
  add_executable(${_exe} ${_exe}.cxx)
  target_link_libraries(${_exe} PRIVATE RhoNNO ${ROOT_LIBRARIES})
endforeach()

# RadonOptimizer needs Geneva (optional, disabled by default)
option(RHONNO_BUILD_RADON_OPTIMIZER "Build RadonOptimizer (needs Geneva)" OFF)
if(RHONNO_BUILD_RADON_OPTIMIZER)
  add_executable(RadonOptimizer RadonOptimizer.cxx)
  target_link_libraries(RadonOptimizer PRIVATE RhoNNO ${ROOT_LIBRARIES})
endif()

# Optional: copy products into source-tree bin/ and lib/ after build
if(DEFINED RHONNO_SOURCE_BIN AND DEFINED RHONNO_SOURCE_LIB)
  add_custom_command(TARGET RhoNNO POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${RHONNO_SOURCE_LIB}
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:RhoNNO> ${RHONNO_SOURCE_LIB}/
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            $<TARGET_FILE_NAME:RhoNNO>
            ${RHONNO_SOURCE_LIB}/$<TARGET_SONAME_FILE_NAME:RhoNNO>
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            $<TARGET_SONAME_FILE_NAME:RhoNNO>
            ${RHONNO_SOURCE_LIB}/$<TARGET_LINKER_FILE_NAME:RhoNNO>
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libRhoNNO_Cint_rdict.pcm
            ${RHONNO_SOURCE_LIB}/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libRhoNNO_Cint.rootmap
            ${RHONNO_SOURCE_LIB}/
    COMMENT "Install libRhoNNO into ${RHONNO_SOURCE_LIB}"
  )
  foreach(_exe strain utrain NetworkTrainer NNOTracker RadonTracker test_xmlpint8 test_xmlpint16)
    add_custom_command(TARGET ${_exe} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory ${RHONNO_SOURCE_BIN}
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${_exe}> ${RHONNO_SOURCE_BIN}/
      COMMENT "Install ${_exe} into ${RHONNO_SOURCE_BIN}"
    )
  endforeach()
endif()

# Installation (pcm/rootmap already installed by ROOT_GENERATE_DICTIONARY)
install(TARGETS RhoNNO strain utrain NetworkTrainer NNOTracker RadonTracker test_xmlpint8 test_xmlpint16
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES ${LIB_HEADERS} TXMLPInt8.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/RhoNNO)

# Config JSON samples (if present)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists.txt)
  add_subdirectory(config)
endif()
