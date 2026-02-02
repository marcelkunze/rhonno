cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(RhoNNO)

# Find ROOT
find_package(ROOT REQUIRED)

# Include ROOT headers
include(${ROOT_USE_FILE})

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler flags similar to the Makefile
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -pthread -m64 -Wunused-parameter -Wshadow")
# Suppress ROOT-related warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow -Wno-#warnings")

# Add optimization flags for better performance
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -ffast-math -funroll-loops")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -march=native -ffast-math -funroll-loops -g")

# Define source files for the library
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

# Generate ROOT dictionary
ROOT_GENERATE_DICTIONARY(RhoNNO_Cint ${LIB_HEADERS} LINKDEF RhoNNO_LinkDef.h)

# Create the library
add_library(RhoNNO SHARED ${LIB_SOURCES} RhoNNO_Cint.cxx)
target_link_libraries(RhoNNO ${ROOT_LIBRARIES})

# Set include directories
target_include_directories(RhoNNO PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Create executables
add_executable(strain strain.cxx)
target_link_libraries(strain RhoNNO ${ROOT_LIBRARIES})

add_executable(utrain utrain.cxx)
target_link_libraries(utrain RhoNNO ${ROOT_LIBRARIES})

add_executable(NetworkTrainer NetworkTrainer.cxx)
target_link_libraries(NetworkTrainer RhoNNO ${ROOT_LIBRARIES})

add_executable(NNOTracker NNOTracker.cxx)
target_link_libraries(NNOTracker RhoNNO ${ROOT_LIBRARIES})

add_executable(RadonTracker RadonTracker.cxx)
target_link_libraries(RadonTracker RhoNNO ${ROOT_LIBRARIES})

# RadonOptimizer requires Geneva library which may not be available
# add_executable(RadonOptimizer RadonOptimizer.cxx)
# target_link_libraries(RadonOptimizer RhoNNO ${ROOT_LIBRARIES})

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Installation
install(TARGETS RhoNNO strain utrain NetworkTrainer NNOTracker RadonTracker
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(FILES ${LIB_HEADERS} DESTINATION include/RhoNNO)

# Install ROOT pcm file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/RhoNNO_Cint_rdict.pcm DESTINATION lib)

# Add subdirectory for config
add_subdirectory(config)