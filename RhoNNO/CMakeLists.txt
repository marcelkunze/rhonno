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
  TXMLPNPU.cxx
  TXMLPGPU.cxx
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

# Optional ONNX Runtime / AMD Ryzen NPU support for TXMLPNPU.
# Point RHONNO_ORT_ROOT at the Ryzen AI C++ install (contains onnxruntime/{include,lib}).
option(RHONNO_WITH_NPU "Build TXMLPNPU against ONNX Runtime + VitisAI EP" ON)
set(RHONNO_ORT_ROOT "$ENV{HOME}/ryzen-ai-cpp/onnxruntime" CACHE PATH
    "Root of the ONNX Runtime SDK shipped with Ryzen AI (has include/ and lib/)")
set(RHONNO_VAIP_CONFIG "$ENV{HOME}/ryzen-ai-cpp/voe-4.0-linux_x86_64/vaip_config.json" CACHE FILEPATH
    "Path to vaip_config.json for the VitisAI EP")
if(NOT RHONNO_WITH_NPU OR NOT EXISTS "${RHONNO_ORT_ROOT}/include/onnxruntime_cxx_api.h")
  # Without ORT the source cannot compile. Drop it before add_library happens.
  list(REMOVE_ITEM LIB_SOURCES TXMLPNPU.cxx)
  message(WARNING "NPU support disabled — TXMLPNPU.cxx excluded from RhoNNO.")
endif()

# Optional HIP/hisBLAS GPU support for TXMLPGPU (AMD ROCm)
option(RHONNO_WITH_GPU "Build TXMLPGPU against HIP + hipBLAS (ROCm)" ON)
if(RHONNO_WITH_GPU)
  find_program(HIPCC hipcc)
  find_path(HIP_INCLUDE_DIR hip/hip_runtime.h HINTS /usr/include /opt/rocm/include)
  find_library(HIP_LIB amdhip64 HINTS /usr/lib/x86_64-linux-gnu /opt/rocm/lib)
  if(HIP_INCLUDE_DIR AND HIP_LIB)
    message(STATUS "GPU support: HIP=${HIP_LIB} (no hipBLAS — custom kernel)")
    set(RHONNO_HAS_GPU TRUE)
  else()
    message(WARNING "GPU support disabled — HIP not found.")
    set(RHONNO_HAS_GPU FALSE)
    list(REMOVE_ITEM LIB_SOURCES TXMLPGPU.cxx)
  endif()
else()
  set(RHONNO_HAS_GPU FALSE)
  list(REMOVE_ITEM LIB_SOURCES TXMLPGPU.cxx)
endif()

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

if(RHONNO_WITH_NPU AND EXISTS "${RHONNO_ORT_ROOT}/include/onnxruntime_cxx_api.h")
  target_include_directories(RhoNNO PRIVATE "${RHONNO_ORT_ROOT}/include")
  target_compile_definitions(RhoNNO PRIVATE RHONNO_WITH_NPU=1
    RHONNO_VAIP_CONFIG_PATH="${RHONNO_VAIP_CONFIG}")
  target_link_directories(RhoNNO PUBLIC "${RHONNO_ORT_ROOT}/lib")
  target_link_libraries(RhoNNO PUBLIC onnxruntime)
  message(STATUS "NPU support: using ONNX Runtime at ${RHONNO_ORT_ROOT}")
endif()

if(RHONNO_HAS_GPU)
  # TXMLPGPU.cxx uses __global__ kernels → needs HIP language, not plain CXX
  enable_language(HIP)
  # Radeon 890M (Strix Point / RDNA3.5) = gfx1150
  set(CMAKE_HIP_ARCHITECTURES "gfx1150" CACHE STRING "HIP target architecture" FORCE)
  set_source_files_properties(TXMLPGPU.cxx PROPERTIES LANGUAGE HIP)
  target_compile_definitions(RhoNNO PRIVATE RHONNO_WITH_GPU=1)
  target_include_directories(RhoNNO PRIVATE "${HIP_INCLUDE_DIR}")
  target_link_libraries(RhoNNO PUBLIC ${HIP_LIB})
endif()

set_target_properties(RhoNNO PROPERTIES
  OUTPUT_NAME RhoNNO
  VERSION 1.0.0
  SOVERSION 1
)

# Executables
foreach(_exe strain utrain NetworkTrainer NNOTracker RadonTracker test_xmlpint8 test_xmlpint16 test_xmlpnpu test_xmlpgpu)
  if(NOT RHONNO_WITH_NPU AND _exe STREQUAL "test_xmlpnpu")
    continue()
  endif()
  if(NOT RHONNO_HAS_GPU AND _exe STREQUAL "test_xmlpgpu")
    continue()
  endif()
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
