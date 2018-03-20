# Create a library called "RhoNNO" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.
#add_subdirectory(DecayTreeFitter)
#add_subdirectory(RhoFitter)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${SYSTEM_INCLUDE_DIRECTORIES} 
  ${BASE_INCLUDE_DIRECTORIES} 
)

Set(INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}/rhonno
)

Include_Directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

set(SRCS 
RhoNNO/VNeuralNetPlotter.cxx
RhoNNO/VNeuralNet.cxx
RhoNNO/VSupervisedNet.cxx
RhoNNO/VUnsupervisedNet.cxx 
RhoNNO/TFD.cxx 
RhoNNO/TSGCS.cxx 
RhoNNO/TSGNG.cxx 
RhoNNO/TPerceptron.cxx 
RhoNNO/TMLP.cxx 
RhoNNO/TXMLP.cxx
RhoNNO/TNeuralNetCell.cxx 
RhoNNO/TGCS.cxx 
RhoNNO/TGNG.cxx 
RhoNNO/TGNGTracker.cxx
RhoNNO/TLVQ.cxx 
RhoNNO/TDataServe.cxx 
RhoNNO/TNNK.cxx
)

set(LINKDEF  RhoNNO/RhoNNO_LinkDef.h)
set(LIBRARY_NAME RhoNNO)

set(DEPENDENCIES Base )

GENERATE_LIBRARY()

# add test programs
add_executable (strain RhoNNO/strain.cxx)
target_link_libraries (strain  RhoNNO ${ROOT_LIBRARIES} TreePlayer)

add_executable (utrain RhoNNO/utrain.cxx)
target_link_libraries (utrain  RhoNNO ${ROOT_LIBRARIES} TreePlayer)

add_executable (NetworkTrainer RhoNNO/NetworkTrainer.cxx)
target_link_libraries (NetworkTrainer  RhoNNO ${ROOT_LIBRARIES} TreePlayer)

add_executable (NNOTracker RhoNNO/NNOTracker.cxx)
target_link_libraries (NNOTracker  RhoNNO ${ROOT_LIBRARIES} TreePlayer)

include(CTest)

# Run Tests
add_test (Teststrain ${CMAKE_BINARY_DIR}/bin/strain)
add_test (Testutrain ${CMAKE_BINARY_DIR}/bin/utrain ${CMAKE_SOURCE_DIR}/rhonno/Data/ppe.dat)
add_test (Testtrainer ${CMAKE_BINARY_DIR}/bin/NetworkTrainer ${CMAKE_BINARY_DIR}/pid.nno)
add_test (Testtracker ${CMAKE_BINARY_DIR}/bin/NNOTracker ${CMAKE_SOURCE_DIR}/rhonno/Data/event2)
file(WRITE ${CMAKE_BINARY_DIR}/pid.nno "
# Training of PIDSelectors with NNO

#define the network topology
xmlp 7 15 10 1
transfer TR_FERMI
momentum 0.2
balance true
plots true
test 10000
start 1
stop 10

#define the data source
datapath ${CMAKE_SOURCE_DIR}/rhonno/Data
networkpath ${CMAKE_SOURCE_DIR}/rhonno/Networks
file PidTuple.root

#set up the input layer (use branch names)
tree PidTuple
input mom:acos(theta):svt:emc:drc:dch:ifr:ifrExp:ifrAdd
autoscale true

#set up the output layer (use branch names)
#Particles pid = {electron=1,muon,pion,kaon,proton}
output abs(pid)==1

#end of file
")
