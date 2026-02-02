set(CONFIGFILES
    GEvolutionaryAlgorithm.json
    GFMinIndividual.json
    Go2.json
    GGradientDescent.json
    GParameterScan.json
    GSimulatedAnnealing.json
    GSwarmAlgorithm.json
)

# Copy the config files to the build tree for easier testing
file(COPY ${CONFIGFILES} DESTINATION ${CMAKE_BINARY_DIR}/config)

# Install config files
install(FILES ${CONFIGFILES} DESTINATION share/RhoNNO/config)
