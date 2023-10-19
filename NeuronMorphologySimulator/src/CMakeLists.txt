cmake_minimum_required(VERSION 3.24.2)
project(NeuronMorphologySimulator)
add_executable(NMorSim main.c)
install(TARGETS NMorSim DESTINATION bin)