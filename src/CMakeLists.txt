# 
# The src folder contains the source code of the project, developed as a library
# called EDDI, short for Electron Density DIstribution calculator.  
#

message(STATUS "Configuring source code")

file(GLOB_RECURSE SOURCES "*.c")

add_library(EDDI SHARED ${SOURCES})
target_include_directories(EDDI PUBLIC ${INCLUDE_DIRECTORIES})
target_link_libraries(EDDI m)
