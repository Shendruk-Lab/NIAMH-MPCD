# CMakeLists.txt for MPCD Simulation program
# Port of the original Makefile done by Timofey Kozhukhov

cmake_minimum_required(VERSION 3.10)
project(MPCD C)

# Set the compiler and compiler flags
set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_FLAGS "-Wall -O3")

# Additional definitions and conditional flags
add_definitions(-DPHASE_COUNT=2 -DTYPE_COUNT=3 -DPOLY_SETS=1 -DQ_SETS=4)
add_definitions(-DTHERMOSTAT_DPD)

# Gather source files and headers
file(GLOB_RECURSE SOURCES
        "mpcd/subroutines/*.c"
        "mpcd/*.c"
        "md/*.c"
        "dependencies/cJson/*.c"
        )
file(GLOB_RECURSE HEADERS
        "mpcd/headers/*.h"
        "md/*.h"
        "dependencies/cJson/*.h"
        )

# Create the executable target and link libraries
add_executable(mpcd ${SOURCES})
target_link_libraries(mpcd m)

### Custom targets (akin to Makefile phonies)
# note no clean command! just delete build directory
# debug
add_custom_target(debug
        COMMAND ${CMAKE_COMMAND} -E env opt=-g make -C ${CMAKE_BINARY_DIR}
        )

# debug+ (valgrind friendly)
add_custom_target(debug+
        COMMAND ${CMAKE_COMMAND} -E env opt=-ggdb3 make -C ${CMAKE_BINARY_DIR}
        )

# profiling
add_custom_target(prof
        COMMAND ${CMAKE_COMMAND} -E env opt=${CMAKE_C_FLAGS} -pg make -C ${CMAKE_BINARY_DIR}
        )

# doxygen
add_custom_target(docs
        COMMAND doxygen doxyfile
        )

# Custom target for building the project
add_custom_target(build_mpcd
        COMMAND make -C ${CMAKE_BINARY_DIR}
        DEPENDS mpcd
        )
