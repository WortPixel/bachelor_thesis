cmake_minimum_required(VERSION 2.6)

project( ICEPLAY )

set(CMAKE_CXX_COMPILER "clang++")

# Find Boost
find_package(Boost COMPONENTS program_options REQUIRED)

# Link root
exec_program("root-config"
       ARGS "--incdir"
       OUTPUT_VARIABLE ROOT_INC)

exec_program("root-config"
       ARGS "--glibs"
       OUTPUT_VARIABLE ROOT_LIB)

# include files (.h files) directories
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories(SYSTEM ${ROOT_INC} ${Boost_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -Wextra -pedantic")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${ROOT_LIB}")

add_library(ICEPLAY_LIB SHARED
	TDetector.cpp
	TSensor.cpp
	TEvent.cpp
	)

# Choose executables
add_executable(main
        main.cpp)

# Link Libraries
target_link_libraries(main ICEPLAY_LIB ${Boost_LIBRARIES})
