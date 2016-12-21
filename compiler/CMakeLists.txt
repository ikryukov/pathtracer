cmake_minimum_required(VERSION 3.1)
set (PROJECT PathTracer)
 
project (${PROJECT})
 
include_directories (
		../src/
		../3rdparty/glm/include/
		)

FILE(GLOB MyHeaders ../src/*.h)
FILE(GLOB MyCSources ../src/*.cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set (HEADERS
		${MyHeaders}
		)

set (SOURCES 
		${MyCSources})

link_directories (
		../3rdparty/GLEW/lib/
		)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_executable (${PROJECT} ${HEADERS} ${SOURCES})
target_link_libraries (${PROJECT})
