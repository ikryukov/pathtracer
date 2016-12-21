cmake_minimum_required(VERSION 2.8)
set (PROJECT PathTracer)
 
project (${PROJECT})
 
include_directories (
		../src/
		../3rdparty/GLEW/include/
		../3rdparty/glm/include/
		)

if (APPLE)
	set (LIBRARIES
			glew
	)
endif (APPLE)

if (WIN32)
	set (LIBRARIES
		glew32
		glu32
		opengl32
	)
endif (WIN32)

if (APPLE)
    find_package (OpenGL REQUIRED)
    find_package (GLUT REQUIRED)
    include_directories(/usr/local/include)
    find_library(COCOA_LIBRARY Cocoa REQUIRED)
    find_library(IOKIT_LIBRARY IOKit REQUIRED)
    find_library(COREVID_LIBRARY CoreVideo REQUIRED)
endif (APPLE)

FILE(GLOB MyHeaders ../src/*.h)
FILE(GLOB MyCSources ../src/*.cpp)

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
target_link_libraries (${PROJECT} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${LIBRARIES} ${COCOA_LIBRARY} ${COREVID_LIBRARY} ${IOKIT_LIBRARY})
