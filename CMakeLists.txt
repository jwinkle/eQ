

# Require CMake 3.5
cmake_minimum_required(VERSION 3.5)


set(PROJECT_NAME eQ)
project(${PROJECT_NAME})

# Set CMake behavior
cmake_policy(SET CMP0004 NEW)

# Get DOLFIN configuration data (DOLFINConfig.cmake must be in
# DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN REQUIRED)
find_package(Boost REQUIRED)

include(${DOLFIN_USE_FILE})

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# Do not throw error for 'multi-line comments' (these are typical in
# rst which includes LaTeX)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-Wno-comment" HAVE_NO_MULTLINE)
if (HAVE_NO_MULTLINE)
  set(CMAKE_CXX_FLAGS "-Wno-comment ${CMAKE_CXX_FLAGS}")
endif()


#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=1")

set(SRC_DIR ./src)
set(ABM_DIR ${SRC_DIR}/abm)

link_directories(./lib)

set(CPM_DIR ./Chipmunk-7.0.1/include)
set(HEADERS
#	fenics/hsl.h fenics/boundary.h
	${CPM_DIR}/chipmunk/chipmunk.h
	${CPM_DIR}/chipmunk/chipmunk_private.h)

# message( ${CPM_DIR} )
# message( ${HEADERS} )

# Executable
add_executable(${PROJECT_NAME}

	${HEADERS}

	${SRC_DIR}/main.cpp
	${SRC_DIR}/fHSL.cpp
	${SRC_DIR}/inputOutput.cpp
	${SRC_DIR}/simulation.cpp
	${SRC_DIR}/Strain.cpp

	${ABM_DIR}/eQabm.cpp
	${ABM_DIR}/eColi.cpp
	${ABM_DIR}/cpmHabitat.cpp
	${ABM_DIR}/cpmTrap.cpp
	${ABM_DIR}/cpmEColi.cpp

	diffuclass.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC ${CPM_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${SRC_DIR} fenics)

target_link_libraries(${PROJECT_NAME} libcp.a)

# Target libraries
target_link_libraries(
#	${PROJECT_NAME} dolfin mshr petsc
	${PROJECT_NAME} dolfin petsc
	${Boost_LIBRARIES} boost_system boost_filesystem
#	/usr/lib/x86_64-linux-gnu/libboost_system.a
#	/usr/lib/x86_64-linux-gnu/libboost_filesystem.a
	)
