# ---------------------------------------------------------------------------------
#     ooo      L   attice-based  |
#   o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
#  o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
# oo---0---oo  A   lgorithm and  |
#  o/./|\.\o   D   evelopment    | Copyright (C) 2013-2021 by
#   o/.|.\o    E   nvironment    | LeMonADE-Viewer Principal Developers (see AUTHORS)
#     ooo                        |
# ---------------------------------------------------------------------------------
#
# This file is part of LeMonADE.
#
# LeMonADE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LeMonADE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.
#
# --------------------------------------------------------------------------------
#
# Project Properties
#
CMAKE_MINIMUM_REQUIRED (VERSION 3.10.0)
PROJECT (LeMonADE-Viewer)
SET (APPLICATION_NAME "LeMonADE-Viewer")
SET (APPLICATION_CODENAME "${PROJECT_NAME}")
SET (APPLICATION_COPYRIGHT_YEARS "2015-2021")
SET (APPLICATION_VERSION_MAJOR 1)
SET (APPLICATION_VERSION_MINOR 3)
SET (APPLICATION_VERSION_PATCH 1)
SET (APPLICATION_VERSION_TYPE SNAPSHOT)
SET (APPLICATION_VERSION_STRING "${APPLICATION_VERSION_MAJOR}.${APPLICATION_VERSION_MINOR}.${APPLICATION_VERSION_PATCH}-${APPLICATION_VERSION_TYPE}")
SET (APPLICATION_ID "${APPLICATION_VENDOR_ID}.${PROJECT_NAME}")


SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fPIC -O3 -std=c++11 -fexpensive-optimizations -Wno-error=narrowing")

#
# Debugging Options
#
SET (CMAKE_VERBOSE_MAKEFILE 0) # Use 1 for debugging, 0 for release
SET (CMAKE_BUILD_TYPE "Release") #default build type is Release
#
# Check if LEMONADE_DIR is given
#
if( NOT LEMONADE_DIR )
set( LEMONADE_SUBMODULE_FOLDER "./extern/LeMonADE" )
option( PULL_LEMONADE "If enabled LeMonADE will be downloaded and compiled in the folder ${LEMONADE_SUBMODULE_FOLDER}" ON )
find_path( LEMONADE_INCLUDE_DIR
    NAMES
        LeMonADE/core/ConfigureSystem.h
        LeMonADE/core/Ingredients.h
        LeMonADE/utility/TaskManager.h
    PATHS
        ENV LEMONADE_ROOT
        ENV LEMONADE_INCLUDE_DIR
        ENV LEMONADE_DIR
        /usr/local/include
        /usr/include
        ${PROJECT_SOURCE_DIR}/extern
    PATH_SUFFIXES
        lemonade
        LeMonADE
        lemonade/include
        LeMonADE/include
        include
    DOC "LeMonADE include location"
)
set( LEMONADE_LIBRARY "LeMonADE" )
find_path( LEMONADE_LIBRARY_DIR
    NAMES
        "lib${LEMONADE_LIBRARY}.a"
    PATHS
        ENV LEMONADE_ROOT
        ENV LEMONADE_LIBRARY_DIR
        ENV LEMONADE_DIR
        /usr/local
        /usr
        ${PROJECT_SOURCE_DIR}/extern
    PATH_SUFFIXES
        LeMonADE/lib
        build/lib
        lib
    DOC "LeMonADE library location"
)
else()
    set( LEMONADE_INCLUDE_DIR "${LEMONADE_DIR}/include" )
    set( LEMONADE_LIBRARY_DIR "${LEMONADE_DIR}/lib" )
endif()
if( NOT LEMONADE_INCLUDE_DIR )
    message( "LEMONADE_INCLUDE_DIR is not provided. If build fails, use -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE/headers/ or install to default location" )
endif()
if( NOT LEMONADE_LIBRARY_DIR )
    message( "LEMONADE_LIBRARY_DIR is not provided. If build fails, use -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE/lib/ or install to default location" )
endif()
if( ( NOT LEMONADE_INCLUDE_DIR ) OR ( NOT LEMONADE_LIBRARY_DIR ) )
    if( PULL_LEMONADE )
        message( "Installing LeMonADE ..." )
        # every command runs in its own process, that's why this workaround
        # of calling sh with a script was used
        execute_process(
            # installing into the source is very unclean, but as 80% is
            # header only it will simply result in only the compiled
            # library being 'installed'
            COMMAND sh -c "\
            git submodule update --init -- ${LEMONADE_SUBMODULE_FOLDER}     && \
            cd ${LEMONADE_SUBMODULE_FOLDER}                                 && \
            mkdir -p build install                                          && \
            cd build                                                        && \
            cmake -DCMAKE_CXX_COMPILER='${CMAKE_CXX_COMPILER}'                 \
                  -DINSTALLDIR_LEMONADE=.. DCMAKE_INSTALL_PREFIX=.. ..      && \
            make install                                                       \
            "
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
        if ( NOT EXISTS "${PROJECT_SOURCE_DIR}/${LEMONADE_SUBMODULE_FOLDER}/lib/lib${LEMONADE_LIBRARY}.a" )
            message( FATAL_ERROR "Something went wrong when trying to install LeMonADE from source!" )
        endif()
        set( LEMONADE_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/${LEMONADE_SUBMODULE_FOLDER}/include" )
        set( LEMONADE_LIBRARY_DIR "${PROJECT_SOURCE_DIR}/${LEMONADE_SUBMODULE_FOLDER}/lib" )
    else()
        message( "You can also call cmake again with the option -DPULL_LEMONADE=ON in order to try to automatically download and compile LeMonADE into a default subfolder." )
    endif()
else()
    message( "Found LeMonADE at ${LEMONADE_INCLUDE_DIR} and ${LEMONADE_LIBRARY_DIR}" )
endif()

#
# Project Output Paths
#
#MESSAGE(${PROJECT_SOURCE_DIR})
#SET (MAINFOLDER ${PROJECT_SOURCE_DIR})
SET (LEMONADEVIEWER_DIR ${PROJECT_SOURCE_DIR})
#MESSAGE(${LEMONADE_DIR})
SET (EXECUTABLE_OUTPUT_PATH "${LEMONADEVIEWER_DIR}/build/")
#SET (LIBRARY_OUTPUT_PATH "${LEMONADE_DIR}/lib")

#
# Project Search Paths
#
LIST (APPEND CMAKE_PREFIX_PATH "${LEMONADEVIEWER_DIR}")
#LIST (APPEND CMAKE_PREFIX_PATH "${LEMONADEVIEWER_DIR}/tools")
#LIST (APPEND CMAKE_PREFIX_PATH "${LEMONADEVIEWER_DIR}/extern")
#SET (CMAKE_MODULE_PATH "${LEMONADE_DIR}/tools/share/cmake")
INCLUDE_DIRECTORIES("${LEMONADEVIEWER_DIR}/include")

INCLUDE_DIRECTORIES( "${LEMONADE_INCLUDE_DIR} ")
link_directories( "${LEMONADE_LIBRARY_DIR}" )


find_package(OpenGL REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS})

set(FLTK_SKIP_FLUID TRUE)
find_package(FLTK REQUIRED)
include_directories( ${FLTK_INCLUDE_DIRS})

#
# Add Build Targets
#
ADD_SUBDIRECTORY(src)



#
# Add Install Targets
#
IF (EXISTS "${LEMONADEVIEWER_DIR}/include/${PROJECT_NAME}" AND IS_DIRECTORY "${LEMONADEVIEWER_DIR}/include/${PROJECT_NAME}")
    INSTALL(DIRECTORY "${MAINFOLDER}/include/${PROJECT_NAME}" DESTINATION "include")
ENDIF(EXISTS "${LEMONADEVIEWER_DIR}/include/${PROJECT_NAME}" AND IS_DIRECTORY "${LEMONADEVIEWER_DIR}/include/${PROJECT_NAME}")


