# ---------------------------------------------------------------------------------
#     ooo      L   attice-based  |
#   o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
#  o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
# oo---0---oo  A   lgorithm and  |
#  o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
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
CMAKE_MINIMUM_REQUIRED (VERSION 2.6.2)
PROJECT (LeMonADE-Viewer)
SET (APPLICATION_NAME "LeMonADE-Viewer")
SET (APPLICATION_CODENAME "${PROJECT_NAME}")
SET (APPLICATION_COPYRIGHT_YEARS "2015")
SET (APPLICATION_VERSION_MAJOR 1)
SET (APPLICATION_VERSION_MINOR 2)
SET (APPLICATION_VERSION_PATCH 1)
SET (APPLICATION_VERSION_TYPE SNAPSHOT)
SET (APPLICATION_VERSION_STRING "${APPLICATION_VERSION_MAJOR}.${APPLICATION_VERSION_MINOR}.${APPLICATION_VERSION_PATCH}-${APPLICATION_VERSION_TYPE}")
SET (APPLICATION_ID "${APPLICATION_VENDOR_ID}.${PROJECT_NAME}")


#
# Debugging Options
#
SET (CMAKE_VERBOSE_MAKEFILE 0) # Use 1 for debugging, 0 for release

#
# Check if LEMONADE_DIR is given
#
if (DEFINED LEMONADE_DIR)
	message("LEMONADE_DIR set to " ${LEMONADE_DIR})
else()
    message( WARNING "LEMONADE_DIR no set! If the LeMonADE library was installed to a custom directory, use -DLEMONADE_DIR=/path/to/LeMonADE/ to set the correct search path")
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

INCLUDE_DIRECTORIES("${LEMONADE_DIR}/include")
link_directories ("${LEMONADE_DIR}/lib/")


find_package(OpenGL REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS})

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


