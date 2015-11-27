# CMakeFiles.txt for volregrid
#
# Andrew Janke - a.janke@gmail.com
#
# do this:
# export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/opt/minc-itk4/lib
#
# to ensure that CMake finds /opt/minc/lib/LIBMINCConfig.cmake

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)


PROJECT(mincfft)

# get current version from git tag
EXECUTE_PROCESS(COMMAND git describe
   OUTPUT_VARIABLE VERSION_STRING
   OUTPUT_STRIP_TRAILING_WHITESPACE)

ADD_DEFINITIONS(-DPACKAGE_NAME="${CMAKE_PROJECT_NAME}")
ADD_DEFINITIONS(-DPACKAGE_VERSION="${VERSION_STRING}")
ADD_DEFINITIONS(-DPACKAGE_BUGREPORT="a.janke@gmail.com")

# load any modules used in this project
SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")

# find the pre-requisites
FIND_PACKAGE(LIBMINC REQUIRED)
FIND_PACKAGE(FFTW REQUIRED)

# set compile options
INCLUDE( ${LIBMINC_USE_FILE} )
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
LINK_LIBRARIES(${LIBMINC_LIBRARIES} ${GSL_LIBRARIES})

ADD_EXECUTABLE(mincfft
   fft_support.h
   fft_support.c
   mincfft.c
   )

# what and where to install
INSTALL( TARGETS mincfft  DESTINATION bin)
