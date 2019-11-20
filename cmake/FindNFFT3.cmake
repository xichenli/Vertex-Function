# CMakeFind module to find (typically, user-installed) NFFT.

# Input variables:
#   NFFT3_ROOT pointing to directory containing NFFT installation
#   NFFT3_INCLUDES : (optional) guess to look for includes
#   NFFT3_LIBRARIES : (optional) location of NFFT libraries, if already known
# Output variables:
#   NFFT3_FOUND : true if NFFT is found
#   NFFT3_LIBRARIES : NFFT libraries
#   NFFT3_INCLUDE_DIRS : Directory with NFFT headers
#
  
include(FindPackageHandleStandardArgs)
include(SelectLibraryConfigurations)

# FIXME: try to use pkg-config hints.
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
find_path(NFFT3_INCLUDE_DIR nfft3.h nfft3_util.h PATHS ${NFFT3_ROOT} ${NFFT3_INCLUDES} PATH_SUFFIXES include)
set(NFFT3_INCLUDE_DIRS ${NFFT3_INCLUDE_DIR})
if (NOT NFFT3_LIBRARIES) 
  find_library(NFFT3_LIBRARY_RELEASE nfft3 PATHS ${NFFT3_ROOT} PATH_SUFFIXES lib)
  select_library_configurations(NFFT3)
endif()

find_package(FFTW3) # FIXME: make QUIET if requested
if (FFTW3_FOUND AND NFFT3_INCLUDE_DIRS AND NFFT3_LIBRARIES)
  list(APPEND NFFT3_LIBRARIES ${FFTW3_LIBRARIES})
  list(APPEND NFFT3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIRS})
endif()

find_package_handle_standard_args(NFFT3 DEFAULT_MSG 
				  NFFT3_LIBRARIES
				  NFFT3_INCLUDE_DIRS
				  FFTW3_FOUND)

