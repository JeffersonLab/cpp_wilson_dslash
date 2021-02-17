# Boilerplate 
cmake_minimum_required(VERSION 3.17)
project(cpp_wilson_dslash VERSION 1.0.0 LANGUAGES CXX C)

# --------------------------------------------------
# Cache file options
#     Things one can set in CCMake or CMakeGUI 
#     or on the command line
# --------------------------------------------------
set(CPP_DSLASH_PARALLEL_ARCH "scalar" CACHE STRING 
	"The parallel arch: values are scalar and parscalar. Default is scalar")

# This part restricts the selection in the GUIs to scalar 
# and parscalar. One can still enter any odd rubbish tho
set_property(CACHE CPP_DSLASH_PARALLEL_ARCH PROPERTY
			STRINGS scalar parscalar)
			
option(CPP_DSLASH_USE_OPENMP "Use OpenMP Threads" OFF) 
option(CPP_DSLASH_USE_SSE2 "USE SSE2" OFF )

option(CPP_DSLASH_NOCOMMS "Disable Communications" OFF)
option(CPP_DSLASH_NOCOMPUTE "Disable Computation"  OFF)
option(CPP_DSLASH_USE_QDPXX "Use QDP" OFF)

# Configuration based on the options: Parallel arch
# Convert user entry to lower so we only need to match
string(TOLOWER ${CPP_DSLASH_PARALLEL_ARCH} dslashParallelArch)
if( dslashParallelArch STREQUAL "scalar")
  message(STATUS "CPP_Dslash: using scalar arch")
  set(CPP_DSLASH_SCALAR 1)      # Macro in Config internal.h
elseif( dslashParallelArch STREQUAL "parscalar")
  message(STATUS "CPP_Dslash: using parscalar arch")
  set(CPP_DSLASH_PARSCALAR 1)   # Macro in  Config internal.h
endif()

# Whether to use SSE2
if( CPP_DSLASH_USE_SSE2 ) 
  message(STATUS "Setting Dslash SSE2")
  set(DSLASH_USE_SSE2 1)    # Convert to value used by header
endif()

# Whether to use OpenMP
if( CPP_DSLASH_USE_OPENMP ) 
  set(DSLASH_USE_OMP_THREADS 1)
  find_package(OpenMP REQUIRED)
  find_package(Threads REQUIRED)
endif()

# Whether to use QDPXX
if( CPP_DSLASH_USE_QDPXX ) 
  set(SSE_USE_QDPXX 1)
  find_package(QDPXX REQUIRED)
  if( NOT QDP_PARALLEL_ARCH STREQUAL dslashParallelArch)
    message(ERROR "QDPXX is used but has different parallel arch. QDPXX Parallel arch is ${QDP_PARALLEL_ARCH}.") 
  endif()
endif()

# Whether to enable nocomute and nocomms 
# options. Do these still do anything useful
if( CPP_DSLASH_NOCOMPUTE )
  message(STATUS "Disabling Compute")
  set(SSEDSLASH_4D_NOCOMPUTE 1)
endif()

if( CPP_DSLASH_NOCOMMS )
  message(STATUS "Disabling Comms")
  set(SSEDSLASH_4D_NOCOMMS 1)
endif()

# Write the config_internal.h file
configure_file(include/dslash_config_internal.h.cmake.in 
			   include/dslash_config_internal.h)

# Install it? (Or make it a public header later?
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/dslash_config_internal.h DESTINATION include)
  
add_subdirectory(lib)
add_subdirectory(tests)

#install the headers 
install(DIRECTORY include DESTINATION . 
	    FILES_MATCHING PATTERN "*.h")
	    
	   # Make the config file
include(CMakePackageConfigHelpers)
configure_package_config_file(Config.cmake.in CPPWilsonDslashConfig.cmake 
							  INSTALL_DESTINATION lib/cmake/CPPWilsonDslash)

# Make the version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    CPPWilsonDslashVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CPPWilsonDslashVersion.cmake 
			  ${CMAKE_CURRENT_BINARY_DIR}/CPPWilsonDslashConfig.cmake 
			  DESTINATION lib/cmake/CPPWilsonDslash)
	    
