add_library(dslash)

#
# An AM_CONDITIONAL if statement to decide which fermion type
#

if( QDP_ARCH_SCALAR ) 
  target_sources(dslash PRIVATE
  	shift_table_scalar.cc 
  	cpp_dslash_scalar_32bit.cc
	cpp_dslash_3d_scalar_32bit.cc
	cpp_dslash_scalar_64bit.cc
	cpp_dslash_3d_scalar_64bit.cc
	cpp_clover_scalar_32bit.cc
	cpp_clover_scalar_64bit.cc
  ) 

  if(CPP_DSLASH_USE_OPENMP) 
    target_sources(dslash PRIVATE dispatch_scalar_openmp.cc)
  else()
	target_sources(dslash PRIVATE dispatch_scalar.cc)
  endif()  
endif()

if( QDP_ARCH_PARSCALAR ) 
  target_sources(dslash PRIVATE
    cpp_dslash_parscalar_32bit.cc  
    cpp_dslash_parscalar_64bit.cc
    cpp_dslash_parscalar_3d_32bit.cc
    cpp_dslash_parscalar_3d_64bit.cc
    cpp_dslash_parscalar_utils_32bit.cc
    cpp_dslash_parscalar_utils_64bit.cc
  )

 if(CPP_DSLASH_USE_OPENMP) 
    target_sources(dslash PRIVATE dispatch_parscalar_openmp.cc)
  else()
	target_sources(dslash PRIVATE dispatch_parscalar.cc)
  endif()

  target_link_libraries(dslash PUBLIC QMP::qmp)
endif()

if(CPP_DSLASH_USE_OPENMP)
	target_link_libraries(dslash PUBLIC OpenMP::OpenMP_CXX)
endif()

if( CPP_DSLASH_USE_QDPXX )
	target_sources(dslash PRIVATE qdp_packer_nopad.cc qdp_packer_3d.cc)
	target_link_libraries(dslash PUBLIC QDPXX::qdp)
endif()

target_include_directories(dslash PUBLIC 
			$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
			$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../include>
			$<INSTALL_INTERFACE:include>)
			



# Install the library and the targets
# into the export for the targets
install(TARGETS dslash EXPORT CPPWilsonDslashTargets
  ARCHIVE DESTINATION lib
  LIBRARY  DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

#install the export description of the targets
install(EXPORT CPPWilsonDslashTargets
                FILE CPPWilsonDslashTargets.cmake
                NAMESPACE CPPWilsonDslash::
                DESTINATION lib/cmake/CPPWilsonDslash)

add_library(CPPWilsonDslash::dslash ALIAS dslash)