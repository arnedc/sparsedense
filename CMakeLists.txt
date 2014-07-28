project(sparsedense)

cmake_minimum_required(VERSION 2.6) 

enable_language(Fortran)

find_package(MPI REQUIRED)
#find_package(HDF5 REQUIRED)
include_directories(${MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/src)
#include(CMakeForceCompiler)
#  CMAKE_FORCE_C_COMPILER(mpicc ${CMAKE_C_COMPILER_ID})
#  CMAKE_FORCE_CXX_COMPILER(mpicxx ${CMAKE_CXX_COMPILER_ID})

aux_source_directory(src helpfunctions)
add_executable(sparsedense ${helpfunctions} main.cpp)

#if(MPI_COMPILE_FLAGS)
#  set_target_properties(pardiso_seq PROPERTIES
#    COMPILE_FLAGS "-DMKL_ILP64 -mpe=mpilog ${MPI_COMPILE_FLAGS}")
#endif()

#if(MPI_LINK_FLAGS)
  set_target_properties(sparsedense PROPERTIES
    LINK_FLAGS "-L/opt/intel/mkl/lib/intel64 -L/opt/intel/lib/intel64 -L/home/arnedc/Libraries/")
#endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

target_link_libraries(sparsedense
                      #${LAPACK_LIBRARIES}
                      #${BLAS_LIBRARIES}
                      -lpardiso500-GNU461-X86-64
                      -lgfortran
                      ${MPI_LIBRARIES}
		      #${HDF5_LIBRARIES}
		      -lmkl_scalapack_lp64   # dit is de MKL variant van ScaLapack
		      -lmkl_blas95_lp64
		      -lmkl_lapack95_lp64
		      -lmkl_intel_lp64
		      -lmkl_sequential
		      -lmkl_core
                      -lmkl_blacs_intelmpi_lp64     # dit is de MKL variant van BLACS
		      #-liomp5
		      -fopenmp
		      -lm -lpthread)



