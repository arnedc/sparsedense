project(sparsedense)

cmake_minimum_required(VERSION 2.6) 

enable_language(Fortran)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/src)


aux_source_directory(src helpfunctions)
add_executable(sparsedense ${helpfunctions} main.cpp)


  set_target_properties(sparsedense PROPERTIES
    LINK_FLAGS "-L/opt/intel/mkl/lib/intel64 -L/opt/intel/lib/intel64 -L/home/arnedc/Libraries/")


find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

target_link_libraries(sparsedense
                      #${LAPACK_LIBRARIES}
                      #${BLAS_LIBRARIES}
                      -lpardiso500-GNU461-X86-64
                      -lgfortran
                      ${MPI_LIBRARIES}
		      -lmkl_scalapack_lp64   # MKL libraries might be replaced with others
		      -lmkl_blas95_lp64
		      -lmkl_lapack95_lp64
		      -lmkl_intel_lp64
		      -lmkl_sequential
		      -lmkl_core
                      -lmkl_blacs_intelmpi_lp64   
		      -fopenmp
		      -lm -lpthread)
