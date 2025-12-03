# load the ff modules first before building
# module load cmake/3.27.9
# module load nvhpc/24.5
# module load hpcx-mpi/2.19
# module load hdf5/1.14.3--hpcx-mpi--2.19--nvhpc--24.5
# module load netcdf-fortran/4.6.1--hpcx-mpi--2.19--nvhpc--24.5
# module load binutils/2.42
# module load python/3.11
# module load parallel-netcdf/1.12.3--hpcx-mpi--2.19--nvhpc--24.5

cmake_minimum_required(VERSION 3.18)
project(AtmosModel LANGUAGES Fortran)

# Option to control OpenACC vs OpenMP build (like USE_OPENACC in the Makefile)
option(USE_OPENACC "Enable OpenACC (NVHPC) instead of OpenMP" OFF)

# ------------------------------------------------------------
# NetCDF (Fortran + C) flags via nf-config / nc-config
# ------------------------------------------------------------
# These commands mirror:
#   FCFLAGS_BASE = -I. `nf-config --fflags`
#   LDFLAGS      = -L`nc-config --libdir` `nf-config --flibs`
#
# Assumes nf-config and nc-config are in PATH on the cluster.
execute_process(
    COMMAND nf-config --fflags
    OUTPUT_VARIABLE NF_FFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND nf-config --flibs
    OUTPUT_VARIABLE NF_FLIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND nc-config --libdir
    OUTPUT_VARIABLE NC_LIBDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Convert the strings to CMake lists
separate_arguments(NF_FFLAGS)
separate_arguments(NF_FLIBS)

# ------------------------------------------------------------
# MPI
# ------------------------------------------------------------
# On your cluster, CMAKE_Fortran_COMPILER will usually be set to mpif90/mpifort
# by the environment (e.g. module system). We still use find_package for safety.
find_package(MPI REQUIRED Fortran)

# ------------------------------------------------------------
# Source files
# ------------------------------------------------------------
set(MODEL_SOURCES
    module_parameters.f90
    module_types.F90
    module_output.F90
    module_physics.f90
    parallel_timer.f90
    model.F90
)

add_executable(model ${MODEL_SOURCES})

# ------------------------------------------------------------
# Compile options (mirroring the Makefile)
# ------------------------------------------------------------
# Base NetCDF flags
target_compile_options(model PRIVATE ${NF_FFLAGS})

# Local include (equivalent to -I.)
target_include_directories(model PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Choose between OpenACC and OpenMP-style build
if(USE_OPENACC)
    message(STATUS "Building with OpenACC (NVHPC)")
    # Equivalent of:
    #   FCFLAGS = -fast -acc -gpu=managed -Minfo=accel -Mpreprocess $(FCFLAGS_BASE) -DUSE_OPENACC
    target_compile_definitions(model PRIVATE USE_OPENACC)
    target_compile_options(model PRIVATE
        -fast
        -acc
        -gpu=managed
	-Minfo=accel
        -Mpreprocess
    )
else()
    message(STATUS "Building with OpenMP-style flags")
    # Equivalent of:
    #   FCFLAGS = -Ofast $(FCFLAGS_BASE) -fPIC -fopenmp -DUSE_OPENMP -cpp
    target_compile_definitions(model PRIVATE USE_OPENMP)
    target_compile_options(model PRIVATE
        -Ofast
        -fPIC
        -fopenmp
        -cpp
    )

    # Optional: also use CMake's OpenMP discovery if available
    find_package(OpenMP)
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(model PRIVATE OpenMP::OpenMP_Fortran)
    endif()
endif()

# ------------------------------------------------------------
# Link options / libraries
# ------------------------------------------------------------
# NetCDF Fortran & C libs from nf-config / nc-config
# NF_FLIBS typically contains things like: -L<dir> -lnetcdff -lnetcdf ...
link_directories(${NC_LIBDIR})
target_link_libraries(model PRIVATE ${NF_FLIBS} MPI::MPI_Fortran)

# ------------------------------------------------------------
# Test setup (replicates `make test`)
# ------------------------------------------------------------
enable_testing()

# If reference NetCDF files are in the source dir, copy them to the build dir
# so the test can find them.
file(GLOB REF_NC_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/output-serial.nc"
    "${CMAKE_CURRENT_SOURCE_DIR}/output-serial-optimized.nc"
)
if(REF_NC_FILES)
    file(COPY ${REF_NC_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

# Test command closely mirrors the Makefile's `test` target:
#   1. Run the model with MPI
#   2. Parse Fractional Delta Mass / Energy
#   3. Check thresholds with awk
#   4. Run python3 nccmp3.py ...
#
# Note: We use bash -c and escape $ as \$ for CMake.
add_test(
    NAME model_test
    COMMAND bash -c
    "echo \"TEST: Running model...\"; \
     OUTPUT=\$(mpirun -np 1 ./model 100 1000 10); \
     MASS=\$(echo \"\${OUTPUT}\" | awk '/Fractional Delta Mass/{print \$NF}'); \
     ENERGY=\$(echo \"\${OUTPUT}\" | awk '/Fractional Delta Energy/{print \$NF}'); \
     awk -v m=\"\${MASS}\" -v e=\"\${ENERGY}\" 'BEGIN { \
       mass = m + 0; energy = e + 0; \
       if (mass < 1e-13 && energy < 1e-3) { \
         print \"TEST PASSED.\"; exit 0; \
       } else { \
         print \"TEST FAILED: thresholds exceeded.\"; \
         printf \"Fractional Delta Mass   = %e  (limit 1e-13)\\n\", mass; \
         printf \"Fractional Delta Energy = %e  (limit 1e-3)\\n\", energy; \
         exit 1; \
       } \
     }'; \
     echo \"TEST: Compare output with reference output file...\"; \
     python3 ${CMAKE_CURRENT_SOURCE_DIR}/nccmp3.py output-serial.nc output-serial-optimized.nc output.nc"
)

# Ensure test runs in the build directory (where ./model and output.nc live)
set_tests_properties(model_test PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)