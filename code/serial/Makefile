# Compiler selection: Use NVHPC for OpenACC support
FC = mpif90

# OpenACC flag - set USE_OPENACC=1 to enable
USE_OPENACC ?= 0

# Base flags
FCFLAGS_BASE = -I. `nf-config --fflags`

ifeq ($(USE_OPENACC),1)
	# NVHPC compiler flags with OpenACC
    FC = mpifort
    FCFLAGS = -fast -acc -gpu=cc80 -Minfo=accel -Mpreprocess $(FCFLAGS_BASE) -DUSE_OPENACC
    LDFLAGS = -L`nc-config --libdir` `nf-config --flibs` -L$(CUDA_HOME)/lib64 -lnvToolsExt
else
	# GNU compiler flags (original)
	FCFLAGS = -Ofast $(FCFLAGS_BASE) -fPIC -fopenmp -DUSE_OPENMP -cpp
    LDFLAGS = -L`nc-config --libdir` `nf-config --flibs`
endif

all: model

model : model.F90 module_parameters.o module_types.F90 module_physics.o \
		module_output.o parallel_timer.o module_nvtx.o
	$(FC) $(FCFLAGS) $< -o $@ *.o $(LDFLAGS)

module_parameters.o : module_parameters.f90
	$(FC) $(FCFLAGS) $< -c

parallel_timer.o : parallel_timer.f90
	$(FC) $(FCFLAGS) $< -c

module_types.o : module_types.F90 module_parameters.o parallel_timer.o
	$(FC) $(FCFLAGS) $< -c

module_output.o : module_output.F90 module_physics.o
	$(FC) $(FCFLAGS) $< -c

module_physics.o : module_physics.f90 module_types.o parallel_timer.o
	$(FC) $(FCFLAGS) $< -c

module_nvtx.o: module_nvtx.F90
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

clean :
	rm -f *.o *.mod model *.nc

test: model
	@echo "TEST: Running model..."
	@OUTPUT=$$(mpirun -np 1 ./model 100 1000 10); \
# 	echo "$$OUTPUT"; \
	MASS=$$(echo "$$OUTPUT" | awk '/Fractional Delta Mass/{print $$NF}'); \
	ENERGY=$$(echo "$$OUTPUT" | awk '/Fractional Delta Energy/{print $$NF}'); \
# 	echo "Fractional Delta Mass   = $$MASS"; \
# 	echo "Fractional Delta Energy = $$ENERGY"; \
	awk -v m="$$MASS" -v e="$$ENERGY" 'BEGIN { \
	  mass = m + 0; energy = e + 0; \
	  if (mass < 1e-13 && energy < 1e-3) { \
		print "TEST PASSED."; exit 0; \
	  } else { \
		print "TEST FAILED: thresholds exceeded."; \
		printf "  mass   = %e  (limit 1e-13)\n", mass; \
		printf "  energy = %e  (limit 1e-3)\n", energy; \
		exit 1; \
	  } \
	}'
	@if [ -f output.nc ]; then \
	  echo "TEST: Compare output with reference output file..."; \
	  python3 nccmp3.py output-serial.nc output-serial-optimized.nc output.nc; \
	else \
	  echo "Output file not found, skipping comparison test."; \
	fi

