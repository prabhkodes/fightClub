FC = mpif90
#FCFLAGS = -O0 -Wall -pedantic -fcheck=all -fbacktrace -I. `nf-config --fflags`
FCFLAGS = -Ofast -I. `nf-config --fflags` -fPIC
LDFLAGS = -L`nc-config --libdir` `nf-config --flibs`

all: model

model : model.f90 module_parameters.o module_types.F90 module_physics.o \
        module_output.o parallel_timer.o
	$(FC) $(FCFLAGS) $< -o $@ *.o $(LDFLAGS)

module_parameters.o : module_parameters.f90
	$(FC) $(FCFLAGS) $< -c

module_types.o : module_types.F90 module_parameters.o
	$(FC) $(FCFLAGS) $< -c

module_output.o : module_output.F90 module_physics.o
	$(FC) $(FCFLAGS) $< -c

module_physics.o : module_physics.f90 module_types.o parallel_timer.o
	$(FC) $(FCFLAGS) $< -c

parallel_timer.o : parallel_timer.f90
	$(FC) $(FCFLAGS) $< -c

clean :
	rm -f *.o *.mod model *.nc


test: model
	@echo "Running model test..."
	@OUTPUT=$$(./model); \
	echo "$$OUTPUT"; \
	MASS=$$(echo "$$OUTPUT" | awk '/Fractional Delta Mass/{print $$NF}'); \
	ENERGY=$$(echo "$$OUTPUT" | awk '/Fractional Delta Energy/{print $$NF}'); \
	echo "Parsed Mass   = $$MASS"; \
	echo "Parsed Energy = $$ENERGY"; \
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
