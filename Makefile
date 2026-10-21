FC := gfortran

PROGRAM := code
BUILD_DIR := build
MOD_DIR := $(BUILD_DIR)/mod

# OpenMP enables the parallel regions already present in the source code.
# Both free-form and fixed-form source line-length limits are disabled.
FFLAGS := -O3 -fopenmp -ffree-line-length-none -ffixed-line-length-none \
	-std=legacy -fallow-argument-mismatch -J$(MOD_DIR) -I$(MOD_DIR)
LDFLAGS := -fopenmp

# model.f90 already contains the parameter modules.  parameters.f90 is an
# alternative parameter set and must not be linked at the same time, because
# it defines the same module names.
SOURCES := model.f90 \
	operations.f90 \
	calculations.f90 \
	integration_method.f90 \
	initial_conditions.f90 \
	source_profile.f90 \
	task.f90 \
	main_program.f90

OBJECTS := $(SOURCES:%.f90=$(BUILD_DIR)/%.o)

.PHONY: all clean run

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(FC) $(OBJECTS) $(LDFLAGS) -o $@

# model.f90 creates all project module files, so compile it first.  Once it is
# ready, the remaining source files can be compiled concurrently with make -j.
$(BUILD_DIR)/model.o: model.f90 | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(filter-out $(BUILD_DIR)/model.o,$(OBJECTS)): $(BUILD_DIR)/model.o

$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR) $(MOD_DIR):
	mkdir -p $@

run: $(PROGRAM)
	./$(PROGRAM)

clean:
	rm -rf $(BUILD_DIR) $(PROGRAM)
