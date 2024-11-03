MPICC = mpicc
CFLAGS_MPI = -fopenmp -O3 -Wall
LINK_FLAGS = -lm

SRC_DIR = src
OBJ_DIR = obj

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

# Executable names
MAIN_EXEC = simulation.out

# Main target that builds both executables
all: $(MAIN_EXEC)

# Build with Score-P profiling enabled
profile: CFLAGS_MPI += -DUSE_SCOREP -g
profile: MPICC := scorep $(MPICC)
profile: $(MAIN_EXEC)

# Build the main test executable with mpicc
$(MAIN_EXEC): $(OBJS)
	$(MPICC) $(CFLAGS_MPI) $(OBJS) -o $(MAIN_EXEC) $(LINK_FLAGS)

# Build object files for the source files in the src directory with mpicc
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(MPICC) $(CFLAGS_MPI) -c $< -o $@

# Clean command that removes object files and both executables
clean:
	rm -rf $(OBJ_DIR) $(MAIN_EXEC)
