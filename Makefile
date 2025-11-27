#NVCC := nvcc
CUDA_PATH := /usr/local/cuda-12.3
NVCC := PATH=$(CUDA_PATH)/bin:$$PATH LD_LIBRARY_PATH=$(CUDA_PATH)/lib64:$$LD_LIBRARY_PATH nvcc

SRC_DIR := src
BIN_DIR := bin
OBJ_DIR := $(BIN_DIR)/obj

SRCS := $(wildcard $(SRC_DIR)/*.cu)
OBJS := $(patsubst $(SRC_DIR)/%.cu,$(OBJ_DIR)/%.o,$(SRCS))
EXEC := $(BIN_DIR)/upc_run

all: $(EXEC)

$(EXEC): $(OBJS)
	@echo "Linking $@"
	$(NVCC) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	@mkdir -p $(@D)
	@echo "Compiling $<"
	$(NVCC) -c $< -o $@

clean:
	rm -rf $(BIN_DIR)

.PHONY: all clean
