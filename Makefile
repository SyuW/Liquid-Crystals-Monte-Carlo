# final program
TARGET_EXEC := lcsim
TESTS_EXEC := tests

# directories
BUILD_DIR := ./build
SRC_DIRS := ./src/cpp

# find all C++ files for compiling the liquid crystal simulation executable
SRCS := $(shell find $(SRC_DIRS) -name '*.cpp' -not -path "./src/cpp/tests/*")

# object files we want to build
TARGET_OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

# compiler options
CPP := g++
CPPFLAGS := -I. $(INC_FLAGS) -Wall -MMD -MP

# build the final program
$(BUILD_DIR)/$(TARGET_EXEC): $(TARGET_OBJS)
	$(CPP) $(TARGET_OBJS) -o $@

# build object files from C++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CPP) $(CPPFLAGS) -c $< -o $@ 

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)