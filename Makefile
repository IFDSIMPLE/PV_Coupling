BUILDS = compile build clean

.PHONY := $(BUILDS)

.DEFAULT_GOAL := compile

compile :
	@echo "Compiling Project"
	cmake --build build
	@echo "Done Compiling Project."

build :
	@echo "Building CMake Project"
	cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -S . -B build
	@echo "Done Building CMake Project"

clean :
	@echo "Cleaning CMake Project"
	rm -rf build bin
	@echo "Done Cleaning CMake Project"
	@echo "You need to run 'make build' to build the project again."

# Credit to Dr James Ferguson for providing us with this script via the Advanced Research Computing lectures