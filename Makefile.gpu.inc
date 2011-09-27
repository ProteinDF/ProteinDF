CUDA_INSTALL_PATH = /usr/local/cuda/
CXXFLAGS = -I${HOME}/NVIDIA_CUDA_SDK/common/inc -I../../include -I../pdf -DHAVE_CONFIG_H -DHAVE_CUDA

LIBS = ${HOME}/NVIDIA_CUDA_SDK/lib/libcutil.a

NVCC := $(CUDA_INSTALL_PATH)/bin/nvcc
CXX := g++
LINK := g++ -fPIC

# Debug/release configuration
ifeq ($(dbg),1)
        COMMONFLAGS += -g
        NVCCFLAGS   += -D_DEBUG --host-compilation c++
        BINSUBDIR   := debug
        LIBSUFFIX   := D
else
        COMMONFLAGS += -O3
        BINSUBDIR   := release
        LIBSUFFIX   :=
        NVCCFLAGS   += --host-compilation c++
        CXXFLAGS    +=
        CFLAGS      +=
endif

.SUFFIXES : .cu .cu_o .cu_dbg.o .c_dbg.o .cpp_dbg.o .cu_rel.o .c_rel.o .cpp_rel.o .cubin .cpp .cxx
.PHONY : clean

.cpp.o : $(SRCDIR)%.cpp $(C_DEPS)
	$(CXX) $(CXXFLAGS) -c $<

.cxx.o : $(SRCDIR)%.cpp $(C_DEPS)
	$(CXX) $(CXXFLAGS) -c $<

.cu.cu_o : $(SRCDIR)%.cu $(CU_DEPS)
	#$(NVCC) $(NVCCFLAGS) --compiler-options "$(CXXFLAGS)" $(SMVERSIONFLAGS) -c $<
	$(NVCC) $(NVCCFLAGS) --compiler-options "$(CXXFLAGS)" $(SMVERSIONFLAGS) -c $< -o $*.cu_o

%.cubin : $(SRCDIR)%.cu cubindirectory
	$(NVCC) $(CUBIN_ARCH_FLAG) $(NVCCFLAGS) $(SMVERSIONFLAGS) -cubin $<

