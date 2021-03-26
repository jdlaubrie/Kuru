PYTHON_VERSION = python2.7
PYTHON_LANG_LEVEL = 2
PYTHON_INCLUDE_PATH = /usr/include/python2.7/
PYTHON_LD_PATH = /usr/lib/python2.7/
EXT_POSTFIX = so
CYTHONFLAGS = -DNPY_NO_DEPRECATED_API

CYTHON = cython --cplus -$(PYTHON_LANG_LEVEL)

ifeq ($(OS),Windows_NT)
REMOVE = del
MOVE = move
else
REMOVE = rm -rf
MOVE = mv
endif
ENTER = cd

ADDITIONAL_FLAGS =

INLINE_THRESHOLD = 100000
CXXSTD=-std=c++14
OPTFLAGS= -O3 -fno-strict-aliasing -DNDEBUG
CXXFLAGS= -fPIC -shared -pthread -Wall $(CYTHONFLAGS) $(OPTFLAGS)

BLAS_VERSION = openblas
BLAS_INCLUDE_PATH = /opt/OpenBLAS/include/
BLAS_LD_PATH = /opt/OpenBLAS/lib/

ifeq ($(findstring mkl,$(BLAS_VERSION)),mkl)
CXXFLAGS += -DHAS_MKL
endif

#NUMPY_INCLUDE_PATH = /usr/local/lib/python2.7/dist-packages/numpy/core/include
NUMPY_INCLUDE_PATH = /usr/include/numpy/

ifeq ($(findstring g++,$(CXX)),g++)
ifneq ($(findstring clang++,$(CXX)),clang++)
CXXFLAGS += -Wno-cpp -Wno-unused-function -Wno-maybe-uninitialized -fabi-version=6
CXXFLAGS += -march=native -ffast-math -fwrapv -funroll-loops -finline-functions -finline-limit=$(INLINE_THRESHOLD) -ffp-contract=fast
else
CXXFLAGS += -Wno-uninitialized -Wno-missing-braces -Wno -warnings
CXXFLAGS += -march=native -ffast-math -mllvm -inline-threshold=$(INLINE_THRESHOLD) -ffp-contract=fast
endif
endif
ifeq ($(findstring c++,$(CXX)),c++)
CXXFLAGS += -Wno-uninitialized -Wno-missing-braces -Wno -NPY_NO_DEPRECATED_API
CXXFLAGS += -march=native -ffast-math -mllvm -inline-threshold=$(INLINE_THRESHOLD) -ffp-contract=fast
endif
ifeq ($(findstring icpc,$(CXX)),icpc)
CXXFLAGS += -march=native -no-ansi-alias -xHost -fp-model fast=1 -inline-forceinline -inline-factor=$(INLINE_THRESHOLD)
endif

ifeq ($(OS),Windows_NT)
CXXFLAGS += -D_hypot=hypot -m64 -Wno-format-extra-args -Wno-format
endif

FASTOR_INCLUDE_PATH = /usr/local/include/Fastor/
FASTOR_OPTFLAGS = -DFASTOR_NO_ALIAS -DFASTOR_USE_VECTORISED_EXPR_ASSIGN -DCONTRACT_OPT=2 -DFASTOR_NO_STATIC_WARNING
CXXFLAGS += $(FASTOR_OPTFLAGS) -I$(FASTOR_INCLUDE_PATH)

CXXFLAGS += $(ADDITIONAL_FLAGS)

all: build_module

build_module:
	@echo "Building low level mechanical constitutive stiffness integrands"
	$(CYTHON) DisplacementApproachIndices.pyx
	$(CXX) $(CXXSTD) $(CXXFLAGS) DisplacementApproachIndices.cpp -o DisplacementApproachIndices.$(EXT_POSTFIX) -I. -I$(PYTHON_INCLUDE_PATH) \
	     -I$(NUMPY_INCLUDE_PATH) -I$(BLAS_INCLUDE_PATH) -L$(PYTHON_LD_PATH) -l$(PYTHON_VERSION) -L$(BLAS_LD_PATH) -l$(BLAS_VERSION)
	@$(MOVE) DisplacementApproachIndices.$(EXT_POSTFIX) ../

	@echo "Building low level mechanical constitutive stiffness dispatcher modules"
	$(CYTHON) _ConstitutiveStiffnessDF_.pyx
	$(CXX) $(CXXSTD) $(CXXFLAGS) _ConstitutiveStiffnessDF_.cpp -o _ConstitutiveStiffnessDF_.$(EXT_POSTFIX) -I. -I$(PYTHON_INCLUDE_PATH) \
	     -I$(NUMPY_INCLUDE_PATH) -I$(BLAS_INCLUDE_PATH) -L$(PYTHON_LD_PATH) -l$(PYTHON_VERSION) -L$(BLAS_LD_PATH) -l$(BLAS_VERSION)
	@$(MOVE) _ConstitutiveStiffnessDF_.$(EXT_POSTFIX) ../

clean:
	$(REMOVE) DisplacementApproachIndices.cpp DisplacementApproachIndices.$(EXT_POSTFIX) ../DisplacementApproachIndices.$(EXT_POSTFIX)
	$(REMOVE) _ConstitutiveStiffnessDF_.cpp _ConstitutiveStiffnessDF_.$(EXT_POSTFIX) ../_ConstitutiveStiffnessDF_.$(EXT_POSTFIX)

source_clean:
	$(REMOVE) DisplacementApproachIndices.cpp
	$(REMOVE) _ConstitutiveStiffnessDF_.cpp
