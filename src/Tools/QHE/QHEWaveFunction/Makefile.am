SUBDIRS=

noinst_LIBRARIES=libQHEWaveFunction.a

libQHEWaveFunction_a_SOURCES=LaughlinOnSphereWaveFunction.cc \
PfaffianOnSphereWaveFunction.cc \
JainCFFilledLevelOnSphereWaveFunction.cc \
JainCFOnSphereOrbitals.cc \
MooreReadOnSphereWaveFunction.cc \
LaughlinOnDiskWaveFunction.cc \
PfaffianOnDiskWaveFunction.cc \
MooreReadOnDiskWaveFunction.cc \
JainCFOnSphereWaveFunction.cc \
LandauSpectrumOnSphere.cc \
UnprojectedJainCFOnSphereWaveFunction.cc \
QHEWaveFunctionManager.cc \
ExactWaveFunction.cc \
DerivativeProduct.cc \
DerivativeProductFactor.cc \
SumDerivativeProduct.cc

bin_PROGRAMS=test

test_SOURCES=test.cc

test_DEPENDENCIES=libQHEWaveFunction.a

test_LDADD=$(top_builddir)/src/GeneralTools/libGeneralTools.a  $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/MathTools/RandomNumber/libRandomNumber.a $(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/GeneralTools/libGeneralTools.a libQHEWaveFunction.a $(top_builddir)/src/MathTools/libMathTools.a  @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@
