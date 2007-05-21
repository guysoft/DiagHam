SUBDIRS=

noinst_LIBRARIES=libQHEWaveFunction.a libMC.a

libQHEWaveFunction_a_SOURCES=LaughlinOnSphereWaveFunction.cc \
PfaffianOnSphereWaveFunction.cc \
PfaffianOnSphereTwoQuasiholeWaveFunction.cc \
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
SumDerivativeProduct.cc \
AbstractFQHEWaveFunctionOneOverR.cc \
LaughlinOnDiskWaveFunctionOneOverR.cc \
HalperinOnDiskWaveFunctionOneOverR.cc \
HalperinOnSphereWaveFunction.cc \
SU4HalperinOnSphereWaveFunction.cc \
PairedCFOnSphereWaveFunction.cc 

libMC_a_SOURCES = \
CFOrbitals.cc \
formula_h6.cc  \
my_ran.cc \
my_recipes.cc \
nrutil.cc \
particle.cc \
OrbitalFactoryNeg.cc \
formula_neg.cc  \
ReferenceWaveFunction.cc \
pfaffian.cc \
cSymMat.cc 


bin_PROGRAMS=test test_neg

test_SOURCES=test.cc

test_DEPENDENCIES=libQHEWaveFunction.a libMC.a $(top_builddir)/src/Matrix/libArray.a

test_LDADD=$(top_builddir)/src/GeneralTools/libGeneralTools.a  $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/MathTools/RandomNumber/libRandomNumber.a $(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/GeneralTools/libGeneralTools.a libQHEWaveFunction.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/Matrix/libArray.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a  libMC.a @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@

test_neg_SOURCES=test_neg.cc

test_neg_DEPENDENCIES=libQHEWaveFunction.a libMC.a $(top_builddir)/src/Matrix/libArray.a

test_neg_LDADD=$(top_builddir)/src/GeneralTools/libGeneralTools.a  $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a  $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/MathTools/RandomNumber/libRandomNumber.a $(top_builddir)/src/Vector/libVector.a $(top_builddir)/src/Output/libOutput.a $(top_builddir)/src/GeneralTools/libGeneralTools.a libQHEWaveFunction.a $(top_builddir)/src/Matrix/libArray.a $(top_builddir)/src/Matrix/libMatrix.a $(top_builddir)/src/MathTools/libMathTools.a $(top_builddir)/src/MathTools/NumericalAnalysis/libNumericalAnalysis.a libMC.a @MATH_LIBS@ @THREAD_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @FORTRAN2C_LIBS@
