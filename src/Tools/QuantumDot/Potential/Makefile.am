SUBDIRS=

noinst_LIBRARIES=libPotential.a

libPotential_a_SOURCES=AbstractPotential.cc \
OneDConstantCellPotential.cc  \
TwoDConstantCellPotential.cc  \
BinaryTwoDConstantCellPotential.cc  \
ThreeDConstantCellPotential.cc  \
BinaryThreeDConstantCellPotential.cc  \
HardBoxPyramidQuantumDotThreeDConstantCellPotential.cc  \
PeriodicPyramidQuantumDotThreeDConstantCellPotential.cc \
DotEmbeddedWellThreeDConstantCellPotential.cc

check_PROGRAMS=TestPotential

TESTS=

TestPotential_SOURCES=TestPotential.cc

TestPotential_DEPENDENCIES=@top_srcdir@/src/Tools/QuantumDot/Potential/libPotential.a

TestPotential_LDADD=-L@top_srcdir@/src/Tools/QuantumDot/Potential -lPotential @MATH_LIBS@ @THREAD_LIBS@
