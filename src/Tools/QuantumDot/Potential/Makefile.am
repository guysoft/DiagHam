SUBDIRS=

noinst_LIBRARIES=libPotential.a

libPotential_a_SOURCES=Potential.cc \
TwoDPotential.cc \
ThreeDPotential.cc \
SegregationThreeDPotential.cc \
PeriodicPyramidQuantumDot.cc

check_PROGRAMS=TestPotential

TESTS=

TestPotential_SOURCES=TestPotential.cc

TestPotential_DEPENDENCIES=@top_srcdir@/src/Tools/QuantumDot/Potential/libPotential.a

TestPotential_LDADD=-L@top_srcdir@/src/Tools/QuantumDot/Potential -lPotential @MATH_LIBS@ @THREAD_LIBS@
