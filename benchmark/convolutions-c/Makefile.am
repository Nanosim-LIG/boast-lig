EXTRA_DIST = \
	tool.h \
	thread_engine.h \
	topology_muscade.cfg \
	placement_muscade.cfg \
	thread_engine.cpp \
	$(SOURCES_CONV_CHECK)

SOURCES_CONV_CHECK =
SOURCES_SHARED =

if HAVE_LIBCONFIG
check = conv_check
else
check =
endif

EXTRA_PROGRAMS = conv_check #convolution_sse convolution_sse_t

CLEANFILES = fort.1 fort.98 *linkinfo

AM_FCFLAGS = -I../modules @LIBABINIT_INCLUDE@ @MPI_INCLUDE@ -I. -I$(srcdir) 
AM_CFLAGS = -Wall -Werror -O2 -ftree-vectorize -ftree-vectorizer-verbose=9 -ffast-math
#AM_CFLAGS = -Wall -Werror -O2 -mtune=core2 -ftree-vectorize -ftree-vectorizer-verbose=9 -ffast-math
#-funroll-loops -mtune=core2 -ftree-vectorize -floop-interchange -floop-strip-mine -floop-block

#convolution_sse_SOURCES = convolution_sse.c
#convolution_sse_LDADD = -lpapi

#convolution_sse_t_SOURCES = convolution_sse_t.c
#convolution_sse_t_LDADD = -lpapi


conv_check_SOURCES = conv_check.f90 MagicFilter.c #MagicFilter_sse.c
conv_check_LDADD = @LIBS_DEPENDENCIES@ -lpthread ~/papi/lib/libpapi.a

check_conv: conv_check
	echo 1 128 128 204900 204900 10 > fort.1
	./conv_check
