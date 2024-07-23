#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([examesh],[1.0],[cfog@mech.ubc.ca])
AC_CONFIG_AUX_DIR([conf])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([refine.cxx])
AC_CONFIG_HEADERS([exa_config.h])

# Checks for programs.
AX_PROG_CXX_MPI
AX_PROG_CC_MPI

LT_INIT

dnl m4_include([conf/m4_ax_openmp.m4])

dnl AC_LANG_PUSH(C++)
dnl AX_OPENMP
dnl AC_SUBST(OPENMP_CXXFLAGS)
dnl AC_LANG_POP(C++)

m4_include([m4/boost.m4])

BOOST_TEST
BOOST_MPI
BOOST_SERIALIZATION

dnl This is only a start on doing configuration for boost better


AC_CHECK_HEADERS(boost/mpi.hpp)
AC_CHECK_LIB([boost_mpi], [main])

# Checks for header files.
AC_CHECK_HEADERS([limits.h locale.h stdint.h stdlib.h string.h unistd.h values.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT32_T

# Checks for library functions.
AC_CHECK_FUNCS([setlocale])
AC_CHECK_LIB(m, sqrt)

################################################################################
#                              CGNS Options 
################################################################################
# Add --with-CGNS-path option to configure script
# Default
AC_ARG_WITH( CGNS-path,
       [AS_HELP_STRING([--with-CGNS-path=PATH],[specify installed location for the CGNS file format library])],,
       [with_CGNS_path="check"] )

CGNS_INCLUDES=""
CGNS_LDFLAGS=""
CGNS_LIBS=""

if (test "x$with_CGNS_path" != "xcheck") ; then
	AC_LANG_PUSH(C++)
#	CPPFLAGS_tmp=$CPPFLAGS
#	CPPFLAGS=$CPPFLAGS -I$with_CGNS_path/include"
	AC_CHECK_HEADER( cgnslib.h,
		[ AC_DEFINE([HAVE_CGNS],[1],["Have CGNS headers"])
		CGNS_INCLUDES=-I$with_CGNS_path/include ],
		[CPPFLAGS=$CPPFLAGS_tmp]
		)
	AC_LANG_POP(C++)
	
	LDFLAGS_tmp=$LDFLAGS
	LDFLAGS="$LDFLAGS -L$with_CGNS_path/lib"
    CGNS_LDFLAGS="-L$with_CGNS_path/lib"
    AC_CHECK_LIB(cgns, cg_close, 
            [AC_DEFINE([HAVE_CGNS_LIB], [1],
               [Define if you find the CGNS library.])
     CGNS_LIBS=-lcgns
     CGNS_LDFLAGS=-L$with_CGNS_path/lib
     ])
	LDFLAGS=$LDFLAGS_tmp
else
	AC_LANG_PUSH(C++)
	AC_CHECK_HEADER( cgnslib.h,
		[ AC_DEFINE([HAVE_CGNS],[1],["Have CGNS headers"])]
		)
	AC_LANG_POP(C++)
	
	AC_CHECK_LIB(cgns, cg_close, 
            [AC_DEFINE([HAVE_CGNS_LIB], [1],
               [Define if you find the CGNS library.])
        CGNS_LIBS=-lcgns
     ])
fi

AC_SUBST(CGNS_INCLUDES)
AC_SUBST(CGNS_LDFLAGS)
AC_SUBST(CGNS_LIBS)
	

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lMeshIO':
AC_ARG_WITH( MeshIO-path,
	     [AS_HELP_STRING([--with-MeshIO-path=PATH],[specify installed location for the MeshIO library])],
		[CPPFLAGS="$CPPFLAGS -I$with_MeshIO_path"])
AC_CHECK_LIB([MeshIO], [libMeshIO_is_present],
		       [LDFLAGS="-Wl,-rpath -Wl,$with_MeshIO_path -L$with_MeshIO_path $LDFLAGS"
		       LIBS="-lMeshIO $LIBS"],,
		       ["-L$with_MeshIO_path"])


AC_CHECK_LIB([sz], [SZ_encoder_enabled])

AC_ARG_WITH( HDF5-path,
	     [AS_HELP_STRING([--with-HDF5-path=PATH],[specify installed location for the HDF5 library])],,
	     [with_HDF5_path="/usr/include/hdf5/serial"] )

CPPFLAGS_tmp=$CPPFLAGS
CPPFLAGS="$CPPFLAGS -I$with_HDF5_path"
AC_CHECK_HEADER(hdf5.h,
	[ AC_DEFINE([HAVE_HDF5],[1],["Have HDF5 headers"])])
CGNS_INCLUDES=-I$with_HDF5_path
AC_CHECK_LIB([hdf5_serial], [H5Ecreate_stack],
			    [LDFLAGS="-L$with_HDF5_path $LDFLAGS"
			     LIBS="-lhdf5_serial $LIBS"],,
			    [-L$with_HDF5_path])
CPPFLAGS=$CPPFLAGS_tmp

AC_SUBST(HAVE_LIBHDF5_SERIAL)

AC_ARG_WITH( CGNS-path,
	     [AS_HELP_STRING([--with-CGNS-path=PATH],[specify installed location for the CGNS file format library])],,
	     [with_CGNS_path="no"] )

if (test "x$with_CGNS_path" == "xno") ; then
   AC_CHECK_HEADER(cgnslib.h,
	[ AC_DEFINE([HAVE_CGNS],[1],["Have CGNS headers"])]
	)
   AC_CHECK_LIB([cgns], [cg_set_file_type])
else 
   CPPFLAGS_tmp=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS -I$with_CGNS_path/include"
   AC_CHECK_HEADER(cgnslib.h,
	[ AC_DEFINE([HAVE_CGNS],[1],["Have CGNS headers"])])
   CGNS_INCLUDES+=" -I$with_CGNS_path/include"
   AC_CHECK_LIB([cgns], [cg_set_file_type],[LDFLAGS="-L$with_CGNS_path/lib $LDFLAGS"
			     LIBS="-lcgns $LIBS"],,[-L$with_CGNS_path/lib])
   CPPFLAGS=$CPPFLAGS_tmp
fi
CPPFLAGS+=" $CGNS_INCLUDES"

AC_SUBST(HAVE_LIBCGNS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
