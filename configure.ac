#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([examesh], [1.0], [cfog@mech.ubc.ca])
AC_CONFIG_AUX_DIR([conf])
AC_CONFIG_MACRO_DIR([conf])
AC_CONFIG_SRCDIR([ExaMesh.h])
AC_CONFIG_HEADERS([exa_config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

m4_include([conf/m4_ax_openmp.m4])

AC_LANG_PUSH(C++)
AX_OPENMP
AC_SUBST(OPENMP_CXXFLAGS)
AC_LANG_POP(C++)

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

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lMeshIO':
AC_ARG_WITH( MeshIO-path,
	     [AC_HELP_STRING([--with-MeshIO-path=PATH],
		[specify installed location for the MeshIO library])],
		[CPPFLAGS="$CPPFLAGS -I$with_MeshIO_path"])
AC_CHECK_LIB([MeshIO], [libMeshIO_is_present],
		       [LDFLAGS="-Wl,-rpath -Wl,$with_MeshIO_path -L$with_MeshIO_path $LDFLAGS"
		       LIBS="-lMeshIO $LIBS"],,
		       ["-L$with_MeshIO_path"])


# FIXME: Replace `main' with a function in `-lboost_unit_test_framework':
AC_CHECK_LIB([boost_unit_test_framework], [main])

AC_CHECK_LIB([sz], [SZ_encoder_enabled])

AC_ARG_WITH( HDF5-path,
	     [AC_HELP_STRING([--with-HDF5-path=PATH],[specify installed location for the HDF5 library])],,
	     [with_HDF5_path="/usr/include/hdf5/serial"] )

CPPFLAGS_tmp=$CPPFLAGS
CPPFLAGS="$CPPFLAGS -I$with_HDF5_path"
AC_CHECK_HEADER(hdf5.h,
	[ AC_DEFINE([HAVE_HDF5],[1],["Have HDF5 headers"])],
	CGNS_INCLUDES=-I$with_HDF5_path
)
AC_CHECK_LIB([hdf5_serial], [H5Ecreate_stack],
			    [LDFLAGS="-L$with_HDF5_path $LDFLAGS"
			     LIBS="-lhdf5_serial $LIBS"],,
			    [-L$with_HDF5_path])
CPPFLAGS=$CPPFLAGS_tmp

AC_SUBST(HAVE_LIBHDF5_SERIAL)

AC_ARG_WITH( CGNS-path,
	     [AC_HELP_STRING([--with-CGNS-path=PATH],[specify installed location for the CGNS file format library])],,
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
	[ AC_DEFINE([HAVE_CGNS],[1],["Have CGNS headers"])],
	CGNS_INCLUDES+=-I$with_CGNS_path/include
	)
   AC_CHECK_LIB([cgns], [cg_set_file_type],,,[$with_CGNS_path])
   CPPFLAGS=$CPPFLAGS_tmp
fi
CPPFLAGS+=$CGNS_INCLUDES

AC_SUBST(HAVE_LIBCGNS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
