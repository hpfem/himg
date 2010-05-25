#
# LIBPNG
#
if(MSVC)
   set(LIBPNG_NAME libpng)
else(MSVC)
   set(LIBPNG_NAME png)
endif(MSVC)

FIND_PATH(LIBPNG_INCLUDE_DIR png.h ${LIBPNG_ROOT}/include /usr/include /usr/include/umfpack /usr/local/include/UMFPACK /usr/include/suitesparse)
FIND_LIBRARY(LIBPNG_LIBRARY ${LIBPNG_NAME} ${LIBPNG_ROOT}/lib /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBPNG DEFAULT_MSG LIBPNG_LIBRARY LIBPNG_INCLUDE_DIR)
