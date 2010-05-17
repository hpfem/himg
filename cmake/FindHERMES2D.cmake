#
# HERMES2D
#

FIND_PATH(HERMES2D_INCLUDE_DIR hermes2d.h ${HERMES2D_SOURCE_ROOT} ${HERMES2D_ROOT}/include /usr/include /usr/include/umfpack /usr/local/include/UMFPACK /usr/include/suitesparse)
FIND_LIBRARY(HERMES2D_LIBRARY_RELEASE hermes2d-real ${HERMES2D_ROOT}/lib /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib)
FIND_LIBRARY(HERMES2D_LIBRARY_DEBUG hermes2d-real-debug ${HERMES2D_ROOT}/lib /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib)

set(HERMES2D_LIBRARY "${HERMES2D_LIBRARY_RELEASE}")
if(DEBUG)
    set(HERMES2D_LIBRARY "${HERMES2D_LIBRARY_DEBUG}")
endif(DEBUG)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HERMES2D DEFAULT_MSG HERMES2D_LIBRARY HERMES2D_INCLUDE_DIR)
