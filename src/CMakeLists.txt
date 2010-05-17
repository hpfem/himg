set (SRC
    himg.cpp pix_by_pix_integrator.cpp h1_adapt_image.cpp h1_image_proj_based_selector.cpp
	shared/bitmap.cpp shared/img2d.cpp
	adapt_storage.cpp
	parameters.cpp shared/input_params.cpp)

# source codes
add_executable(${PROJECT_NAME} ${SRC})

# definitions
# if(NOT WITH_JPEG)
	add_definitions("-DHIMG_NO_JPEG")
# endif(NOT WITH_JPEG)

# include directories
include_directories(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
include_directories(${HERMES2D_INCLUDE_DIR})
include_directories(${LIBPNG_INCLUDE_DIR})

# linking
target_link_libraries(${PROJECT_NAME} ${HERMES2D_LIBRARY} ${UMFPACK_LIBRARY} ${AMD_LIBRARY} ${LIBPNG_LIBRARY})
