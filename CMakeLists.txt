cmake_minimum_required(VERSION 3.10)

project(jmtx)
set(CMAKE_C_STANDARD 23)

enable_testing()
add_library(jmtx
    include/jmtx/common.h
    source/matrix_base_internal.h
    source/matrix_base.c
    include/jmtx/matrix_base.h
    include/jmtx/solver_base.h
    )

include(source/float/CMakeLists.txt)
include(source/double/CMakeLists.txt)
include(source/cfloat/CMakeLists.txt)
include(source/cdouble/CMakeLists.txt)
include(source/conversion/CMakeLists.txt)



target_sources(jmtx PRIVATE
        ${MATRIX_TYPE_SOURCE_FILES} ${MATRIX_TYPE_HEADER_FILES}
        ${SOLVER_SOURCE_FILES} ${SOLVER_HEADER_FILES}
        ${INTERNAL_SOURCE_FILES} ${INTERNAL_HEADER_FILES}
        ${CONVERSION_SOURCE_FILES} ${CONVERSION_HEADER_FILES}
        )



if (CMAKE_BUILD_TYPE STREQUAL Debug)
    target_compile_options(jmtx PRIVATE -Wall -Wpedantic -Werror -Wextra)
elseif(CMAKE_BUILD_TYPE STREQUAL Release)
    #    target_compile_definitions(jmtx PRIVATE JMTX_NO_VERIFY_PARAMS=1)
endif ()

if (CMAKE_BUILD_TYPE STREQUAL Release-Generate-Profile)
    target_compile_options(jmtx PRIVATE -fprofile-use)
    target_link_libraries(jmtx gcov)
endif ()

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(jmtx gomp m)
    target_compile_options(jmtx PRIVATE -fopenmp -fanalyzer -march=native -mtune=native)
endif ()
