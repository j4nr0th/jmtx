cmake_minimum_required(VERSION 3.10)

project(jmtx)
set(CMAKE_C_STANDARD 99)

find_package(OpenMP REQUIRED)
if (CMAKE_PROJECT_NAME STREQUAL "jmtx")
    enable_testing()
endif()
add_library(jmtx
    include/jmtx/common.h
    source/matrix_base_internal.h
    source/matrix_base.c
    include/jmtx/matrix_base.h
    include/jmtx/solver_base.h
)

include(source/float/CMakeLists.txt)
include(source/double/CMakeLists.txt)

if(NOT CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    include(source/cfloat/CMakeLists.txt)
    include(source/cdouble/CMakeLists.txt)
endif()
include(source/conversion/CMakeLists.txt)



target_sources(jmtx PRIVATE
        ${MATRIX_TYPE_SOURCE_FILES} ${MATRIX_TYPE_HEADER_FILES}
        ${SOLVER_SOURCE_FILES} ${SOLVER_HEADER_FILES}
        ${INTERNAL_SOURCE_FILES} ${INTERNAL_HEADER_FILES}
        ${CONVERSION_SOURCE_FILES} ${CONVERSION_HEADER_FILES}
        ${DECOMPOSITIONS_SOURCE_FILES} ${DECOMPOSITIONS_HEADER_FILES}
        )

set_property(TARGET jmtx PROPERTY POSITION_INDEPENDENT_CODE ON)

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_compile_options(jmtx PRIVATE -Wall -Wpedantic -Werror -Wextra)
    if (CMAKE_BUILD_TYPE STREQUAL Debug)
        target_compile_options(jmtx PRIVATE -Wall -Wpedantic -Werror -Wextra -fanalyzer)
    endif ()
endif ()
#   MSVC support for OpenMP is only for C89, so don't include it for MSVC, since the project is written for C99
if (OpenMP_C_FOUND AND NOT CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    target_link_libraries(jmtx OpenMP::OpenMP_C ${OpenMP_C_LIBRARIES})
    target_compile_options(jmtx PRIVATE ${OpenMP_C_FLAGS})
    target_include_directories(jmtx PRIVATE OpenMP_C_INCLUDE_DIRS)
endif ()

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(jmtx m)
#elseif(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
#    add_compile_definitions(jmtx JMTX_NO_ARRAY_ATTRIBS)
#    add_compile_definitions(jmtx JMTX_MSVC)
endif ()
