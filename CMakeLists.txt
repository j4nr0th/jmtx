cmake_minimum_required(VERSION 3.10)

project(jmtx)
set(CMAKE_C_STANDARD 99)

find_package(OpenMP REQUIRED)
if (CMAKE_PROJECT_NAME STREQUAL "jmtx")
    enable_testing()
endif ()
add_library(jmtx
        source/common.h
        source/matrix_base_internal.h
        source/matrix_base.c
        source/matrix_base.h
        source/solver_base.h
        source/modules/float_module.c
        source/modules/double_module.c
        source/modules/cfloat_module.c
        source/modules/cdouble_module.c
)

message("Added extra sources: ${MATRIX_MODULE_SOURCE_FILES}")

set_property(TARGET jmtx PROPERTY POSITION_INDEPENDENT_CODE ON)

include(CheckIPOSupported)
check_ipo_supported(RESULT result LANGUAGES C)
if (result)
    set_property(TARGET jmtx PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

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

enable_testing()


set(TEST_INCLUDE_HEADER "../source/modules/jmtx_float_defines.h")
set(TEST_TYPE_NAME "float")
add_subdirectory(tests ${CMAKE_BINARY_DIR}/float-tests)

set(TEST_INCLUDE_HEADER "../source/modules/jmtx_cfloat_defines.h")
set(TEST_TYPE_NAME "cfloat")
add_subdirectory(tests ${CMAKE_BINARY_DIR}/cfloat-tests)
