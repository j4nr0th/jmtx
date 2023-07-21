cmake_minimum_required(VERSION 3.10)

project(jmtx)
set(CMAKE_C_STANDARD 23)

list(APPEND MATRIX_TYPE_SOURCE_FILES
        source/matrices/sparse_row_compressed.c
        source/matrices/sparse_column_compressed.c
        source/matrices/dense_row_major.c
        source/matrices/dense_col_major.c
        source/matrices/matrix_base.c
        )

list(APPEND MATRIX_TYPE_HEADER_FILES
        source/matrices/sparse_row_compressed.h
        source/matrices/sparse_column_compressed.h
        source/matrices/dense_row_major.h
        source/matrices/dense_col_major.h
        source/matrices/matrix_base.h
        )

list(APPEND INTERNAL_HEADER_FILES
        source/matrices/sparse_row_compressed_internal.h)

list(APPEND INTERNAL_SOURCE_FILES)

list(APPEND SOLVER_SOURCE_FILES
        source/solvers/jacobi_point_iteration.c
        source/solvers/gauss_seidel_iteration.c
        source/solvers/bicgstab_iteration.c
        )

list(APPEND SOLVER_HEADER_FILES
        source/solvers/jacobi_point_iteration.h
        source/solvers/gauss_seidel_iteration.h
        source/solvers/bicgstab_iteration.h
        )

add_library(jmtx
        ${MATRIX_TYPE_SOURCE_FILES} ${MATRIX_TYPE_HEADER_FILES}
        ${SOLVER_SOURCE_FILES} ${SOLVER_HEADER_FILES}
        ${INTERNAL_SOURCE_FILES} ${INTERNAL_HEADER_FILES}
          )

#add_executable(test_mtx source/matrices/matrix_base.h source/matrices/sparse_row_compressed.h source/matrices/sparse_row_compressed.c tests/test_common.c tests/test_common.h source/matrices/sparse_column_compressed.c source/matrices/sparse_column_compressed.h source/solvers/jacobi_point_iteration.c source/solvers/jacobi_point_iteration.h iter_test.c source/solvers/gauss_seidel_iteration.c source/solvers/gauss_seidel_iteration.h source/solvers/bicgstab_iteration.c source/solvers/bicgstab_iteration.h tests/base_crs_test.c)
#target_link_libraries(test_mtx m pthread)
target_link_libraries(jmtx pthread jdm)

#target_compile_definitions(test_mtx PRIVATE MTX_ERROR_MESSAGES MTX_MATRIX_CHECKS)
target_compile_options(jmtx PRIVATE -march=native)

enable_testing()

add_executable(base_crs_test tests/base_crs_test.c tests/test_common.c)
target_link_libraries(base_crs_test jmtx m)
add_test(NAME base_crs COMMAND base_crs_test)

add_executable(base_ccs_test tests/base_ccs_test.c tests/test_common.c)
target_link_libraries(base_ccs_test jmtx m)
add_test(NAME base_ccs COMMAND base_ccs_test)

add_executable(old_crs_demo run_crs_tests.c)
target_link_libraries(old_crs_demo PRIVATE jmtx m)
add_test(NAME old_crs COMMAND old_crs_demo)

add_executable(old_ccs_demo run_ccs_tests.c)
target_link_libraries(old_ccs_demo PRIVATE jmtx m)
add_test(NAME old_ccs COMMAND old_ccs_demo)

add_executable(old_compare_demo compare_crs_and_ccs.c tests/test_common.c)
target_link_libraries(old_compare_demo PRIVATE jmtx m)
add_test(NAME old_compare COMMAND old_compare_demo)

add_executable(add_crs_test tests/crs_addition_test.c tests/test_common.c)
target_link_libraries(add_crs_test PRIVATE jmtx m)
add_test(NAME addition_crs COMMAND add_crs_test)

#add_subdirectory(cplot)

add_executable(jacobi_point_test tests/jacobi_test.c tests/test_common.c source/solvers/jacobi_point_iteration.c)
#target_include_directories(jacobi_point_test PRIVATE cplot/include)
target_link_libraries(jacobi_point_test PRIVATE jmtx m
        #cplot
        )
add_test(NAME jacobi COMMAND jacobi_point_test)

if (CMAKE_BUILD_TYPE STREQUAL Debug)
    target_compile_options(jmtx PRIVATE -Wall -Wpedantic -Werror -Wextra)
elseif(CMAKE_BUILD_TYPE STREQUAL Release)
#    target_compile_definitions(jmtx PRIVATE JMTX_NO_VERIFY_PARAMS=1)
endif ()
