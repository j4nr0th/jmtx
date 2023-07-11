cmake_minimum_required(VERSION 3.10)

project(mtxlib)
set(CMAKE_C_STANDARD 23)

list(APPEND MATRIX_TYPE_SOURCE_FILES
        source/matrices/sparse_row_compressed.c
        source/matrices/sparse_column_compressed.c
        # source/matrices/dense_matrix.c
        source/matrices/matrix_base.c
        )

list(APPEND MATRIX_TYPE_HEADER_FILES
        source/matrices/sparse_row_compressed.h
        source/matrices/sparse_column_compressed.h
        # source/matrices/dense_matrix.h
        source/matrices/matrix_base.h
        )

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

add_library(mtxlib
        ${MATRIX_TYPE_SOURCE_FILES} ${MATRIX_TYPE_HEADER_FILES}
        ${${SOLVER_SOURCE_FILE}} ${SOLVER_HEADER_FILE})

add_executable(test_mtx source/matrices/matrix_base.h source/matrices/sparse_row_compressed.h source/matrices/sparse_row_compressed.c common.c common.h source/matrices/sparse_column_compressed.c source/matrices/sparse_column_compressed.h source/solvers/jacobi_point_iteration.c source/solvers/jacobi_point_iteration.h iter_test.c source/solvers/gauss_seidel_iteration.c source/solvers/gauss_seidel_iteration.h source/solvers/bicgstab_iteration.c source/solvers/bicgstab_iteration.h)
target_link_libraries(test_mtx m pthread)
target_link_libraries(mtxlib pthread jdm)

#target_compile_definitions(test_mtx PRIVATE MTX_ERROR_MESSAGES MTX_MATRIX_CHECKS)
target_compile_options(mtxlib PRIVATE -march=native)

if (CMAKE_BUILD_TYPE STREQUAL Debug)
    target_compile_options(mtxlib PRIVATE -Wall -Wpedantic -Werror -Wextra)
endif ()
