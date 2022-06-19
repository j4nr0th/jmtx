cmake_minimum_required(VERSION 3.10)

project(mtxlib)
set(CMAKE_C_STANDARD 17)
add_library(mtxlib matrix_base.h sparse_row_compressed.c sparse_row_compressed.h sparse_column_compressed.c sparse_column_compressed.h jacobi_point_iteration.c jacobi_point_iteration.h gauss_seidel_iteration.c gauss_seidel_iteration.h bicgstab_iteration.c bicgstab_iteration.h errors.c errors.h)
add_executable(test_mtx matrix_base.h sparse_row_compressed.h sparse_row_compressed.c common.c common.h sparse_column_compressed.c sparse_column_compressed.h jacobi_point_iteration.c jacobi_point_iteration.h iter_test.c gauss_seidel_iteration.c gauss_seidel_iteration.h bicgstab_iteration.c bicgstab_iteration.h errors.c errors.h)
target_link_libraries(test_mtx m pthread)
target_link_libraries(mtxlib pthread)

target_compile_definitions(test_mtx PRIVATE MTX_ERROR_MESSAGES MTX_MATRIX_CHECKS)
# target_compile_options(test_mtx PRIVATE -march=native -O3)