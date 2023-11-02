cmake_minimum_required(VERSION 3.10)

project(jmtx)
enable_language(Fortran)
set(CMAKE_C_STANDARD 23)
set(CMAKE_Fortran_STANDARD 95)

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

list(APPEND FORTRAN_SOURCE_FILES
        source/f_exper/f_crs_base.f90
        source/f_exper/fjacobi.f90)

add_library(jmtx
        ${MATRIX_TYPE_SOURCE_FILES} ${MATRIX_TYPE_HEADER_FILES}
        ${SOLVER_SOURCE_FILES} ${SOLVER_HEADER_FILES}
        ${INTERNAL_SOURCE_FILES} ${INTERNAL_HEADER_FILES}
        ${FORTRAN_SOURCE_FILES}
        source/solvers/conjugate_gradient_iteration.c
        source/solvers/conjugate_gradient_iteration.h
        source/solvers/incomplete_lu_decomposition.c
        source/solvers/incomplete_lu_decomposition.h
)

#add_executable(test_mtx source/matrices/matrix_base.h source/matrices/sparse_row_compressed.h source/matrices/sparse_row_compressed.c tests/test_common.c tests/test_common.h source/matrices/sparse_column_compressed.c source/matrices/sparse_column_compressed.h source/solvers/jacobi_point_iteration.c source/solvers/jacobi_point_iteration.h iter_test.c source/solvers/gauss_seidel_iteration.c source/solvers/gauss_seidel_iteration.h source/solvers/bicgstab_iteration.c source/solvers/bicgstab_iteration.h tests/base_crs_test.c)
#target_link_libraries(test_mtx m pthread)
target_link_libraries(jmtx pthread jdm m)

#target_compile_definitions(test_mtx PRIVATE MTX_ERROR_MESSAGES MTX_MATRIX_CHECKS)
target_compile_options(jmtx PRIVATE -march=native)

option(JMTX_USE_OPENMP "Use OpenMP for parallel functions" ON)

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_compile_options(jmtx PRIVATE -fopenmp)
endif ()
if (JMTX_USE_OPENMP)
    target_compile_definitions(jmtx PRIVATE JMTX_OPENMP)
endif ()

enable_testing()

add_executable(base_crs_test tests/matrix_ops_tests/base_crs_test.c tests/test_common.c)
target_link_libraries(base_crs_test jmtx m)
add_test(NAME base_crs COMMAND base_crs_test)

add_executable(base_ccs_test tests/matrix_ops_tests/base_ccs_test.c tests/test_common.c)
target_link_libraries(base_ccs_test jmtx m)
add_test(NAME base_ccs COMMAND base_ccs_test)

add_executable(old_crs_demo run_crs_tests.c tests/test_common.c)
target_link_libraries(old_crs_demo PRIVATE jmtx m)
add_test(NAME old_crs COMMAND old_crs_demo)

add_executable(old_ccs_demo run_ccs_tests.c tests/test_common.c)
target_link_libraries(old_ccs_demo PRIVATE jmtx m)
add_test(NAME old_ccs COMMAND old_ccs_demo)

add_executable(old_compare_demo compare_crs_and_ccs.c tests/test_common.c)
target_link_libraries(old_compare_demo PRIVATE jmtx m)
add_test(NAME old_compare COMMAND old_compare_demo)

add_executable(add_crs_test tests/matrix_ops_tests/crs_addition_test.c tests/test_common.c)
target_link_libraries(add_crs_test PRIVATE jmtx m)
add_test(NAME addition_crs COMMAND add_crs_test)

#add_subdirectory(cplot)

add_executable(jacobi_point_test tests/solver_tests/jacobi_test.c tests/test_common.c source/solvers/jacobi_point_iteration.c)
#target_include_directories(jacobi_point_test PRIVATE cplot/include)
target_link_libraries(jacobi_point_test PRIVATE jmtx m
        #cplot
        )
add_test(NAME jacobi COMMAND jacobi_point_test)

add_executable(test_fortan tests/fortran_tests/ftest_multiply.c tests/test_common.c tests/test_common.h source/f_exper/f_interface.h)
target_link_libraries(test_fortan PRIVATE jmtx m)
add_test(NAME test_fortran_interface COMMAND test_fortan)

add_executable(fjacobi_point_test tests/fortran_tests/fjacobi_test.c tests/test_common.c source/solvers/jacobi_point_iteration.c)
#target_include_directories(jacobi_point_test PRIVATE cplot/include)
target_link_libraries(fjacobi_point_test PRIVATE jmtx m
        #cplot
)
add_test(NAME fortran_jacobi COMMAND fjacobi_point_test)


if (CMAKE_BUILD_TYPE STREQUAL Debug)
    target_compile_options(jmtx PRIVATE -Wall -Wpedantic -Werror -Wextra)
elseif(CMAKE_BUILD_TYPE STREQUAL Release)
#    target_compile_definitions(jmtx PRIVATE JMTX_NO_VERIFY_PARAMS=1)
endif ()

add_executable(omp_jacobi_test tests/solver_tests/omp_jacobi_test.c tests/test_common.c)
target_link_libraries(omp_jacobi_test PRIVATE jmtx gomp)
target_compile_definitions(omp_jacobi_test PRIVATE JMTX_OPENMP)
target_compile_options(omp_jacobi_test PRIVATE -fopenmp)
add_test(NAME omp_jacobi COMMAND omp_jacobi_test)

add_executable(omp_gauss_test tests/solver_tests/omp_gauss_test.c tests/test_common.c)
target_link_libraries(omp_gauss_test PRIVATE jmtx gomp)
target_compile_definitions(omp_gauss_test PRIVATE JMTX_OPENMP)
target_compile_options(omp_gauss_test PRIVATE -fopenmp)
add_test(NAME omp_gauss COMMAND omp_gauss_test)

add_executable(omp_build_test tests/matrix_ops_tests/omp_build_test.c tests/test_common.c)
target_link_libraries(omp_build_test PRIVATE jmtx gomp)
target_compile_definitions(omp_build_test PRIVATE JMTX_OPENMP)
target_compile_options(omp_build_test PRIVATE -fopenmp)
add_test(NAME omp_build COMMAND omp_build_test)

add_executable(omp_cg_test tests/solver_tests/omp_cg_test.c tests/test_common.c)
target_link_libraries(omp_cg_test PRIVATE jmtx gomp)
target_compile_definitions(omp_cg_test PRIVATE JMTX_OPENMP)
target_compile_options(omp_cg_test PRIVATE -fopenmp)
add_test(NAME omp_cg COMMAND omp_cg_test)

