cmake_minimum_required(VERSION 3.10)

project(jmtx)
set(CMAKE_C_STANDARD 23)

list(APPEND MATRIX_TYPE_SOURCE_FILES
        source/matrices/sparse_row_compressed.c
        source/matrices/sparse_column_compressed.c
        source/matrices/dense_row_major.c
        source/matrices/dense_col_major.c
        source/matrices/matrix_base.c
        source/matrices/sparse_conversion.c
        source/matrices/sparse_multiplication.c
        source/matrices/band_row_major.c
        source/matrices/sparse_diagonal_compressed.c
        )

list(APPEND MATRIX_TYPE_HEADER_FILES
        include/jmtx/matrices/sparse_row_compressed.h
        include/jmtx/matrices/sparse_column_compressed.h
        source/matrices/dense_row_major.h
        source/matrices/dense_col_major.h
        include/jmtx/matrices/matrix_base.h
        include/jmtx/matrices/sparse_conversion.h
        include/jmtx/matrices/sparse_multiplication.h
        include/jmtx/matrices/band_row_major.h
        include/jmtx/matrices/band_row_major_safe.h
        source/matrices/sparse_diagonal_compressed.h
        source/matrices/sparse_diagonal_compressed_safe.h
        )

list(APPEND INTERNAL_HEADER_FILES
        source/matrices/matrix_base_internal.h
        source/matrices/sparse_row_compressed_internal.h
        source/matrices/sparse_column_compressed_internal.h
        source/matrices/band_row_major_internal.h
        source/matrices/sparse_diagonal_compressed_internal.h
        )

list(APPEND INTERNAL_SOURCE_FILES )

list(APPEND SOLVER_SOURCE_FILES
        source/solvers/jacobi_point_iteration.c
        source/solvers/gauss_seidel_iteration.c
        source/solvers/bicgstab_iteration.c
        source/solvers/conjugate_gradient_iteration.c
        source/solvers/incomplete_lu_decomposition.c
        source/solvers/lu_solving.c
        source/solvers/cholesky_solving.c
        source/solvers/incomplete_cholesky_decomposition.c
        source/solvers/band_lu_decomposition.c
        )

list(APPEND SOLVER_HEADER_FILES
        include/jmtx/solvers/jacobi_point_iteration.h
        include/jmtx/solvers/gauss_seidel_iteration.h
        source/solvers/bicgstab_iteration.h
        include/jmtx/solvers/conjugate_gradient_iteration.h
        include/jmtx/solvers/incomplete_lu_decomposition.h
        include/jmtx/solvers/lu_solving.h
        include/jmtx/solvers/cholesky_solving.h
        include/jmtx/solvers/incomplete_cholesky_decomposition.h
        include/jmtx/solvers/solver_base.h
        include/jmtx/solvers/band_lu_decomposition.h
        )


add_library(jmtx
        ${MATRIX_TYPE_SOURCE_FILES} ${MATRIX_TYPE_HEADER_FILES}
        ${SOLVER_SOURCE_FILES} ${SOLVER_HEADER_FILES}
        ${INTERNAL_SOURCE_FILES} ${INTERNAL_HEADER_FILES}
        ${FORTRAN_SOURCE_FILES}
        include/jmtx/common.h
        source/matrices/basic_io.c
        source/matrices/basic_io.h
        )

target_link_libraries(jmtx gomp m)

target_compile_options(jmtx PRIVATE -march=native)

option(JMTX_USE_OPENMP "Use OpenMP for parallel functions" ON)

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_compile_options(jmtx PRIVATE -fopenmp -fanalyzer)
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

add_executable(old_crs_demo tests/really_old_tests/run_crs_tests.c tests/test_common.c)
target_link_libraries(old_crs_demo PRIVATE jmtx m)
add_test(NAME old_crs COMMAND old_crs_demo)

add_executable(old_ccs_demo tests/really_old_tests/run_ccs_tests.c tests/test_common.c)
target_link_libraries(old_ccs_demo PRIVATE jmtx m)
add_test(NAME old_ccs COMMAND old_ccs_demo)

add_executable(old_compare_demo tests/really_old_tests/compare_crs_and_ccs.c tests/test_common.c)
target_link_libraries(old_compare_demo PRIVATE jmtx m)
add_test(NAME old_compare COMMAND old_compare_demo)

add_executable(add_crs_test tests/matrix_ops_tests/crs_addition_test.c tests/test_common.c)
target_link_libraries(add_crs_test PRIVATE jmtx m)
add_test(NAME addition_crs COMMAND add_crs_test)

add_executable(jacobi_point_test tests/solver_tests/jacobi_test.c tests/test_common.c source/solvers/jacobi_point_iteration.c)
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

add_executable(incomplete_lu_test tests/solver_tests/inclu_test.c tests/test_common.c)
target_link_libraries(incomplete_lu_test PRIVATE jmtx gomp)
target_compile_definitions(incomplete_lu_test PRIVATE JMTX_OPENMP)
target_compile_options(incomplete_lu_test PRIVATE -fopenmp)
add_test(NAME solve_ilu COMMAND incomplete_lu_test)

add_executable(omp_ilu_test tests/solver_tests/omp_inclu_test.c tests/test_common.c)
target_link_libraries(omp_ilu_test PRIVATE jmtx gomp)
target_compile_definitions(omp_ilu_test PRIVATE JMTX_OPENMP)
target_compile_options(omp_ilu_test PRIVATE -fopenmp)

add_executable(icho_test tests/solver_tests/incho_test.c tests/test_common.c)
target_link_libraries(icho_test PRIVATE jmtx gomp)
target_compile_definitions(icho_test PRIVATE JMTX_OPENMP)
target_compile_options(icho_test PRIVATE -fopenmp)
add_test(NAME solve_icho COMMAND icho_test)

add_executable(icpcg_test tests/solver_tests/icpcg_test.c tests/test_common.c)
target_link_libraries(icpcg_test PRIVATE jmtx gomp)
target_compile_definitions(icpcg_test PRIVATE JMTX_OPENMP)
target_compile_options(icpcg_test PRIVATE -fopenmp)
add_test(NAME test_icpcg COMMAND icpcg_test)

add_executable(lu_solve_test tests/matrix_ops_tests/lu_solve_test.c tests/test_common.c)
target_link_libraries(lu_solve_test PRIVATE jmtx gomp)
target_compile_definitions(lu_solve_test PRIVATE JMTX_OPENMP)
target_compile_options(lu_solve_test PRIVATE -fopenmp)
add_test(NAME lu_solve_test COMMAND lu_solve_test)

add_executable(base_io_test tests/matrix_ops_tests/io_test.c tests/test_common.c)
target_link_libraries(base_io_test PRIVATE jmtx)
add_test(NAME base_io_test COMMAND base_io_test)

add_executable(base_brm_test tests/matrix_ops_tests/base_brm_test.c tests/test_common.c)
target_link_libraries(base_brm_test PRIVATE jmtx)
add_test(NAME base_brm_test COMMAND base_brm_test)

add_executable(lu_solve_test_brm tests/matrix_ops_tests/lu_solve_test_brm.c tests/test_common.c)
target_link_libraries(lu_solve_test_brm PRIVATE jmtx)
add_test(NAME lu_solve_test_brm COMMAND lu_solve_test_brm)

add_executable(brm_decompose_test tests/solver_tests/brm_decompose_test.c tests/test_common.c)
target_link_libraries(brm_decompose_test PRIVATE jmtx)
add_test(NAME brm_decompose_test COMMAND brm_decompose_test)

add_executable(brm_lu_improvement_test tests/solver_tests/brm_lu_improvement_test.c tests/test_common.c)
target_link_libraries(brm_lu_improvement_test PRIVATE jmtx)
add_test(NAME brm_lu_improvement_test COMMAND brm_lu_improvement_test)
