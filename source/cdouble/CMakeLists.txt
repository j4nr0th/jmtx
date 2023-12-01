
list(APPEND MATRIX_TYPE_SOURCE_FILES
        source/cdouble/matrices/sparse_row_compressed.c
        source/cdouble/matrices/sparse_column_compressed.c
        source/cdouble/matrices/sparse_conversion.c
        source/cdouble/matrices/sparse_multiplication.c
        source/cdouble/matrices/band_row_major.c
        source/cdouble/matrices/sparse_diagonal_compressed.c
)

list(APPEND MATRIX_TYPE_HEADER_FILES
        include/jmtx/cdouble/matrices/sparse_row_compressed.h
        include/jmtx/cdouble/matrices/sparse_column_compressed.h
        include/jmtx/cdouble/matrices/sparse_conversion.h
        include/jmtx/cdouble/matrices/sparse_multiplication.h
        include/jmtx/cdouble/matrices/band_row_major.h
        include/jmtx/cdouble/matrices/band_row_major_safe.h
        include/jmtx/cdouble/matrices/sparse_diagonal_compressed.h
        include/jmtx/cdouble/matrices/sparse_diagonal_compressed_safe.h
)

list(APPEND INTERNAL_HEADER_FILES
        source/cdouble/matrices/sparse_row_compressed_internal.h
        source/cdouble/matrices/sparse_column_compressed_internal.h
        source/cdouble/matrices/band_row_major_internal.h
        source/cdouble/matrices/sparse_diagonal_compressed_internal.h
)

list(APPEND INTERNAL_SOURCE_FILES
        source/cdouble/matrices/basic_io.c
        source/cdouble/matrices/basic_io.h)

list(APPEND SOLVER_SOURCE_FILES
        source/cdouble/solvers/jacobi_point_iteration.c
        source/cdouble/solvers/gauss_seidel_iteration.c
        source/cdouble/solvers/incomplete_lu_decomposition.c
        source/cdouble/solvers/lu_solving.c
        source/cdouble/solvers/cholesky_solving.c
        source/cdouble/solvers/incomplete_cholesky_decomposition.c
        source/cdouble/solvers/band_lu_decomposition.c
)

list(APPEND SOLVER_HEADER_FILES
        include/jmtx/cdouble/solvers/jacobi_point_iteration.h
        include/jmtx/cdouble/solvers/gauss_seidel_iteration.h
        include/jmtx/cdouble/solvers/incomplete_lu_decomposition.h
        include/jmtx/cdouble/solvers/lu_solving.h
        include/jmtx/cdouble/solvers/cholesky_solving.h
        include/jmtx/cdouble/solvers/incomplete_cholesky_decomposition.h
        include/jmtx/cdouble/solvers/band_lu_decomposition.h
)


add_executable(base_crs_test_cdouble tests/cdouble/matrix_ops_tests/base_crs_test.c tests/cdouble/test_common.c)
target_link_libraries(base_crs_test_cdouble jmtx m)
add_test(NAME base_crs_cdouble COMMAND base_crs_test_cdouble)

add_executable(base_ccs_test_cdouble tests/cdouble/matrix_ops_tests/base_ccs_test.c tests/cdouble/test_common.c)
target_link_libraries(base_ccs_test_cdouble jmtx m)
add_test(NAME base_ccs_cdouble COMMAND base_ccs_test_cdouble)

add_executable(old_crs_demo_cdouble tests/cdouble/really_old_tests/run_crs_tests.c tests/cdouble/test_common.c)
target_link_libraries(old_crs_demo_cdouble PRIVATE jmtx m)
add_test(NAME old_crs_cdouble COMMAND old_crs_demo_cdouble)

add_executable(old_ccs_demo_cdouble tests/cdouble/really_old_tests/run_ccs_tests.c tests/cdouble/test_common.c)
target_link_libraries(old_ccs_demo_cdouble PRIVATE jmtx m)
add_test(NAME old_ccs_cdouble COMMAND old_ccs_demo_cdouble)

add_executable(add_crs_test_cdouble tests/cdouble/matrix_ops_tests/crs_addition_test.c tests/cdouble/test_common.c)
target_link_libraries(add_crs_test_cdouble PRIVATE jmtx m)
add_test(NAME addition_crs_cdouble COMMAND add_crs_test_cdouble)

add_executable(jacobi_point_test_cdouble tests/cdouble/solver_tests/jacobi_test.c tests/cdouble/test_common.c source/cdouble/solvers/jacobi_point_iteration.c)
target_link_libraries(jacobi_point_test_cdouble PRIVATE jmtx m)
add_test(NAME jacobi_cdouble COMMAND jacobi_point_test_cdouble)

add_executable(omp_jacobi_test_cdouble tests/cdouble/solver_tests/omp_jacobi_test.c tests/cdouble/test_common.c)
target_link_libraries(omp_jacobi_test_cdouble PRIVATE jmtx)
add_test(NAME omp_jacobi_cdouble COMMAND omp_jacobi_test_cdouble)

add_executable(omp_gauss_test_cdouble tests/cdouble/solver_tests/omp_gauss_test.c tests/cdouble/test_common.c)
target_link_libraries(omp_gauss_test_cdouble PRIVATE jmtx)
add_test(NAME omp_gauss_cdouble COMMAND omp_gauss_test_cdouble)

add_executable(omp_build_test_cdouble tests/cdouble/matrix_ops_tests/omp_build_test.c tests/cdouble/test_common.c)
target_link_libraries(omp_build_test_cdouble PRIVATE jmtx)
add_test(NAME omp_build_cdouble COMMAND omp_build_test_cdouble)

add_executable(incomplete_lu_test_cdouble tests/cdouble/solver_tests/inclu_test.c tests/cdouble/test_common.c)
target_link_libraries(incomplete_lu_test_cdouble PRIVATE jmtx)
add_test(NAME solve_ilu_cdouble COMMAND incomplete_lu_test_cdouble)

add_executable(omp_ilu_test_cdouble tests/cdouble/solver_tests/omp_inclu_test.c tests/cdouble/test_common.c)
target_link_libraries(omp_ilu_test_cdouble PRIVATE jmtx)
target_compile_options(omp_ilu_test_cdouble PRIVATE -fopenmp)

add_executable(icho_test_cdouble tests/cdouble/solver_tests/incho_test.c tests/cdouble/test_common.c)
target_link_libraries(icho_test_cdouble PRIVATE jmtx)
add_test(NAME solve_icho_cdouble COMMAND icho_test_cdouble)

add_executable(lu_solve_test_cdouble tests/cdouble/matrix_ops_tests/lu_solve_test.c tests/cdouble/test_common.c)
target_link_libraries(lu_solve_test_cdouble PRIVATE jmtx)
add_test(NAME lu_solve_test_cdouble COMMAND lu_solve_test_cdouble)

add_executable(base_io_test_cdouble tests/cdouble/matrix_ops_tests/io_test.c tests/cdouble/test_common.c)
target_link_libraries(base_io_test_cdouble PRIVATE jmtx)
add_test(NAME base_io_test_cdouble COMMAND base_io_test_cdouble)

add_executable(base_brm_test_cdouble tests/cdouble/matrix_ops_tests/base_brm_test.c tests/cdouble/test_common.c)
target_link_libraries(base_brm_test_cdouble PRIVATE jmtx)
add_test(NAME base_brm_test_cdouble COMMAND base_brm_test_cdouble)

add_executable(lu_solve_test_brm_cdouble tests/cdouble/matrix_ops_tests/lu_solve_test_brm.c tests/cdouble/test_common.c)
target_link_libraries(lu_solve_test_brm_cdouble PRIVATE jmtx)
add_test(NAME lu_solve_test_brm_cdouble COMMAND lu_solve_test_brm_cdouble)

add_executable(brm_decompose_test_cdouble tests/cdouble/solver_tests/brm_decompose_test.c tests/cdouble/test_common.c)
target_link_libraries(brm_decompose_test_cdouble PRIVATE jmtx)
add_test(NAME brm_decompose_test_cdouble COMMAND brm_decompose_test_cdouble)

add_executable(brm_lu_improvement_test_cdouble tests/cdouble/solver_tests/brm_lu_improvement_test.c tests/cdouble/test_common.c)
target_link_libraries(brm_lu_improvement_test_cdouble PRIVATE jmtx)
add_test(NAME brm_lu_improvement_test_cdouble COMMAND brm_lu_improvement_test_cdouble)

add_executable(base_cds_test_cdouble tests/cdouble/matrix_ops_tests/base_cds_test.c tests/cdouble/test_common.c)
target_link_libraries(base_cds_test_cdouble PRIVATE jmtx)
add_test(NAME base_cds_test_cdouble COMMAND base_cds_test_cdouble)

