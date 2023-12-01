
list(APPEND MATRIX_TYPE_SOURCE_FILES
        source/cfloat/matrices/sparse_row_compressed.c
        source/cfloat/matrices/sparse_column_compressed.c
        source/cfloat/matrices/sparse_conversion.c
        source/cfloat/matrices/sparse_multiplication.c
        source/cfloat/matrices/band_row_major.c
        source/cfloat/matrices/sparse_diagonal_compressed.c
)

list(APPEND MATRIX_TYPE_HEADER_FILES
        include/jmtx/cfloat/matrices/sparse_row_compressed.h
        include/jmtx/cfloat/matrices/sparse_column_compressed.h
        include/jmtx/cfloat/matrices/sparse_conversion.h
        include/jmtx/cfloat/matrices/sparse_multiplication.h
        include/jmtx/cfloat/matrices/band_row_major.h
        include/jmtx/cfloat/matrices/band_row_major_safe.h
        include/jmtx/cfloat/matrices/sparse_diagonal_compressed.h
        include/jmtx/cfloat/matrices/sparse_diagonal_compressed_safe.h
)

list(APPEND INTERNAL_HEADER_FILES
        source/cfloat/matrices/sparse_row_compressed_internal.h
        source/cfloat/matrices/sparse_column_compressed_internal.h
        source/cfloat/matrices/band_row_major_internal.h
        source/cfloat/matrices/sparse_diagonal_compressed_internal.h
)

list(APPEND INTERNAL_SOURCE_FILES
        source/cfloat/matrices/basic_io.c
        source/cfloat/matrices/basic_io.h)

list(APPEND SOLVER_SOURCE_FILES
        source/cfloat/solvers/jacobi_point_iteration.c
        source/cfloat/solvers/gauss_seidel_iteration.c
        source/cfloat/solvers/incomplete_lu_decomposition.c
        source/cfloat/solvers/lu_solving.c
        source/cfloat/solvers/cholesky_solving.c
        source/cfloat/solvers/incomplete_cholesky_decomposition.c
        source/cfloat/solvers/band_lu_decomposition.c
)

list(APPEND SOLVER_HEADER_FILES
        include/jmtx/cfloat/solvers/jacobi_point_iteration.h
        include/jmtx/cfloat/solvers/gauss_seidel_iteration.h
        include/jmtx/cfloat/solvers/incomplete_lu_decomposition.h
        include/jmtx/cfloat/solvers/lu_solving.h
        include/jmtx/cfloat/solvers/cholesky_solving.h
        include/jmtx/cfloat/solvers/incomplete_cholesky_decomposition.h
        include/jmtx/cfloat/solvers/band_lu_decomposition.h
)


add_executable(base_crs_test_cfloat tests/cfloat/matrix_ops_tests/base_crs_test.c tests/cfloat/test_common.c)
target_link_libraries(base_crs_test_cfloat jmtx m)
add_test(NAME base_crs_cfloat COMMAND base_crs_test_cfloat)

add_executable(base_ccs_test_cfloat tests/cfloat/matrix_ops_tests/base_ccs_test.c tests/cfloat/test_common.c)
target_link_libraries(base_ccs_test_cfloat jmtx m)
add_test(NAME base_ccs_cfloat COMMAND base_ccs_test_cfloat)

add_executable(old_crs_demo_cfloat tests/cfloat/really_old_tests/run_crs_tests.c tests/cfloat/test_common.c)
target_link_libraries(old_crs_demo_cfloat PRIVATE jmtx m)
add_test(NAME old_crs_cfloat COMMAND old_crs_demo_cfloat)

add_executable(old_ccs_demo_cfloat tests/cfloat/really_old_tests/run_ccs_tests.c tests/cfloat/test_common.c)
target_link_libraries(old_ccs_demo_cfloat PRIVATE jmtx m)
add_test(NAME old_ccs_cfloat COMMAND old_ccs_demo_cfloat)

add_executable(add_crs_test_cfloat tests/cfloat/matrix_ops_tests/crs_addition_test.c tests/cfloat/test_common.c)
target_link_libraries(add_crs_test_cfloat PRIVATE jmtx m)
add_test(NAME addition_crs_cfloat COMMAND add_crs_test_cfloat)

add_executable(jacobi_point_test_cfloat tests/cfloat/solver_tests/jacobi_test.c tests/cfloat/test_common.c source/cfloat/solvers/jacobi_point_iteration.c)
target_link_libraries(jacobi_point_test_cfloat PRIVATE jmtx m)
add_test(NAME jacobi_cfloat COMMAND jacobi_point_test_cfloat)

add_executable(omp_jacobi_test_cfloat tests/cfloat/solver_tests/omp_jacobi_test.c tests/cfloat/test_common.c)
target_link_libraries(omp_jacobi_test_cfloat PRIVATE jmtx)
add_test(NAME omp_jacobi_cfloat COMMAND omp_jacobi_test_cfloat)

add_executable(omp_gauss_test_cfloat tests/cfloat/solver_tests/omp_gauss_test.c tests/cfloat/test_common.c)
target_link_libraries(omp_gauss_test_cfloat PRIVATE jmtx)
add_test(NAME omp_gauss_cfloat COMMAND omp_gauss_test_cfloat)

add_executable(omp_build_test_cfloat tests/cfloat/matrix_ops_tests/omp_build_test.c tests/cfloat/test_common.c)
target_link_libraries(omp_build_test_cfloat PRIVATE jmtx)
add_test(NAME omp_build_cfloat COMMAND omp_build_test_cfloat)

add_executable(incomplete_lu_test_cfloat tests/cfloat/solver_tests/inclu_test.c tests/cfloat/test_common.c)
target_link_libraries(incomplete_lu_test_cfloat PRIVATE jmtx)
add_test(NAME solve_ilu_cfloat COMMAND incomplete_lu_test_cfloat)

add_executable(omp_ilu_test_cfloat tests/cfloat/solver_tests/omp_inclu_test.c tests/cfloat/test_common.c)
target_link_libraries(omp_ilu_test_cfloat PRIVATE jmtx)
target_compile_options(omp_ilu_test_cfloat PRIVATE -fopenmp)

add_executable(icho_test_cfloat tests/cfloat/solver_tests/incho_test.c tests/cfloat/test_common.c)
target_link_libraries(icho_test_cfloat PRIVATE jmtx)
add_test(NAME solve_icho_cfloat COMMAND icho_test_cfloat)

add_executable(lu_solve_test_cfloat tests/cfloat/matrix_ops_tests/lu_solve_test.c tests/cfloat/test_common.c)
target_link_libraries(lu_solve_test_cfloat PRIVATE jmtx)
add_test(NAME lu_solve_test_cfloat COMMAND lu_solve_test_cfloat)

add_executable(base_io_test_cfloat tests/cfloat/matrix_ops_tests/io_test.c tests/cfloat/test_common.c)
target_link_libraries(base_io_test_cfloat PRIVATE jmtx)
add_test(NAME base_io_test_cfloat COMMAND base_io_test_cfloat)

add_executable(base_brm_test_cfloat tests/cfloat/matrix_ops_tests/base_brm_test.c tests/cfloat/test_common.c)
target_link_libraries(base_brm_test_cfloat PRIVATE jmtx)
add_test(NAME base_brm_test_cfloat COMMAND base_brm_test_cfloat)

add_executable(lu_solve_test_brm_cfloat tests/cfloat/matrix_ops_tests/lu_solve_test_brm.c tests/cfloat/test_common.c)
target_link_libraries(lu_solve_test_brm_cfloat PRIVATE jmtx)
add_test(NAME lu_solve_test_brm_cfloat COMMAND lu_solve_test_brm_cfloat)

add_executable(brm_decompose_test_cfloat tests/cfloat/solver_tests/brm_decompose_test.c tests/cfloat/test_common.c)
target_link_libraries(brm_decompose_test_cfloat PRIVATE jmtx)
add_test(NAME brm_decompose_test_cfloat COMMAND brm_decompose_test_cfloat)

add_executable(brm_lu_improvement_test_cfloat tests/cfloat/solver_tests/brm_lu_improvement_test.c tests/cfloat/test_common.c)
target_link_libraries(brm_lu_improvement_test_cfloat PRIVATE jmtx)
add_test(NAME brm_lu_improvement_test_cfloat COMMAND brm_lu_improvement_test_cfloat)

add_executable(base_cds_test_cfloat tests/cfloat/matrix_ops_tests/base_cds_test.c tests/cfloat/test_common.c)
target_link_libraries(base_cds_test_cfloat PRIVATE jmtx)
add_test(NAME base_cds_test_cfloat COMMAND base_cds_test_cfloat)

