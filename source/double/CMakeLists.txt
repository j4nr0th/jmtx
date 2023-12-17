
list(APPEND MATRIX_TYPE_SOURCE_FILES
        source/double/matrices/sparse_row_compressed.c
        source/double/matrices/sparse_column_compressed.c
        source/double/matrices/sparse_conversion.c
        source/double/matrices/sparse_multiplication.c
        source/double/matrices/band_row_major.c
        source/double/matrices/sparse_diagonal_compressed.c
)

list(APPEND MATRIX_TYPE_HEADER_FILES
        include/jmtx/double/matrices/sparse_row_compressed.h
        include/jmtx/double/matrices/sparse_column_compressed.h
        include/jmtx/double/matrices/sparse_conversion.h
        include/jmtx/double/matrices/sparse_conversion_safe.h
        include/jmtx/double/matrices/sparse_multiplication.h
        include/jmtx/double/matrices/sparse_multiplication_safe.h
        include/jmtx/double/matrices/band_row_major.h
        include/jmtx/double/matrices/band_row_major_safe.h
        include/jmtx/double/matrices/sparse_diagonal_compressed.h
        include/jmtx/double/matrices/sparse_diagonal_compressed_safe.h
)

list(APPEND INTERNAL_HEADER_FILES
        source/double/matrices/sparse_row_compressed_internal.h
        source/double/matrices/sparse_column_compressed_internal.h
        source/double/matrices/band_row_major_internal.h
        source/double/matrices/sparse_diagonal_compressed_internal.h
)

list(APPEND INTERNAL_SOURCE_FILES
        source/double/matrices/basic_io.c
        source/double/matrices/basic_io.h)

list(APPEND SOLVER_SOURCE_FILES
        source/double/solvers/jacobi_point_iteration.c
        source/double/solvers/gauss_seidel_iteration.c
        source/double/solvers/bicgstab_iteration.c
        source/double/solvers/conjugate_gradient_iteration.c
        source/double/solvers/incomplete_lu_decomposition.c
        source/double/solvers/lu_solving.c
        source/double/solvers/cholesky_solving.c
        source/double/solvers/incomplete_cholesky_decomposition.c
        source/double/solvers/band_lu_decomposition.c
)

list(APPEND SOLVER_HEADER_FILES
        include/jmtx/double/solvers/jacobi_point_iteration.h
        include/jmtx/double/solvers/gauss_seidel_iteration.h
        include/jmtx/double/solvers/bicgstab_iteration.h
        include/jmtx/double/solvers/conjugate_gradient_iteration.h
        include/jmtx/double/solvers/incomplete_lu_decomposition.h
        include/jmtx/double/solvers/lu_solving.h
        include/jmtx/double/solvers/cholesky_solving.h
        include/jmtx/double/solvers/incomplete_cholesky_decomposition.h
        include/jmtx/double/solvers/band_lu_decomposition.h
)


add_executable(base_crs_test_double tests/double/matrix_ops_tests/base_crs_test.c tests/double/test_common.c)
target_link_libraries(base_crs_test_double jmtx m)
add_test(NAME base_crs_double COMMAND base_crs_test_double)

add_executable(base_ccs_test_double tests/double/matrix_ops_tests/base_ccs_test.c tests/double/test_common.c)
target_link_libraries(base_ccs_test_double jmtx m)
add_test(NAME base_ccs_double COMMAND base_ccs_test_double)

add_executable(old_crs_demo_double tests/double/really_old_tests/run_crs_tests.c tests/double/test_common.c)
target_link_libraries(old_crs_demo_double PRIVATE jmtx m)
add_test(NAME old_crs_double COMMAND old_crs_demo_double)

add_executable(old_ccs_demo_double tests/double/really_old_tests/run_ccs_tests.c tests/double/test_common.c)
target_link_libraries(old_ccs_demo_double PRIVATE jmtx m)
add_test(NAME old_ccs_double COMMAND old_ccs_demo_double)

add_executable(add_crs_test_double tests/double/matrix_ops_tests/crs_addition_test.c tests/double/test_common.c)
target_link_libraries(add_crs_test_double PRIVATE jmtx m)
add_test(NAME addition_crs_double COMMAND add_crs_test_double)

add_executable(jacobi_point_test_double tests/double/solver_tests/jacobi_test.c tests/double/test_common.c source/double/solvers/jacobi_point_iteration.c)
target_link_libraries(jacobi_point_test_double PRIVATE jmtx m)
add_test(NAME jacobi_double COMMAND jacobi_point_test_double)

add_executable(omp_jacobi_test_double tests/double/solver_tests/omp_jacobi_test.c tests/double/test_common.c)
target_link_libraries(omp_jacobi_test_double PRIVATE jmtx)
add_test(NAME omp_jacobi_double COMMAND omp_jacobi_test_double)

add_executable(gauss_test_double tests/double/solver_tests/gauss_test.c tests/double/test_common.c)
target_link_libraries(gauss_test_double PRIVATE jmtx)
add_test(NAME gauss_double COMMAND gauss_test_double)

add_executable(omp_build_test_double tests/double/matrix_ops_tests/omp_build_test.c tests/double/test_common.c)
target_link_libraries(omp_build_test_double PRIVATE jmtx)
add_test(NAME omp_build_double COMMAND omp_build_test_double)

add_executable(omp_cg_test_double tests/double/solver_tests/omp_cg_test.c tests/double/test_common.c)
target_link_libraries(omp_cg_test_double PRIVATE jmtx)
add_test(NAME omp_cg_double COMMAND omp_cg_test_double)

add_executable(incomplete_lu_test_double tests/double/solver_tests/inclu_test.c tests/double/test_common.c)
target_link_libraries(incomplete_lu_test_double PRIVATE jmtx)
add_test(NAME solve_ilu_double COMMAND incomplete_lu_test_double)

add_executable(ilu_test_double tests/double/solver_tests/ilu_test.c tests/double/test_common.c)
target_link_libraries(ilu_test_double PRIVATE jmtx)
target_compile_options(ilu_test_double PRIVATE -fopenmp)

add_executable(icho_test_double tests/double/solver_tests/incho_test.c tests/double/test_common.c)
target_link_libraries(icho_test_double PRIVATE jmtx)
add_test(NAME solve_icho_double COMMAND icho_test_double)

add_executable(icpcg_test_double tests/double/solver_tests/icpcg_test.c tests/double/test_common.c)
target_link_libraries(icpcg_test_double PRIVATE jmtx)
add_test(NAME test_icpcg_double COMMAND icpcg_test_double)

add_executable(lu_solve_test_double tests/double/matrix_ops_tests/lu_solve_test.c tests/double/test_common.c)
target_link_libraries(lu_solve_test_double PRIVATE jmtx)
add_test(NAME lu_solve_test_double COMMAND lu_solve_test_double)

add_executable(base_io_test_double tests/double/matrix_ops_tests/io_test.c tests/double/test_common.c)
target_link_libraries(base_io_test_double PRIVATE jmtx)
add_test(NAME base_io_test_double COMMAND base_io_test_double)

add_executable(base_brm_test_double tests/double/matrix_ops_tests/base_brm_test.c tests/double/test_common.c)
target_link_libraries(base_brm_test_double PRIVATE jmtx)
add_test(NAME base_brm_test_double COMMAND base_brm_test_double)

add_executable(lu_solve_test_brm_double tests/double/matrix_ops_tests/lu_solve_test_brm.c tests/double/test_common.c)
target_link_libraries(lu_solve_test_brm_double PRIVATE jmtx)
add_test(NAME lu_solve_test_brm_double COMMAND lu_solve_test_brm_double)

add_executable(brm_decompose_test_double tests/double/solver_tests/brm_decompose_test.c tests/double/test_common.c)
target_link_libraries(brm_decompose_test_double PRIVATE jmtx)
add_test(NAME brm_decompose_test_double COMMAND brm_decompose_test_double)

add_executable(brm_lu_improvement_test_double tests/double/solver_tests/brm_lu_improvement_test.c tests/double/test_common.c)
target_link_libraries(brm_lu_improvement_test_double PRIVATE jmtx)
add_test(NAME brm_lu_improvement_test_double COMMAND brm_lu_improvement_test_double)

add_executable(base_cds_test_double tests/double/matrix_ops_tests/base_cds_test.c tests/double/test_common.c)
target_link_libraries(base_cds_test_double PRIVATE jmtx)
add_test(NAME base_cds_test_double COMMAND base_cds_test_double)

add_executable(cg_cds_test_double tests/double/solver_tests/cg_cds_test.c tests/double/test_common.c)
target_link_libraries(cg_cds_test_double PRIVATE jmtx)
add_test(NAME cg_cds_double COMMAND cg_cds_test_double)

add_executable(icho_cds_test_double tests/double/solver_tests/incho_cds_test.c tests/double/test_common.c)
target_link_libraries(icho_cds_test_double PRIVATE jmtx)
add_test(NAME icho_cds_test_double COMMAND icho_cds_test_double)

add_executable(bicgstab_cds_test_double tests/double/solver_tests/bicgstab_cds_test.c tests/double/test_common.c)
target_link_libraries(bicgstab_cds_test_double PRIVATE jmtx)
add_test(NAME bicgstab_cds_test_double COMMAND bicgstab_cds_test_double)

