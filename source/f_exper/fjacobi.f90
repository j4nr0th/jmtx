! Created by jan on 29.9.2023.

module fjacobi
    use f_crs_base
    use, intrinsic :: ISO_C_Binding
    implicit none

    private
    public fortran_point_jacobi

contains

    subroutine fortran_point_jacobi(n_elements, dims, end_of_row_offsets, indices, values, tolerance, max_iterations, &
            x, y) bind(C, name="f_pt_jacobi")
        integer(C_INT32_t), intent(in), value :: n_elements, dims
        integer(C_INT32_t), intent(in) :: end_of_row_offsets(dims), indices(n_elements)
        real(C_FLOAT), intent(in), value :: tolerance
        real(C_FLOAT), target, intent(inout) :: x(dims)
        real(C_FLOAT), intent(in) :: values(n_elements)
        integer(C_INT32_T), intent(in), value :: max_iterations
        real(C_FLOAT), intent(in) :: y(dims)

        real(C_FLOAT) :: div_factors(dims)
        real(C_FLOAT) :: x_out(dims), x_swap(dims)

        integer(C_INT32_T) :: n, iter_count, i_end(dims), i_begin(dims)
        real(C_FLOAT) :: residual

        logical :: not_empty(dims)
!        allocate(div_factors(dims))
        do n = 1,dims
            div_factors(n) = 1.0 / crs_get_element(n_elements, dims, end_of_row_offsets, indices, values, n, n)
            x(n) = div_factors(n) * y(n)
        end do

        residual = 1.0
        iter_count = 1

        i_end = end_of_row_offsets(2:)
        i_begin(2:) = end_of_row_offsets(:dims - 1) + 1
        i_begin(1) = 1

        not_empty = (i_begin - 1 /= i_end)


        do while (residual >= tolerance .and. iter_count <= max_iterations)
            do n = 1,dims
                if (not_empty(n)) then
                    x_out(n) = (y(n) - crs_multiply_with_vector_row_exclude(i_end(n) - i_begin(n) + 1, dims,&
                            indices(i_begin(n):i_end(n)), values(i_begin(n):i_end(n)), x, n)) * div_factors(n)
                else
                    x_out(n) = 0
                end if
            end do

            do n = 1,dims
                x_swap(n) = crs_multiply_with_vector_row(i_end(n) - i_begin(n) + 1, dims,&
                        indices(i_begin(n):i_end(n)), values(i_begin(n):i_end(n)), x_out) ** 2
            end do

            residual = sqrt(sum(x_swap))

            x_swap = x_out
            x = x_out
            x_out = x_swap
            iter_count = iter_count + 1
        end do
    end subroutine fortran_point_jacobi

end module fjacobi