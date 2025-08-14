module polynomial_class

    use numPrecision
    use genericProcedures,      only : fatalError

    implicit none
    private

    type, public :: polynomial
        private
        real(defReal), allocatable :: coefficients(:)
        integer(shortInt)          :: order = -1
    contains

        procedure :: build
        generic   :: evaluate => evaluate_scalar, evaluate_array
        procedure :: deg
        procedure :: evaluate_scalar
        procedure :: evaluate_array
        procedure :: kill
    end type polynomial

contains

    ! Initializes the polynomial with given coefficients
    subroutine build(self, coeffs)
        class(polynomial), intent(inout) :: self
        real(defReal), intent(in)        :: coeffs(:)

        if (allocated(self % coefficients)) deallocate(self % coefficients)

        allocate(self % coefficients, source=coeffs)
        self % order = size(self % coefficients) - 1
    end subroutine build

    ! Evaluates the polynomial at a given x-value
    pure function evaluate_scalar(self, x) result(val)
        class(polynomial), intent(in) :: self
        real(defReal), intent(in)     :: x
        real(defReal)                 :: val
        integer(shortInt)             :: i

        val = ZERO
        do i = 1, self % order + 1
            val = val + self % coefficients(i) * x**(i - 1)
        end do
    end function evaluate_scalar

    pure function evaluate_array(self, x) result(val)
        class(polynomial), intent(in) :: self
        real(defReal), intent(in)     :: x(:)
        real(defReal)                 :: val(size(x))
        integer(shortInt)             :: i

        val = ZERO
        do i = 1, self % order + 1
            val = val + self % coefficients(i) * x**(i - 1)
        end do
    end function evaluate_array

    ! Returns the polynomial order
    pure function deg(self) result(val)
        class(polynomial), intent(in)   :: self
        integer(shortInt)               :: val

        val = self % order
    end function deg 

    ! Deallocates the polynomial and resets it
    elemental subroutine kill(self)
        class(polynomial), intent(inout) :: self

        if (allocated(self % coefficients)) deallocate(self % coefficients)
        self % order = -1
    end subroutine kill

end module polynomial_class