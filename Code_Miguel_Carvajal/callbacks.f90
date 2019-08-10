module callbacks
    use types
    implicit none
     ! declare an abstract type called Functor
    
    type, abstract :: Functor
        contains
        procedure(func), deferred :: eval
    end type

    type, abstract :: VectorFunctor
        contains
        procedure(vecfunc), deferred :: eval
    end type

    abstract interface
        function func(this, x) result(fx)
            import :: Functor, dp
            class(Functor) :: this
            real(dp), intent(in) :: x(:)
            real(dp) :: fx
        end function func
        !receive a vector and resturns a vector
        function vecfunc(this, x) result(fx)
            import :: VectorFunctor, dp
            class(VectorFunctor) :: this
            real(dp), intent(in) :: x(:)
            real(dp),dimension(:), allocatable :: fx
        end function vecfunc

    end interface
end module callbacks