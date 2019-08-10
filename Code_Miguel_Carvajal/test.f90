module tests
    use spin_system
    use utils
    use types
    use callbacks
    use eigen
    implicit none

    type, extends(VectorFunctor) :: HeisenbergHamiltonianOperator
        real(dp),allocatable :: M(:,:)
    contains
        procedure :: eval => VectorMatrixMult
    end type

    

    contains

    function VectorMatrixMult(this,x) result(fx)
        class(HeisenbergHamiltonianOperator) :: this
        real(dp),intent(in) :: x(:)
        


        real(dp),dimension(:), allocatable :: fx
        integer :: i
        
        allocate(fx(size(x,1)))
        
        !standard matrix vec mutliplication
        fx = matmul(this%M,x)
    end function VectorMatrixMult

    subroutine test()

        integer, dimension(4) :: b2
        integer, dimension(8) :: b3
        integer, dimension(16) :: b4
        real(dp) :: Psi(2),NPsi(2)
        integer ::f1,f2
        logical :: valid
        real(dp) :: weight
        b2 = (/1,2,0,0/)
        b4 = (/3,5,6,9,10,12,0,0,0,0,0,0,0,0,0,0/)

        print *,"Testing..."
        call InitSpinSystem(2)
        
        call assert(ALL(Basis == b2))
        ! call OperatorSiSiNext(1,weight)
        
        Psi = 0
        Psi(1) = 1
        print *,Basis
        call OperateHamiltonian(Psi,NPsi)
        print *,NPsi

        print *, "All tests pass"



        call InitSpinSystem(4)

        call assert(ALL(Basis == b4))

        ! test spin flip
        call OperatorSpinFlip(0,1,f2,valid)
        call assert(valid)!, "not valid")
        call assert(f2 == 2)!, "not 2")

        call OperatorSpinFlip(3,8,f2,valid)
        call assert(valid)!, "not valid")
        call assert(f2 == 1)!, "not 1")
        call OperatorSpinFlip(1,1,f2,valid)
        call assert(.not. valid)!, "valid")
        

    end subroutine test

    subroutine test_eigen
        type(HeisenbergHamiltonianOperator) :: op
        real(dp) :: M1(2,2), M2(5,5), M3(100,100)
        real(dp) :: E1(2), v1(2,2), E2(5),V2(5,5), E3(100), V3(100,100)
        integer :: i
        M3 = 0
        do i = 1,100
            M3(i,i) = -20*i
            if ( i >1) M3(i,i - 1) = 0.1
            if (i < 100) M3(i,i +1) = 0.1
        enddo


        M1(1,:) = [1,1]
        M1(2,:) = [1,3]
        ! M(3,:) = [1,1]
        op%M = M1
        ! call symetric_eigenvalues(2,M,v)


        
     ! call lanczos(2,[ 1.0_dp,0.0_dp],op,E1,v1)
     ! print *,"EigenValues = ",E1
     ! call symetric_eigenvalues(2,M1,E1)
     ! print *,"EigenValues ExactDiag = ",E1

     M2 = reshape([1.96,  -6.49, -0.47,  -7.20,  -0.65,&
                  -6.49,  3.80,  -6.39,  1.50, -6.34,&
                  -0.47, -6.39,  4.17,  -1.51,  2.67,&
                   -7.20,  1.50, -1.51,  5.70,  1.80,&
                    -0.65, -6.34,  2.67,  1.80, -7.10],[5,5])

     op%M = M2
     
    ! call random_number(E2) 
    ! print *, E2
    ! ! call lanczos(5,[1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp],op,E2,v2)
    ! call lanczos(5,E2,op,E2,v2)
    ! print *,"EigenValues = ",E2
    ! call symetric_eigenvalues(5,M2,E2)
    ! print *,"EigenValues  ExactDiag = ",E2
    ! call random_number(E2) 
    ! print *, E2
    ! ! call lanczos(5,[1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp],op,E2,v2)
    ! call lanczos(5,E2,op,E2,v2)
    ! print *,"EigenValues = ",E2
    ! call symetric_eigenvalues(5,M2,E2)
    ! print *,"EigenValues  ExactDiag = ",E2

    op%M = M3
    print *, "M3"
    call random_number(E3)
    
    call lanczos(100,E3,op,E3,v3)
    print *,"EigenValues = ",E3(1)
    call symetric_eigenvalues(100,M3,E3)
    print *,"EigenValues  ExactDiag = ",minval(E3)

    print *, "result",E3(1)


    end subroutine test_eigen

end module tests