module eigen
    use types
    use utils !assert
    use callbacks
    implicit none
    private 
    public lanczos,symetric_eigenvalues,matvecmul

contains


    subroutine lanczos(Dim, StartVec, Evaluator, EigenValue,EigenVector)
        !direct algorithm to find the lowest eigenvalue of a system
        ! http://www.netlib.org/utk/people/JackDongarra/etemplates/node103.html

        integer,intent(in) :: Dim
        real(dp),intent(in) :: StartVec(Dim)
        class(VectorFunctor) :: Evaluator
        real(dp),intent(out) :: EigenVector(Dim)
        real(dp),intent(out) :: EigenValue
        real(dp),allocatable :: Tj(:,:), RitzValues(:),RitzVectors(:,:)
        real(dp),allocatable :: RitzNorms(:)
        
        !private variables
        real(dp) :: V(Dim,Dim) ! the matrix V that holds the orthonormal basis {v,Av,A^2v,...}
        real(dp) :: r(Dim) ! aux vector for computation
        integer :: j
        real(dp) :: Beta(1:Dim),Alpha(1:Dim) !banded elements of the matrix

        real(dp),parameter :: eps = tiny(1.0_dp), thres = sqrt(eps)
        real(dp) :: NormA !aproximate norm of matrix A


        Beta(1) =  norm2(StartVec)
        
        
        V(:,1) = StartVec /Beta(1)
        



        do j=1,Dim - 1

            r = Evaluator%eval(V(:,j))
            

            alpha(j) = DOT_PRODUCT(V(:,j),r)
            if (j > 1)r = r - V(:,j-1)*Beta(j)
            r =  r - V(:,j)*alpha(j)
            
            Beta(j + 1) = NORM2(r)
            V(:,j+1) = r /Beta(j+1)

            ! Reorthogonalize
            ! Here we do full reorthogonalization
            ! http://www.netlib.org/utk/people/JackDongarra/etemplates/node109.html
            call FullReorthogonalization(j,V)
            if(j > 2) then
                
                call ComputeRitzNorms(j)
                if( RitzNorms(1) < 1e-5) then
                     EigenValue = RitzValues(1)
                     EigenVector = matvecmul(Dim,j,V(:,1:j),RitzVectors(:,1) )
                     
                     print *,"Iters" , j
                    return
                endif
            endif
            ! call SelectiveReorthogonalization(j)
            
            
        enddo

        r = Evaluator%eval(V(:,Dim))
        alpha(Dim) =   dot_product(r,V(:,Dim))
        call ComputeRitzNorms(Dim)

        EigenValue = RitzValues(1)

    contains

            subroutine FullReorthogonalization(j, V)
                integer,intent(in) :: j
                real(dp),intent(inout) :: V(:,:)
                real(dp),dimension(j,size(V,1)) :: TransposedV
                real(dp),dimension(size(V,1)) :: x
                real(dp),dimension(j) :: h
                real(dp) :: Norm
                integer :: i,k
                ! print *, "Starting full reorthonormalization"
                x = V(:,j+1)
                do i =1,j
                    do k = 1,size(V,1)
                        TransposedV(i,k) = V(k,i)
                    enddo
                enddo
                Norm = norm2(x)
                ! compute reortogonalization coeff
                h = matmul(TransposedV,x)
                V(:,j+1) = x - matmul(V(:,1:j),h)
            end subroutine FullReorthogonalization


            subroutine SelectiveReorthogonalization(j)
                ! performs selective reortonormalization using 
                ! a recurrence relation described here
                ! http://www.netlib.org/utk/people/JackDongarra/etemplates/node110.html
                !
                integer,intent(in) :: j
                integer :: i
                integer :: k !loop variable
                real(dp),save,dimension(:), allocatable :: WCurrent, WNext, WPrev
                logical, save :: FirstCall = .true.
                ! print *, "Starting selective reorthogonalization"

                if (FirstCall) then
                    WPrev = [0.0_dp, 0.0_dp]
                    WCurrent = [0.0_dp, eps, eps]
                    allocate(WNext(4))
                    FirstCall = .false.
                endif
            
                ! print *, "Calling with j ",j
                !compute WNext using WCurrent && WPrev
                WNext(1) = 0
                do k = 2, j
                    WNext(k) = beta(k)* WCurrent(k+1) + (alpha(k)- alpha(j)) * WCurrent(k) + &
                    beta(k-1)*WCurrent(k-1) - beta(j-1) * WPrev(k)
                enddo

                call UpdateNormA(j)
                WNext(2:j) = (WNext(2:j) + sign([(1.0_dp,i=2,j)],WNext(2:j))* eps * NormA)/Beta( j + 1)
                WNext(j + 1 ) = NormA * eps
                WNext(j + 2 ) = eps 

                if (maxval(WNext) > thres) then
                    ! reorthogonalize
                    ! print *, "Reorthogonalize"
                    call FullReorthogonalization(j,V)

                endif

                WPrev = WCurrent
                WCurrent = WNext
                deallocate(WNext)
                allocate(WNext(j+3))

            end subroutine SelectiveReorthogonalization

            subroutine UpdateNormA(j)
                ! Computes an aproximation of the norm of A
                ! using the norm of T
                integer :: j
                ! if (j == 1) then
                !     NormA = abs(alpha(1) + beta(2))
                ! else
                !     NormA = max(NormA, abs(beta(j) + alpha(j) + beta(j+1)))
                ! endif
            end subroutine UpdateNormA


            subroutine ComputeRitzNorms(j)
                integer,intent(in) :: j
                integer :: i
                if (allocated(Tj)) deallocate(Tj)
                if(allocated(RitzValues)) deallocate(RitzValues)
                if(allocated(RitzVectors)) deallocate(RitzVectors)
                allocate(Tj(j,j))
                Tj = 0

                allocate(RitzValues(j))
                allocate(RitzVectors(j,j))
                ! fill T matrix
                
                do i=1,j
                    Tj(i,i) = alpha(i)
                    if (i < j ) Tj(i,i + 1) = beta(i + 1)
                enddo
                
                
                ! solve eigen problem
                call symetric_eigen(j,Tj,RitzValues,RitzVectors)

                if(allocated(RitzNorms)) deallocate(RitzNorms)
                allocate(RitzNorms(j))

                RitzNorms = abs(Beta(1:j)*RitzVectors(j,:))
                ! print *, "RitzNorms", RitzNorms

            end subroutine ComputeRitzNorms


    end subroutine lanczos

    subroutine symetric_eigenvalues(Dim,A, Values)
        integer :: Dim
        real(dp),intent(in) :: A(Dim,Dim)
        real(dp) :: A_(Dim,Dim)
        real(dp),intent(out) :: Values(Dim)
        integer,parameter :: LWORK_MAX = 1000
        real(dp) :: Work(LWORK_MAX)
        integer :: lWork, info
        A_ = A
        !query optimal lwork
        lWork = -1
        
        call dsyev("N", "L",Dim,A_,Dim,Values,Work,lWork,info)
        lWork = min(LWORK_MAX, int(Work(1)))
        
        !solve eigenproblem
        call dsyev("N", "U", Dim, A_, Dim, Values,Work,lWork,info)
        
        call assert(info == 0)!, "symetric_eigenvalues: diagonalization failed")


    end subroutine symetric_eigenvalues



    subroutine symetric_eigen(Dim,A, Values,Vectors)
        ! symetric_eigen solves the eigebproblem for
        ! a symetric real matrix A using LAPACK dsyev
        ! the eignevalues are stored in Values
        ! and the eigenvectors are stored in Vectors
        integer :: Dim
        real(dp) :: A(Dim,Dim)
        real(dp) :: A_(Dim,Dim)
        real(dp),intent(out) :: Values(Dim), Vectors(Dim,Dim)
        integer,parameter :: LWORK_MAX = 1000
        real(dp) :: Work(LWORK_MAX)
        integer :: lWork, info
        integer :: i
        
        A_ = A
        !query optimal lwork
        lWork = -1
        call dsyev("V", "L",Dim,A_,Dim,Values,Work,lWork,info)
        lWork = min(LWORK_MAX, int(Work(1)))
        
        !solve eigenproblem
        call DSYEV("V", "U", Dim, A_, Dim, Values,Work,lWork,info)
        
        
        call assert(info == 0)!, "symetric_eigenvalues: diagonalization failed")

     
        Vectors = A_

    end subroutine symetric_eigen

    function matvecmul(M,N,A,v) result(rs)
        implicit none
        integer,intent(in) :: M,N
        real(dp):: A(M,N)
        real(dp):: v(N)
        real(dp):: rs(M)
        integer :: i,j
        
        do i = 1,M
            rs(i) = 0
            do j= 1,N
                rs(i) =  A(i,j)*v(j) + rs(i)
            enddo
        enddo

    end function matvecmul

    subroutine LanczosTimeEvolution(InitialState, Time,op,FinalState)
    end subroutine LanczosTimeEvolution




end module eigen