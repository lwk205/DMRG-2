module problem 
    use spin_system
    use types
    use eigen
    use constants
    use utils
    use callbacks
    ! implement the task for the assignment
    implicit none

    type, extends(VectorFunctor) :: HeisenbergHamiltonianOperator

    contains
        procedure :: eval => VectorMatrixMult
    end type


private 
    
public prob1, VectorMatrixMult
contains

    function VectorMatrixMult(this,x) result(fx)
        class(HeisenbergHamiltonianOperator) :: this
        real(dp),intent(in) :: x(:)
        real(dp),dimension(:), allocatable :: fx
        allocate(fx(BasisDimension))
        call OperateHamiltonian(x, fx)
    end function VectorMatrixMult

    subroutine prob1
        type(HeisenbergHamiltonianOperator) :: op
        real(dp)  :: Energy
        real(dp),allocatable,dimension(:) :: StateVector,NewStateVector
        integer :: NumSites
        integer :: FileH

        open(newunit=FileH,file="out.txt",status="replace")

        do NumSites = 4,16,2
            call InitSpinSystem(NumSites)  
            
            
            if(allocated(StateVector))deallocate(StateVector)
            allocate(StateVector(BasisDimension))
            call random_number(StateVector)
            print *, "====================="
            print *,"Starting lanczos"
            print *,"Numero de sitios: ",NumSites 
            print *, "BasisDimension:",BasisDimension 


            call lanczos(BasisDimension,StateVector,op,Energy,StateVector)


            !normalize 
            StateVector = StateVector /norm2(StateVector)

            call prob1c(BasisDimension,StateVector,NumSites)
        

            ! NewStateVector = StateVector
            
            print *,NumSites,Energy,Energy/NumSites
            ! print *, "State Vector", StateVector
            
            !check diagonalization
            ! call OperateHamiltonian(StateVector,NewStateVector)
            ! print *, "NewStateVector",NewStateVector/StateVector
            

            write (FileH,*) NumSites,Energy,Energy/NumSites
            
            call ComputeCorrelation(NumSites, StateVector)
        enddo


        contains 
            subroutine ComputeCorrelation(N,state)
                real(dp) :: state(:)
                integer,intent(in)  :: N !num sites
                real(dp),allocatable :: correlation(:),aux(:)
                integer :: l
                integer :: CORRELATION_UNIT
                allocate(correlation(0:N))


                open(newunit = CORRELATION_UNIT, file="correlation"//str(N)//".dat",status="replace")
                do l = 0,N-1
                    aux = state
                    call OperateSz(aux,l)
                    call OperateSz(aux,0)
                    correlation(l) = dot_product(aux,state)   

                    write (CORRELATION_UNIT,*) l,correlation(l)

                enddo
                call ComputeStructureFactor(N, correlation)
                close(CORRELATION_UNIT)

                
            end subroutine ComputeCorrelation

            subroutine ComputeStructureFactor(N, correlation)
                implicit none
                integer,intent(in) :: N
                real(dp) :: correlation(N)
                real(dp) :: Rl(N)
                real(dp) :: q
                integer :: i
                complex :: img_number = (0.0,1.0)
                complex(dp) :: S(N)
                integer :: STRUCTURE_UNIT

                open(newunit = STRUCTURE_UNIT, file= "structure"//str(N)//".dat",status="replace")

                print *,"Structure factor, N = " , N
                print *, "q", "                     S(q)"
                Rl = (/ (i, i=0,N-1)/)
                do i=0,N-1
                    q = 2 * pi * i/N;
                    S(i + 1) =  sum( exp(- img_number * q * Rl)*correlation)
                    print*, q, S(i+1)
                    write (STRUCTURE_UNIT,*) q, REALPART(S(i+1)), IMAGPART(S(i+1))
                enddo
                close(STRUCTURE_UNIT)

            end subroutine ComputeStructureFactor

            subroutine prob1c(N,GroundState,Sites)
                implicit none
                integer,intent(in) :: N,Sites
                real(dp) :: GroundState(N)
                real(dp) :: ExitedState(N), FinalState(N)
                integer :: i,j
                type(HeisenbergHamiltonianOperator) :: op
                real(dp) :: StartTime, EndTime, TimeStep, TimePoints(100)
                real(dp) :: EvolutionPoints(100)
                integer :: EVOLUTION_UNIT

                open(newunit=EVOLUTION_UNIT, file="evolution"//str(Sites)//".dat",status="replace");

                StartTime = 0;
                EndTime = 0.5
                TimeStep = 0.005
                ! apply operator $S^z_{q=\pi}$
                ExitedState = GroundState
                do i=0,Sites-1
                    call OperateSz(ExitedState,i)
                    ExitedState = ExitedState*(-1)**i
                enddo   
                ! now i have the excited state i must 
                TimePoints = [(StartTime + i*TimeStep,i=0,100-1)]
            
                call LanczosTimeEvolution(N,GroundState,ExitedState, TimePoints,op,EvolutionPoints)
                print *,"Writing evolution to a file"
                do i=1,100
                    write(EVOLUTION_UNIT,*) TimePoints(i),EvolutionPoints(i)
                enddo
                close(EVOLUTION_UNIT)

            end subroutine prob1c

    end subroutine prob1




end module problem