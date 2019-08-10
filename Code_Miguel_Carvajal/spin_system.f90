module spin_system
    use types
    implicit none
    integer :: NumSites
    integer,dimension(:),allocatable :: Basis, BasisLookup
    integer :: BasisDimension
    private 
    public InitSpinSystem,Basis,BasisDimension,OperatorSpinFlip,OperateHamiltonian,OperateSz


contains

    subroutine OperateHamiltonian(PsiOld, PsiNew)
        ! performs matrix vector multiplication PsiNew = H * PsiOld
        ! the hamitonian is the Hisemberg hamitonian for spin 1/2 systems
        real(dp),intent(in)  :: PsiOld(:)
        real(dp),intent(out) :: PsiNew(:)
        integer :: BasisFunction, NewBasisFunc, NewBasisFuncPos
        real(dp) :: sum
        logical :: valid

        integer:: i,j
        PsiNew = 0
        do i = 1,BasisDimension
            BasisFunction = Basis(i) 
            call OperatorSzSzNext(BasisFunction, sum)
            PsiNew(i) =  PsiNew(i) + PsiOld(i) * sum
            do j= 0,NumSites-2
                call OperatorSpinFlip(j,BasisFunction,NewBasisFunc, valid)
                if (valid) then
                    NewBasisFuncPos = BasisLookup(NewBasisFunc)
                    PsiNew(NewBasisFuncPos) = PsiNew(NewBasisFuncPos) + 0.5 * PsiOld(i)
                endif
            enddo
        enddo
    end subroutine OperateHamiltonian

    subroutine OperatorSpinFlip(Site,BasisFunc,NewBasisFunc, valid)
        ! Flips the spins between the given site and the next
        ! periodic boundary conditions are used
        ! valid is false if operator gives the empty state
        integer,intent(in) :: site
        integer,intent(in) :: BasisFunc
        integer,intent(out) :: NewBasisFunc
        logical,intent(out) :: valid
        integer :: SiteNext 
        valid  = .true.
        SiteNext = MOD(Site + 1, NumSites)
        NewBasisFunc = BasisFunc
        if (btest(BasisFunc,Site) .and. .not. btest(BasisFunc,SiteNext) )then
            NewBasisFunc =  ibclr(NewBasisFunc,Site)
            NewBasisFunc = ibset(NewBasisFunc,SiteNext)
            return
        endif
        if (.not. btest(BasisFunc,Site) .and.  btest(BasisFunc,SiteNext) )then
            NewBasisFunc =  ibset(NewBasisFunc,Site)
            NewBasisFunc = ibclr(NewBasisFunc,SiteNext)
            return
        endif
        valid = .false.

    end subroutine OperatorSpinFlip

    
    subroutine OperatorSzSzNext(base, weight)
        integer, intent(in) :: base
        real(dp), intent(out) :: weight

        integer j
        weight  = 0
        do j = 0,NumSites-2
            
            weight = weight + (IBITS(base,j,1) - 0.5) * (IBITS(base, MOD(j + 1, NumSites),1)  - 0.5)
       enddo

    end subroutine
 
    subroutine InitSpinSystem(sites)
        integer :: sites

        NumSites = sites
        BasisDimension = 2**NumSites

        if( ALLOCATED(Basis)) then
             DEALLOCATE(Basis)
             DEALLOCATE(BasisLookup)
        endif
        ALLOCATE(Basis(BasisDimension))
        ALLOCATE(BasisLookup(0:BasisDimension-1))
        
        ! generate basis
        call GenBasis
    end subroutine InitSpinSystem

    subroutine OperateSz(state, site)
        real(dp), intent(inout) :: state(:)
        integer :: site, i,N
        integer :: basis_element
        N = size(state)
        do i=1,N
            basis_element = Basis(i)
            ! print *, i, basis_element, ibits(basis_element,site,1)
            state(i) = state(i) * (ibits(basis_element,site,1) - 0.5)
        enddo
        

    end subroutine

    subroutine GenBasis
        implicit none
        integer :: i, j
        real :: Sz
        integer :: count 
        count = 0
        Basis = 0
        BasisLookup = -1 ! not present
        
        do i = 0, BasisDimension-1
            Sz = 0
            do j = 0, NumSites-1
                Sz = Sz +  (IBITS(i,j,1) - 0.5)
            enddo
            if (Sz == 0) then
                 count = count  + 1
                 Basis(count) = i
                 BasisLookup(i) = count
            endif
        enddo
        BasisDimension = count
        ! eliminate zeros
    end subroutine GenBasis

end module spin_system