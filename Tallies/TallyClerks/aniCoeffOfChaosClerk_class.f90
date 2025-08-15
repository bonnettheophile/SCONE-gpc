module aniCoeffOfChaosClerk_class

    use numPrecision
    use tallyCodes
    use universalVariables
    use dictionary_class,      only : dictionary
    use genericProcedures,     only : fatalError
    use particle_class,        only : particle, particleState
    use particleDungeon_class, only : particleDungeon
    use outputFile_class,      only : outputFile
    use legendrePoly_func,     only : evaluateLegendre
  
    use scoreMemory_class,     only : scoreMemory
    use tallyResult_class,     only : tallyResult, tallyResultEmpty, polyResult
    use tallyClerk_inter,      only : tallyClerk, kill_super => kill

    ! Tally Maps
    use tallyMap_inter,             only : tallyMap
    use tallyMapFactory_func,       only : new_tallyMap

    use linearAlgebra_func
  
    implicit none
    private
    
    !! 
    !! Estimator for polynomial chaos coefficients when using virtual densities
    !! Output is a 2D array of size (order+1) ** number. We start assuming a single random parameter
    !! Additional assumption: population control is enforcing starting weight is constant

    !! Private members:
    !!   chaosOfPop   -> Array for end of generation chaos model
    !!   map          -> Space to store tally map
    !!   P            -> Order of chaos model
    !!   startPop     -> Starting population for each cycle

    !! Interface
    !!   tallyClerk Interface

    !! SAMPLE DICTIONARY INPUT
    !! 
    !! myaniCoeffOfChaosClerk {
    !!    type aniCoeffOfChaosClerk;
    !!    order 2;
    !!    # map { <tallyMap definition>} #
    !! }
    !!

    type, public, extends(tallyClerk) :: aniCoeffOfChaosClerk
      private
        real(defReal), dimension(:,:,:), allocatable :: chaosOfPop   ! Coefficients for the end of generation pop
        class(tallyMap), allocatable                 :: map
        integer(shortInt)                            :: P            ! Order of gpc model
        integer(shortInt)                            :: fitOrder
        real(defReal)                                :: startPop
        real(defReal), allocatable                   :: histogram(:)
        real(defReal), allocatable                   :: binCentre(:)
        real(defReal), allocatable                   :: fitCoeff(:,:)
        real(defReal)                                :: dx
        logical(defBool)                             :: firstCycle = .true.

      
    contains
      ! Procedures used during build
      procedure :: init
      procedure :: kill
      procedure :: validReports
      procedure :: getSize

      procedure :: reportCycleStart
      procedure :: reportCycleEnd

      procedure :: print
      procedure :: display
      procedure :: getResult
      procedure :: computeChaosPop
    end type aniCoeffOfChaosClerk

contains

    !! Initialise from dictionary and name
    subroutine init(self, dict, name)
      class(aniCoeffOfChaosClerk), intent(inout) :: self
      class(dictionary), intent(in)           :: dict
      character(nameLen), intent(in)          :: name
      integer(shortInt)                       :: i
      character(100),parameter :: Here = 'init (aniCoeffOfChaosClerk.f90)'

      ! Needs no settings, just load name
      call self % setName(name)
    
      self % startPop = ZERO

      ! Order keyword must be present
      if (dict % isPresent('order')) then
        call dict % get(self % P, ' order')
        allocate(self % chaosOfPop(self % P + 1,self % P + 1,self % P + 1))
        self % chaosOfPop = ZERO
      else 
        call fatalError(Here, "Order must by provided") 
      end if

      ! fitOrder keyword must be present
      if (dict % isPresent('fitOrder')) then
        call dict % get(self % fitOrder, ' fitOrder')
        allocate(self % fitCoeff(3, self % fitOrder + 1))
        self % fitCoeff = ZERO
      else 
        call fatalError(Here, "fitOrder must by provided") 
      end if

      ! Map is for following gpc coefficients
      if( dict % isPresent('map')) then
        call new_tallyMap(self % map, dict % getDictPtr('map'))
        allocate(self % histogram(product(self % map % binArrayShape())))
        allocate(self % binCentre(product(self % map % binArrayShape())))


        self % histogram = ZERO

        self % dx = TWO / size(self % histogram)
        self % binCentre(1) = - ONE + self % dx / TWO

        ! Initialize binCentre, assume interval is [-1,1]
        do i = 2, size(self % binCentre)
          self % binCentre(i) = self % binCentre(i-1) + self % dx 
        end do
      end if

    end subroutine init

    !! Return to uninitialised state
    elemental subroutine kill(self)
      class(aniCoeffOfChaosClerk), intent(inout) :: self
      
      ! Kill superclass
      call kill_super(self)

      if (allocated(self % chaosOfPop)) deallocate(self % chaosOfPop)
      if (allocated(self % fitCoeff)) deallocate(self % fitCoeff)
      if (allocated(self % map)) then
        call self % map % kill()
        deallocate(self % map)
      end if
      if (allocated(self % histogram)) deallocate(self % histogram)
      if (allocated(self % binCentre)) deallocate(self % binCentre)    
    end subroutine kill


    ! Only needed at start/end of generation
    function validReports(self) result(validCodes)
      class(aniCoeffOfChaosClerk), intent(in)          :: self
      integer(shortInt), dimension(:), allocatable  :: validCodes

      validCodes = [ cycleStart_CODE, cycleEnd_CODE ]
    end function validReports
    
    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!

    elemental function getSize(self) result(S)
      class(aniCoeffOfChaosClerk), intent(in) :: self
      integer(shortInt)                    :: S
    
      ! Order to the power of number of random parameters
      S = (self % P + 1)**3
    end function getSize

    !! Process beginning of a cycle
    !! While this allows non-constant generation weight, it should be ensured 
    !! for this implementation to give sensible results.
    subroutine reportCycleStart(self, start, mem)
      class(aniCoeffOfChaosClerk), intent(inout) :: self
      class(particleDungeon), intent(in)      :: start
      type(scoreMemory), intent(inout)        :: mem

      self % startPop = start % popWeight()
    end subroutine reportCycleStart

    !! Process end of cycle and compute chaotic coefficients 
    !! Assuming we have a single parameter
    !! 
    subroutine reportCycleEnd(self, end, mem)
      class(aniCoeffOfChaosClerk), intent(inout) :: self
      class(particleDungeon), intent(in)         :: end
      type(scoreMemory), intent(inout)           :: mem
      real(defReal), dimension(3, self % P + 1)  :: legendrePol
      real(defReal)                              :: tmp_score(self % P + 1, self % P + 1, self % P + 1)
      real(defReal), allocatable                 :: gaussPoints(:,:)
      real(defReal)                              :: chaosPop
      type(particle)                             :: p
      type(particleState)                        :: state
      integer(shortInt)                          :: i, l, binIdx
      integer(shortInt)                          :: j1, j2, j3, k1, k2, k3, G, O
      real(defReal)                              :: gaussWgt, polProduct, score
      integer(longInt)                           :: address
      real(defReal), dimension(size(self % binCentre))                      :: x, b
      real(defReal), dimension(size(self % binCentre), self % fitOrder+1)   :: A 
      character(100),parameter :: Here = 'reportCycleEnd (aniCoeffOfChaosClerk.f90)'
      
      O = self % P
      ! Get adequate quadrature parameters
      if (O == 1) then
        gaussPoints = G2
        G = 2
      else if (O == 2) then
        gaussPoints = G3
        G = 3
      else if (O == 3) then
        gaussPoints = G4
        G = 4
      else if (O == 4) then
        gaussPoints = G5
        G = 5
      else
        call fatalError(Here, "Gauss quadrature order not supported")
      end if

      ! Reinialize histogram array
      do l = 1, 3
        self % histogram = ZERO
        do i = 1, end % popSize()
          p = end % get(i)
          state = p 

          ! Cheating to not have to struggle with multiple maps indexing
          state % X(1) = p % X(l)

          ! Find bin index
          if (allocated(self % map)) then
            binIdx = self % map % map(state)
          else
            binIdx = 1
          end if
          ! Return if invalid bin index
          if (binIdx == 0) return
          
          ! Fill histogram of particles wrt their uncertain parameter
          self % histogram(binIdx) = self % histogram(binIdx) + state % wgt
        end do

        ! Set x array for linear fitting
        do i = 1, size(A, 2)
          A(:,i) = self % binCentre**(i-1)
        end do 

        ! Set y for linear fitting
        b = self % histogram
        ! Perform least square linear fitting using LAPACK 
        call solveLeastSquare(A, x, b)
        ! Save fit results
        !self % fitCoeff = self % fitCoeff + (x(1:self % fitOrder+1) - self % fitCoeff) / (mem % cycles+1)
        self % fitCoeff(l,:) = x(1:self % fitOrder+1)
      end do
      call kill_linearAlgebra()

      ! Initialise temporary score
      tmp_score = ZERO
      ! First we need to build the gpc model for the end population
      do i = 1, end % popSize()
        ! Reinitialise temporary score
        ! Evaluate Legendre polynomials up to right order 
        p = end % get(i)
        legendrePol(1,:) = evaluateLegendre(O, p % X(1)) 
        legendrePol(2,:) = evaluateLegendre(O, p % X(2)) 
        legendrePol(3,:) = evaluateLegendre(O, p % X(3)) 
        !legendrePol = evaluateLegendre(O, ZERO)
        do j3 = 1, O + 1
          do j2 = 1, O + 1
            do j1 = 1, O + 1
              tmp_score(j1,j2,j3) = tmp_score(j1,j2,j3) + (TWO*(real(j1,defReal)-ONE) + ONE) &
                        * (TWO*(real(j2,defReal)-ONE) + ONE) * (TWO*(real(j3,defReal)-ONE) + ONE) &
                        * legendrePol(1,j1) * legendrePol(2,j2) * legendrePol(3,j3) * p % w * p % k_eff
            end do
          end do
        end do
      end do

      ! Update chaos model for new population size
      self % chaosOfPop = self % chaosOfPop + (tmp_score - self % chaosOfPop) / (mem % cycles+1)

      ! Reinitialise temporary score
      tmp_score = ZERO
      ! Use chaos model to determine population at quadrature points

      ! Loop on gpc order
      do j3 = 1, O + 1
        do j2 = 1, O + 1
          do j1 = 1, O + 1
            ! Loop for 3d gauss quadrature
            score = ZERO
            do k3 = 1, G
              legendrePol(1,:) = evaluateLegendre(O, gaussPoints(1, k3))
              do k2 = 1, G
                legendrePol(2,:) = evaluateLegendre(O, gaussPoints(1, k2))
                do k1 = 1, G
                  legendrePol(3,:) = evaluateLegendre(O, gaussPoints(1, k1))
                  gaussWgt = gaussPoints(2, k1) * gaussPoints(2, k2) * gaussPoints(2, k3) / (TWO**3)
                  polProduct = legendrePol(1,j1) * legendrePol(2,j2) * legendrePol(3,j3)
                  chaosPop = self % computeChaosPop(gaussPoints(1, k1), gaussPoints(1, k2), gaussPoints(1, k3))
                  score = score + (TWO*(real(j1,defReal)-ONE) + ONE)*(TWO*(real(j2,defReal)-ONE) + ONE) &
                            *(TWO*(real(j3,defReal)-ONE) + ONE) * polProduct * gaussWgt &
                            * chaosPop / self % startPop
                end do
              end do
            end do
            tmp_score(j1,j2,j3) = tmp_score(j1,j2,j3) + score
          end do
        end do
      end do

      ! Score coefficients for chaotic keff
      do j3 = 1, O + 1
        do j2 = 1, O + 1
          do j1 = 1, O + 1
            address = self % getMemAddress() + j1 + (j2-1) * (O+1) + (j3-1) * (O+1)**2 - 1
            call mem % accumulate(tmp_score(j1, j2, j3), address)
          end do
        end do
      end do 

      deallocate(gaussPoints)
    
    end subroutine reportCycleEnd

  !!
  !! Evaluate the gpc model for the end of generation population size
  !! Input is the vector at which it needs to be evaluated
  !! Return the population size
  !!

  function computeChaosPop(self, x, y, z) result(pop)
    class(aniCoeffOfChaosClerk), intent(in)   :: self
    real(defReal), intent(in)                 :: x, y, z
    real(defReal)                             :: pop
    real(defReal), dimension(3, self % P + 1) :: legendrePol
    integer(shortInt)                         :: i, j, k, P

    P = self % P
    pop = ZERO
    legendrePol(1,:) = evaluateLegendre(P, x)
    legendrePol(2,:) = evaluateLegendre(P, y)
    legendrePol(3,:) = evaluateLegendre(P, z)
    do k = 1, P + 1
      do j = 1, P + 1
        do i = 1, P + 1
          pop = pop + self % chaosOfPop(i,j,k) * legendrePol(1,i) * legendrePol(2, j) * legendrePol(3, k)
        end do
      end do
    end do
  end function computeChaosPop

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
    subroutine display(self, mem)
        class(aniCoeffOfChaosClerk), intent(in)  :: self
        type(scoreMemory), intent(in)      :: mem
    
        print *, 'aniCoeffOfChaosClerk does not support display yet'
    
      end subroutine display
    
  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
      subroutine print(self, outFile, mem)
        class(aniCoeffOfChaosClerk), intent(in)       :: self
        class(outputFile), intent(inout)           :: outFile
        type(scoreMemory), intent(in)              :: mem
        real(defReal)                              :: val, std
        integer(shortInt)                          :: j1, j2, j3, O
        integer(shortInt),dimension(:),allocatable :: resArrayShape
        integer(longInt)                           :: address
        character(nameLen)                         :: name
        
        O = self % P

        ! Begin block
        call outFile % startBlock(self % getName())
    
        ! If collision clerk has map print map information
        if (allocated(self % map)) then
          call self % map % print(outFile)
        end if
    
        ! Write results.
        ! Get shape of result array
          resArrayShape = [self % getSize()]
    
        ! Start array
        name ='Res'
        call outFile % startArray(name, resArrayShape)
    
        ! Print results to the file
        !do i = 1, product(resArrayShape)
        !  call mem % getResult(val, std, self % getMemAddress() - 1 + i)
        !  call outFile % addResult(val,std)
        !end do

        do j3 = 1, O + 1
          do j2 = 1, O + 1
            do j1 = 1, O + 1
              address = self % getMemAddress() + j3 + (j2-1) * (O+1) + (j1-1) * (O+1)**2 - 1
              call mem % getResult(val, std, address)
              call outFile % addResult(val,std)
            end do
          end do
        end do
    
        call outFile % endArray()
        call outFile % endBlock()
    
      end subroutine print

  !! 
  !! Get fitting model for end of generation probability distribution of X
  !!

    pure subroutine getResult(self, res, mem)
      class(aniCoeffOfChaosClerk), intent(in)              :: self
      class(tallyResult), allocatable, intent(inout)    :: res
      type(scoreMemory), intent(in)                     :: mem

      allocate(res, source= polyResult(self % fitCoeff))

    end subroutine
    
    end module aniCoeffOfChaosClerk_class    



