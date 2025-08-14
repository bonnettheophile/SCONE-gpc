module coeffOfChaosClerk_class

    use numPrecision
    use tallyCodes
    use universalVariables
    use dictionary_class,      only : dictionary
    use genericProcedures,     only : fatalError
    use particle_class,        only : particle, particleState
    use particleDungeon_class, only : particleDungeon
    use outputFile_class,      only : outputFile
    use legendrePoly_func,     only : evaluateLegendre
    use polynomial_class,      only : polynomial
    use genericProcedures,     only : binarySearch, interpolate, quickSort

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
    !! myCoeffofChaosClerk {
    !!    type coeffOfChaosClerk;
    !!    order 2;
    !!    # map { <tallyMap definition>} #
    !! }
    !!

    type, public, extends(tallyClerk) :: coeffOfChaosClerk
      private
        real(defReal), dimension(:), allocatable :: chaosOfPop, lastChaosOfPop   ! Coefficients for the end of generation pop
        class(tallyMap), allocatable             :: map
        integer(shortInt)                        :: P            ! Order of gpc model
        integer(shortInt)                        :: fitOrder     ! Order of polyfit for pdf
        real(defReal)                            :: startPop
        real(defReal), allocatable, dimension(:) :: histogram, binEdges, binCentre, cumLaw, values, fitCoeff
        real(defReal)                            :: dx
        logical(defBool)                         :: firstCycle = .true.


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
    end type coeffOfChaosClerk

contains

    !! Initialise from dictionary and name
    subroutine init(self, dict, name)
      class(coeffOfChaosClerk), intent(inout) :: self
      class(dictionary), intent(in)           :: dict
      character(nameLen), intent(in)          :: name
      integer(shortInt)                       :: i  
      character(100),parameter :: Here = 'init (coeffOfChaosClerk.f90)'

      ! Needs no settings, just load name
      call self % setName(name)
    
      self % startPop = ZERO

      ! Order keyword must be present
      if (dict % isPresent('order')) then
        call dict % get(self % P, ' order')
        allocate(self % chaosOfPop(self % P + 1))
        allocate(self % lastChaosOfPop(self % P + 1))
        self % chaosOfPop = ZERO
        self % lastChaosOfPop = ZERO
      else 
        call fatalError(Here, "Order must by provided") 
      end if

      ! Order keyword must be present
      if (dict % isPresent('fitOrder')) then
        call dict % get(self % fitOrder, ' fitOrder')
        allocate(self % fitCoeff(self % fitOrder + 1))
        self % fitCoeff = ZERO
      else 
        call fatalError(Here, "fitOrder must by provided") 
      end if

      ! Map is for following gpc coefficients
      if( dict % isPresent('map')) then
        call new_tallyMap(self % map, dict % getDictPtr('map'))
        allocate(self % histogram(product(self % map %binArrayShape())))
        allocate(self % binCentre(product(self % map %binArrayShape())))
        allocate(self % binEdges(product(self % map %binArrayShape())+1))
        allocate(self % cumLaw(product(self % map %binArrayShape())))

        self % histogram = ZERO
        self % cumLaw = ZERO

        self % dx = TWO / size(self % histogram)
        self % binEdges(1) = - ONE
        self % binCentre(1) = - ONE + self % dx / TWO

        ! Initialize binCentre, assume interval is [-1,1]
        do i = 2, size(self % binCentre)
          self % binCentre(i) = self % binCentre(i-1) + self % dx 
        end do
        do i = 2, size(self % binEdges)
          self % binEdges(i) = - ONE + i * self % dx
        end do

      end if


    end subroutine init

    !! Return to uninitialised state
    elemental subroutine kill(self)
      class(coeffOfChaosClerk), intent(inout) :: self
      
      ! Kill superclass
      call kill_super(self)

      if (allocated(self % chaosOfPop)) deallocate(self % chaosOfPop)
      if (allocated(self % lastChaosOfPop)) deallocate(self % lastChaosOfPop)
      if (allocated(self % fitCoeff)) deallocate(self % fitCoeff)
      if (allocated(self % map)) then
        call self % map % kill()
        deallocate(self % map)
      end if
      if (allocated(self % histogram)) deallocate(self % histogram)
      if (allocated(self % values)) deallocate(self % values)
      if (allocated(self % binEdges)) deallocate(self % binEdges)
      if (allocated(self % cumLaw)) deallocate(self % cumLaw)

      self % fitOrder = 0
      self % P = 0

    end subroutine kill


    ! Only needed at start/end of generation
    function validReports(self) result(validCodes)
      class(coeffOfChaosClerk), intent(in)          :: self
      integer(shortInt), dimension(:), allocatable  :: validCodes

      validCodes = [ cycleStart_CODE, cycleEnd_CODE ]
    end function validReports
    
    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!

    elemental function getSize(self) result(S)
      class(coeffOfChaosClerk), intent(in) :: self
      integer(shortInt)                    :: S
    
      S = self % P + 1
      !if (allocated(self % map)) S = self % map % bins(0)
    end function getSize

    !! Process beginning of a cycle
    !! While this allows non-constant generation weight, it should be ensured 
    !! for this implementation to give sensible results.
    subroutine reportCycleStart(self, start, mem)
      class(coeffOfChaosClerk), intent(inout)                             :: self
      class(particleDungeon), intent(in)                                  :: start
      type(scoreMemory), intent(inout)                                    :: mem

      self % startPop = start % popWeight()

    end subroutine reportCycleStart

    !! Process end of cycle and compute chaotic coefficients 
    !! Assuming we have a single parameter
    !! 
    subroutine reportCycleEnd(self, end, mem)
      class(coeffOfChaosClerk), intent(inout)     :: self
      class(particleDungeon), intent(in)          :: end
      type(scoreMemory), intent(inout)            :: mem
      real(defReal), dimension(self % P + 1)      :: legendrePol, tmp_score
      real(defReal), dimension(:,:), allocatable  :: gaussPoints
      real(defReal)                               :: chaosPop, lastChaosPop
      type(particle)                              :: p
      type(particleState)                         :: state
      integer(shortInt)                           :: i, j, G, binIdx
      real(defReal), dimension(size(self % binCentre))                      :: x, b
      real(defReal), dimension(size(self % binCentre), self % fitOrder+1)   :: A 
      character(100),parameter :: Here = 'reportCycleEnd (coeffOfChaosClerk.f90)'

      ! Get adequate quadrature parameters
      if (self % P == 1) then
        gaussPoints = G2
        G = 2
      else if (self % P == 2) then
        gaussPoints = G3
        G = 3
      else if (self % P == 3) then
        gaussPoints = G4
        G = 4
      else if (self % P == 4) then
        gaussPoints = G5
        G = 5
      else
        call fatalError(Here, "Gauss quadrature order not supported")
      end if

      ! Reinialize histogram array
      self % histogram = ZERO
      do i = 1, end % popSize()
        p = end % get(i)
        state = p 

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
      self % fitCoeff = self % fitCoeff + (x(1:self % fitOrder+1) - self % fitCoeff) / (mem % cycles+1)
      call kill_linearAlgebra()

      ! Initialise temporary score
      tmp_score = ZERO
      ! First we need to build the gpc model for the end population
      do i = 1, end % popSize()
        p = end % get(i)
        ! Evaluate Legendre polynomials up to right order 
        legendrePol = evaluateLegendre(self % P, p % X(1)) 
        do j = 1, self % P + 1
          tmp_score(j) = tmp_score(j) + (2*(j-1) + 1) * legendrePol(j) * p % w * p % k_eff
        end do
      end do

      ! Update chaos model for new population size
      self % chaosOfPop = self % chaosOfPop + (tmp_score - self % chaosOfPop) / (mem % cycles+1)

      ! Reinitialise temporary score
      tmp_score = ZERO
      ! Use chaos model to determine population at quadrature points
      do i = 1, G
        legendrePol = evaluateLegendre(self % P, gaussPoints(1, i))
        chaosPop = sum(legendrePol * self % chaosOfPop)
        lastChaosPop = self % startPop
        do j = 1, self % P + 1
          tmp_score(j) = tmp_score(j) + &
               (2*(j-1) + 1) * legendrePol(j) * chaosPop / lastChaosPop * gaussPoints(2,i) / TWO
        end do
      end do

      ! Score coefficients for chaotic keff
      do i = 1, self % P + 1
        call mem % accumulate(tmp_score(i), self % getMemAddress() + i - 1)
      end do 
      deallocate(gaussPoints)
    end subroutine reportCycleEnd

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
    subroutine display(self, mem)
        class(coeffOfChaosClerk), intent(in)  :: self
        type(scoreMemory), intent(in)      :: mem
    
        print *, 'coeffOfChaosClerk does not support display yet'
    
      end subroutine display
    
  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
      subroutine print(self, outFile, mem)
        class(coeffOfChaosClerk), intent(in)       :: self
        class(outputFile), intent(inout)           :: outFile
        type(scoreMemory), intent(in)              :: mem
        real(defReal)                              :: val, std
        integer(shortInt)                          :: i
        integer(shortInt),dimension(:),allocatable :: resArrayShape
        character(nameLen)                         :: name
    
        ! Begin block
        call outFile % startBlock(self % getName())
    
        resArrayShape = [self % getSize()]
    
        ! Start array
        name ='Res'
        call outFile % startArray(name, resArrayShape)
    
        ! Print results to the file
        do i = 1, product(resArrayShape)
          call mem % getResult(val, std, self % getMemAddress() - 1 + i)
          call outFile % addResult(val,std)
    
        end do
    
        call outFile % endArray()
        call outFile % endBlock()
    
      end subroutine print
    
  !! 
  !! Get fitting model for end of generation probability distribution of X
  !!

    pure subroutine getResult(self, res, mem)
      class(coeffOfChaosClerk), intent(in)              :: self
      class(tallyResult), allocatable, intent(inout)    :: res
      type(scoreMemory), intent(in)                     :: mem

      allocate(res, source= polyResult(self % fitCoeff))

    end subroutine
      
  end module coeffOfChaosClerk_class    



