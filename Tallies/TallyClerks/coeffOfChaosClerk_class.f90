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
    use genericProcedures,     only : binarySearch
  
    use scoreMemory_class,     only : scoreMemory
    use tallyResult_class,     only : tallyResult, tallyResultEmpty
    use tallyClerk_inter,      only : tallyClerk, kill_super => kill

    ! Tally Maps
    use tallyMap_inter,             only : tallyMap
    use tallyMapFactory_func,       only : new_tallyMap
  
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
        real(defReal), dimension(:), allocatable :: chaosOfPop   ! Coefficients for the end of generation pop
        class(tallyMap), allocatable             :: map
        integer(shortInt)                        :: P            ! Order of gpc model
        real(defReal)                            :: startPop
        real(defReal), allocatable, dimension(:) :: histogram, binEdges, binCentre, cumLaw
        real(defReal)                            :: a, b, norm

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
    end type coeffOfChaosClerk

contains

    !! Initialise from dictionary and name
    subroutine init(self, dict, name)
      class(coeffOfChaosClerk), intent(inout) :: self
      class(dictionary), intent(in)           :: dict
      character(nameLen), intent(in)          :: name
      real(defReal)                           :: dx
      integer(shortInt)                       :: i  
      character(100),parameter :: Here = 'init (coeffOfChaosClerk.f90)'

      ! Needs no settings, just load name
      call self % setName(name)
    
      self % startPop = ZERO

      ! Order keyword must be present
      if (dict % isPresent('order')) then
        call dict % get(self % P, ' order')
        allocate(self % chaosOfPop(self % P + 1))
        self % chaosOfPop = ZERO
      else 
        call fatalError(Here, "Order must by provided") 
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

        dx = TWO / size(self % histogram)
        self % binEdges(1) = - ONE

        ! Initialize binCentre, assume interval is [-1,1]
        do i = 1, size(self % binEdges)
          self % binEdges(i+1) = - ONE + i*dx 
          self % binCentre(i) = - ONE + i*dx/TWO 
        end do
      end if


    end subroutine init

    !! Return to uninitialised state
    elemental subroutine kill(self)
      class(coeffOfChaosClerk), intent(inout) :: self
      
      ! Kill superclass
      call kill_super(self)

      if (allocated(self % chaosOfPop)) deallocate(self % chaosOfPop)
      if (allocated(self % map)) then
        call self % map % kill()
        deallocate(self % map)
      end if
      if (allocated(self % histogram)) deallocate(self % histogram)
      if (allocated(self % binEdges)) deallocate(self % binEdges)
      if (allocated(self % cumLaw)) deallocate(self % cumLaw)

      self % a = ZERO
      self % b = ZERO
      self % norm = ZERO
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
      class(coeffOfChaosClerk), intent(inout) :: self
      class(particleDungeon), intent(in)      :: start
      type(scoreMemory), intent(inout)        :: mem
      real(defReal)                         :: Sw, Sx, Sy, Sxx, Sxy, S
      type(particleState)                   :: state
      type(particle)                        :: p
      integer(shortInt) :: binIdx, i

      self % startPop = start % popWeight()

      self % histogram = ZERO
      self % cumLaw = ZERO
      S = ZERO
      
      do i = 1, start % popSize()
        p = start % get(i)
        state = p 

        ! Find bin index
        if (allocated(self % map)) then
          binIdx = self % map % map(state)
        else
          binIdx = 1
        end if
        ! Return if invalid bin index
        if (binIdx == 0) return
        
        self % histogram(binIdx) = self % histogram(binIdx) + state % wgt
        S = S + state % wgt
      end do

      self % cumLaw(1) = self % histogram(1)
      do i = 2, size(self % cumLaw)
        self % cumLaw(i) = self % cumLaw(i-1) + self % histogram(i)
      end do

      self % cumLaw = self % cumLaw / S

      ! Do weighted least square linear fitting on the histogram
      Sw = sum(self % histogram)
      Sx = sum(self % histogram * self % binCentre)
      Sy = sum(self % histogram**2)
      Sxy = sum(self % histogram**2 * self % binCentre)
      Sxx = sum(self % binCentre**2 * self % histogram)

      ! Compute linear fit coefficients
      self % a = (Sw * Sxy - Sx * Sy) / (Sw * Sxx - Sx**2)
      self % b = (Sy - self % a * Sx) / Sw
      self % norm = self % a / TWO + self % b + self % b**2 / ( 2 * self % a)

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
      real(defReal)                               :: chaosPop = ZERO, val
      type(particle)                              :: p
      integer(shortInt)                           :: i, j, G
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

      ! Initialise temporary score
      tmp_score = ZERO
      ! First we need to build the gpc model for the end population
      do i = 1, end % popSize()
        ! Reinitialise temporary score
        p = end % get(i)
        ! Remap from start of gen distribution to U[0,1)
        val = self % cumLaw(binarySearch(self % binEdges, p % X(1)))
        ! Remap from U[0,1) to U[-1,1)
        val = 2*val - ONE
        ! Evaluate Legendre polynomials up to right order 
        legendrePol = evaluateLegendre(self % P, val) 

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
        do j = 1, self % P + 1
          tmp_score(j) = tmp_score(j) + &
               (2*(j-1) + 1) * legendrePol(j) * chaosPop / self % startPop * gaussPoints(2,i)/TWO
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
    
    end module coeffOfChaosClerk_class    



