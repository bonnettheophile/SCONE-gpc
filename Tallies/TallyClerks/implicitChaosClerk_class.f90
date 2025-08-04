module implicitChaosClerk_class

    use numPrecision
    use tallyCodes
    use universalVariables
    use endfConstants
    use dictionary_class,      only : dictionary
    use genericProcedures,     only : fatalError
    use particle_class,        only : particle, particleState
    use particleDungeon_class, only : particleDungeon
    use outputFile_class,      only : outputFile
    use legendrePoly_func,     only : evaluateLegendre
    use genericProcedures,     only : binarySearch
  
  ! Nuclear Data Interfaces
    use nuclearDataReg_mod,         only : ndReg_get => get
    use nuclearDatabase_inter,      only : nuclearDatabase
    use neutronMaterial_inter,      only : neutronMaterial,neutronMaterial_CptrCast
    use neutronXSPackages_class,    only : neutronMacroXSs

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
    !!    type implicitChaosClerk;
    !!    order 2;
    !!    # map { <tallyMap definition>} #
    !! }
    !!

    type, public, extends(tallyClerk) :: implicitChaosClerk
      private
        real(defReal), dimension(:), allocatable :: chaosOfPop   ! Coefficients for the end of generation pop
        integer(shortInt)                        :: P            ! Order of gpc model
        class(tallyMap), allocatable             :: map
        real(defReal)                        :: startPop
        real(defReal)                        :: a, b, norm
        real(defReal), dimension(:), allocatable :: histogram, binCentre, binSup, cumLaw, binEdges
      
    contains
      ! Procedures used during build
      procedure :: init
      procedure :: kill
      procedure :: validReports
      procedure :: getSize

      procedure :: reportCycleStart
      procedure :: reportHist

      procedure :: print
      procedure :: display

    end type implicitChaosClerk

contains

    !! Initialise from dictionary and name
    subroutine init(self, dict, name)
      class(implicitChaosClerk), intent(inout) :: self
      class(dictionary), intent(in)           :: dict
      character(nameLen), intent(in)          :: name
      character(nameLen)                      :: chr
      real(defReal)                           :: dx
      integer(shortInt)                       :: i
      character(100),parameter :: Here = 'init (implicitChaosClerk.f90)'

      ! Needs no settings, just load name
      call self % setName(name)
    
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
        allocate(self % binSup(product(self % map %binArrayShape())))
        allocate(self % binEdges(product(self % map %binArrayShape())+1))
        allocate(self % cumLaw(product(self % map %binArrayShape())))

        self % histogram = ZERO
        self % cumLaw = ZERO

        dx = TWO / size(self % histogram)

        self % binEdges(1) = - ONE
        ! Initialize binCentre, assume interval is [-1,1]
        do i = 1, size(self % binCentre)
          self % binCentre(i) = - ONE + i*dx/TWO 
          self % binSup(i)    = - ONE + i*dx
          self % binEdges(i+1) = - ONE + i*dx
        end do
      end if


      self % startPop = ZERO
    end subroutine init

    !! Return to uninitialised state
    elemental subroutine kill(self)
      class(implicitChaosClerk), intent(inout) :: self
      
      ! Kill superclass
      call kill_super(self)

      if (allocated(self % chaosOfPop)) deallocate(self % chaosOfPop)
      if (allocated(self % map)) then
        call self % map % kill()
        deallocate(self % map)
      end if
      if (allocated(self % histogram)) deallocate(self % histogram)
      if (allocated(self % binCentre)) deallocate(self % binCentre)
      if (allocated(self % binEdges)) deallocate(self % binEdges)
      if (allocated(self % binSup)) deallocate(self % binSup)

    end subroutine kill


    ! Only needed at start/end of generation
    function validReports(self) result(validCodes)
      class(implicitChaosClerk), intent(in)          :: self
      integer(shortInt), dimension(:), allocatable  :: validCodes

      validCodes = [hist_CODE, cycleStart_CODE]
    end function validReports
    
    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!

    elemental function getSize(self) result(S)
      class(implicitChaosClerk), intent(in) :: self
      integer(shortInt)                    :: S
    
      S = self % P + 1
    end function getSize

    subroutine reportCycleStart(self, start, mem)
      class(implicitChaosClerk), intent(inout) :: self
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

      !print *, self % cumLaw

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

      print *, self % a, self % b
      ! Normalize linear law
      !self % a = self % a / (2 * self % b)


    
    end subroutine reportCycleStart
  
  
    !!
    !! Process history report
    !! Gets fate code from the particle
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportHist(self, p, xsData, mem)
      class(implicitChaosClerk), intent(inout) :: self
      class(particle), intent(in)             :: p
      class(nuclearDatabase),intent(inout)    :: xsData
      type(scoreMemory), intent(inout)        :: mem
      real(defReal)                           :: score, val
      real(defReal), dimension(self % P + 1)  :: legendrePol
      integer(shortInt)                       :: j


      if (p % fate /= leak_FATE) then
        ! Evaluate Legendre polynomials up to right order
        !val = ((self % a * p % X(1)**2) / TWO + self % b * p % X(1) + self % b**2 / (2 * self % a))
        val = self % cumLaw(binarySearch(self % binEdges, p % X(1)))
        !do j = 1, size(self % cumLaw) 
        !  if (p % X(1) <= self % binSup(j)) then 
        !    val = self % cumLaw(j)
        !    exit
        !  end if
        !end do
        val = 2*val - ONE
        if (abs(val) > 1) then
          !print *, val
          if (val > ONE) val = ONE
          if (val < - ONE) val = - ONE
        end if
        
        legendrePol = evaluateLegendre(self % P, val) 
        do j = 1, self % P + 1
          score = (2*(j-1) + 1) * legendrePol(j) * p % w / self % startPop
          call mem % score(score, self % getMemAddress() + j - 1)
        end do
      end if
    end subroutine reportHist

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
    subroutine display(self, mem)
        class(implicitChaosClerk), intent(in)  :: self
        type(scoreMemory), intent(in)      :: mem
    
        print *, 'implicitChaosClerk does not support display yet'
    
      end subroutine display
    
  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
      subroutine print(self, outFile, mem)
        class(implicitChaosClerk), intent(in)       :: self
        class(outputFile), intent(inout)           :: outFile
        type(scoreMemory), intent(in)              :: mem
        real(defReal)                              :: val, std
        integer(shortInt)                          :: i
        integer(shortInt),dimension(:),allocatable :: resArrayShape
        character(nameLen)                         :: name
    
        ! Begin block
        call outFile % startBlock(self % getName())
    
        ! Write results.
        ! Get shape of result array
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
    
    end module implicitChaosClerk_class    



