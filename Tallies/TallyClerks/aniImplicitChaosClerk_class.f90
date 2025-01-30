module aniImplicitChaosClerk_class

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
    !!    type aniImplicitChaosClerk;
    !!    order 2;
    !!    # map { <tallyMap definition>} #
    !! }
    !!

    type, public, extends(tallyClerk) :: aniImplicitChaosClerk
      private
        real(defReal), dimension(:), allocatable :: chaosOfPop   ! Coefficients for the end of generation pop
        integer(shortInt)                        :: P            ! Order of gpc model
        real(defReal)                        :: startPop
      
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

    end type aniImplicitChaosClerk

contains

    !! Initialise from dictionary and name
    subroutine init(self, dict, name)
      class(aniImplicitChaosClerk), intent(inout) :: self
      class(dictionary), intent(in)           :: dict
      character(nameLen), intent(in)          :: name
      character(100),parameter :: Here = 'init (aniImplicitChaosClerk.f90)'

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


      self % startPop = ZERO
    end subroutine init

    !! Return to uninitialised state
    elemental subroutine kill(self)
      class(aniImplicitChaosClerk), intent(inout) :: self
      
      ! Kill superclass
      call kill_super(self)

      if (allocated(self % chaosOfPop)) deallocate(self % chaosOfPop)

    end subroutine kill


    ! Only needed at start/end of generation
    function validReports(self) result(validCodes)
      class(aniImplicitChaosClerk), intent(in)          :: self
      integer(shortInt), dimension(:), allocatable  :: validCodes

      validCodes = [hist_CODE, cycleStart_CODE]
    end function validReports
    
    !!
    !! Return memory size of the clerk
    !!
    !! See tallyClerk_inter for details
    !!

    elemental function getSize(self) result(S)
      class(aniImplicitChaosClerk), intent(in) :: self
      integer(shortInt)                    :: S
    
      S = (self % P + 1)**3
    end function getSize

    subroutine reportCycleStart(self, start, mem)
      class(aniImplicitChaosClerk), intent(inout) :: self
      class(particleDungeon), intent(in)      :: start
      type(scoreMemory), intent(inout)        :: mem

      self % startPop = start % popWeight()
    
    end subroutine reportCycleStart
  
  
    !!
    !! Process history report
    !! Gets fate code from the particle
    !!
    !! See tallyClerk_inter for details
    !!
    subroutine reportHist(self, p, xsData, mem)
      class(aniImplicitChaosClerk), intent(inout) :: self
      class(particle), intent(in)             :: p
      class(nuclearDatabase),intent(inout)    :: xsData
      type(scoreMemory), intent(inout)        :: mem
      real(defReal), dimension(3, self % P + 1)  :: legendrePol
      real(defReal)                           :: score
      integer(shortInt)                       :: i,j,k, nb_coeffs


      nb_coeffs = self % P + 1
  
      if (p % fate /= leak_FATE) then
        ! Evaluate Legendre polynomials up to right order 
        legendrePol(1,:) = evaluateLegendre(self % P, p % X(1)) 
        legendrePol(2,:) = evaluateLegendre(self % P, p % X(2)) 
        legendrePol(3,:) = evaluateLegendre(self % P, p % X(3)) 
        do k = 1, self % P + 1
          do j = 1, self % P + 1
            do i = 1, self % P + 1
              score = (2*(i-1) + 1) * (2*(j-1) + 1) * (2*(k-1) + 1) * legendrePol(1,i) * legendrePol(2,j) * legendrePol(3,k) &
                        * p % w / self % startPop
              call mem % score(score, self % getMemAddress() + i + (j-1) * nb_coeffs + (k-1) * (nb_coeffs**2) - 1)
            end do
          end do
        end do
      end if
    end subroutine reportHist

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
    subroutine display(self, mem)
        class(aniImplicitChaosClerk), intent(in)  :: self
        type(scoreMemory), intent(in)      :: mem
    
        print *, 'aniImplicitChaosClerk does not support display yet'
    
      end subroutine display
    
  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
      subroutine print(self, outFile, mem)
        class(aniImplicitChaosClerk), intent(in)       :: self
        class(outputFile), intent(inout)           :: outFile
        type(scoreMemory), intent(in)              :: mem
        real(defReal)                              :: val, std
        integer(shortInt)                          :: i, j, k, nb_coeffs
        integer(shortInt),dimension(:),allocatable :: resArrayShape
        character(nameLen)                         :: name
        
        nb_coeffs = self % P + 1
        ! Begin block
        call outFile % startBlock(self % getName())
        
        ! Write results.
        ! Get shape of result array
        resArrayShape = [self % getSize()]
    
        ! Start array
        name ='Res'
        call outFile % startArray(name, resArrayShape)
    
        ! Print results to the file
        do k = 1, nb_coeffs
          do j = 1, nb_coeffs
            do i = 1, nb_coeffs
              call mem % getResult(val, std, self % getMemAddress() + i + (j-1) * nb_coeffs + (k-1) * nb_coeffs**2 - 1)
              call outFile % addResult(val,std)
            end do
          end do
        end do
    
        call outFile % endArray()
        call outFile % endBlock()
    
      end subroutine print
    
    end module aniImplicitChaosClerk_class    



