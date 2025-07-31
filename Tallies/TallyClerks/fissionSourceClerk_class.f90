module fissionSourceClerk_class

  use numPrecision
  use tallyCodes
  use dictionary_class,      only : dictionary
  use genericProcedures,     only : fatalError
  use particle_class,        only : particle, particleState
  use particleDungeon_class, only : particleDungeon
  use outputFile_class,      only : outputFile

  use scoreMemory_class,     only : scoreMemory
  use tallyResult_class,     only : tallyResult, tallyResultEmpty
  use tallyClerk_inter,      only : tallyClerk, kill_super => kill

  use tallyFilter_inter,     only : tallyFilter
  use tallyFilterFactory_func, only : new_tallyFilter

  use tallyMap_inter,        only : tallyMap
  use tallyMapFactory_func,  only : new_tallyMap

  implicit none
  private

  !!
  !! Simplest possible analog estimator of the fission source distribution
  !! 
  !!
  !! Private Members:
  !! 
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myClerk {
  !!   type fissionSourceClerk;
  !!   map { type any map }
  !! }
  !!
  type, public,extends(tallyClerk) :: fissionSourceClerk
    private
    class(tallyFilter), allocatable                       :: filter
    class(tallyMap), allocatable                          :: map

    integer(shortInt)   :: width
    
  contains
    ! Procedures used during build
    procedure :: init
    procedure :: kill
    procedure :: validReports
    procedure :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: reportCycleEnd

    ! Output procedures
    procedure  :: display
    procedure  :: print
  end type fissionSourceClerk

contains

  !!
  !! Initialise from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(fissionSourceClerk), intent(inout) :: self
    class(dictionary), intent(in)        :: dict
    character(nameLen), intent(in)       :: name
    character(100),parameter :: Here = 'init (fissionSourceClerk.f90)'

    ! Assign name
    call self % setName(name)

    ! Load filetr
    if( dict % isPresent('filter')) then
      call new_tallyFilter(self % filter, dict % getDictPtr('filter'))
    end if

    ! Load map
    if( dict % isPresent('map')) then
      call new_tallyMap(self % map, dict % getDictPtr('map'))
    end if

    ! Set width
    self % width = 1

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(fissionSourceClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill and deallocate filter
    if (allocated(self % filter)) then
      deallocate(self % filter)
    end if

    ! Kill and deallocate map
    if (allocated(self % map)) then
      call self % map % kill()
      deallocate(self % map)
    end if

    self % width   = 0

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(fissionSourceClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [ cycleEnd_CODE ]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(fissionSourceClerk), intent(in) :: self
    integer(shortInt)                  :: S

    S = self % width
    if(allocated(self % map)) S = S * self % map % bins(0)
  end function getSize

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(fissionSourceClerk), intent(inout) :: self
    class(particleDungeon), intent(in)    :: end
    type(scoreMemory), intent(inout)      :: mem
    real(defReal)                         :: endPop, scoreVal
    integer(shortInt)                     :: i
    type(particle)                        :: p
    type(particleState)                   :: state
    integer(longInt)                      :: addr, binIdx

    ! Update end population weight
    endPop = end % popWeight()

    ! Close batch
    if( mem % lastCycle() ) then
      do i = 1, end % popSize()
        call end % copy(p, i)
        state = p
        state % X = state % Xold

          ! Check if within filter
        if (allocated(self % filter)) then
          if (self % filter % isFail(state)) return
        end if

        ! Find bin index
        if (allocated(self % map)) then
          binIdx = self % map % map(state)
        else
          binIdx = 1
        end if

        ! Return if invalid bin index
        if (binIdx == 0) return
        
        ! Calculate bin address
        addr = self % getMemAddress() + self % width * (binIdx - 1)  - 1
        
        ! Append all bins
        scoreVal = state % wgt / endPop
        call mem % score(scoreVal, addr + 1)
      end do
    end if

  end subroutine reportCycleEnd

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(fissionSourceClerk), intent(in)  :: self
    type(scoreMemory), intent(in)       :: mem

    print *, 'collisionClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk in the slot to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(fissionSourceClerk), intent(in)          :: self
    class(outputFile), intent(inout)           :: outFile
    type(scoreMemory), intent(in)              :: mem
    real(defReal)                              :: val, std
    integer(shortInt)                          :: i
    integer(shortInt),dimension(:),allocatable :: resArrayShape
    character(nameLen)                         :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If collision clerk has map print map information
    if (allocated(self % map)) then
      call self % map % print(outFile)
    end if

    ! Write results.
    ! Get shape of result array
    if (allocated(self % map)) then
      resArrayShape = [self % width, self % map % binArrayShape()]
    else
      resArrayShape = [self % width]
    end if

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

end module fissionSourceClerk_class
