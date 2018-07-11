!Aidan Fike
!March 2, 2017
!Program to simulate diatomic oxygen molecules in a cube based on lennard jones 
!iteractions. Program will simulate molecules for 100,000 4fs timesteps in a NVT
!ensemble. Temperature will begin at 320K and be gradually scaled to 77K

program equilSim
    use Constants
    use mdSubroutines
    use analysisSubroutines
    use testMod

    implicit none

    !Iterators
    integer :: p,l,m,k

    !Position and Velocity information about atoms in the system 
    real(dp), dimension(:, :, :), allocatable :: pos, vel !pos: [m], vel: [m/s]
    real(dp), dimension(:, :, :), allocatable :: oldPos, oldVel !oldPos: [m]
                                                                !oldVel: [m/s]
    real(dp), dimension(numDimensions) :: avgPos !Temp variable to measure 
                                                 !pos of an O2 molecule
    integer, dimension(numMolecules, numMolecules) :: bondedList
    integer, dimension(numMolecules) :: bondedLength
    

    !Force exerted on atoms at a given timestep
    real(dp), dimension(:, :, :), allocatable :: force ![N] 
    real(dp), dimension(numDimensions) :: dim![m] Size of each wall 
                                             !of the cubic enclosure

    real(dp) :: potentialEnergy ![J] Potential of the entire system
    real(dp) :: totEnergy ![J] Total energy of the system
    real(dp) :: potential ![J] the total potential energy of the system
    real(dp) :: kineticEnergy ![J] The total kinetic energy of the system
    real(dp) :: kineticEnergySum![J] Take a sum of KE of all steps

    real(dp), dimension(CvvStep) :: Cvv
    real(dp), dimension(MSDStep) :: MSD
    real(dp), dimension(:, :, :), allocatable :: vStore
    real(dp), dimension(:, :, :), allocatable :: pStore
    integer :: nCorr = 0!Used to count the number of times Cvv is added to 
    integer :: nCorr_MSD = 0!Used to count the number of times MSD is added to 
    integer :: nCorr_Rot = 0!Used to count the number of times CvvRot is added to 

    !Set up data for finding g(r)
    integer :: numBins 
    integer :: numFullBins 

    !Will have dimension (numCombos, numBins) where (1, *) is for OO, (2, *) is
    !for HH, (3, *) is for OH
    integer, dimension(:, :), allocatable :: bins

    !Speed Distribution
    integer, dimension(numVelBox) :: velBoxes

    !Used to time simulation
    real :: start_time
    real :: end_time

    real(dp) :: desiredTemperature = temperature

    CALL cpu_time(start_time)

    !Allocate space on the heap for position, velocity, and force information
    allocate(pos(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(vel(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(oldPos(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(oldVel(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(force(numMolecules, numAtomsPerMolecule, numDimensions))
    allocate(vStore(numMolecules, numDimensions, CvvStep))
    allocate(pStore(numMolecules, numDimensions, MSDStep))

    call openEquilFiles()
    call readIn(pos, vel, dim)

    call testLJ(sigma / (2.0**(1.0/6.0)), 0.1553D0 * 4.184D0)

    !Adjust the velocities in the system such that the net velocity 
    !in each direction is zero. This prevents wandering ice-cube problem
    call zeroNetMomentum(vel)

    kineticEnergySum = 0.0

    do p = 1, numEquilSteps
        if(mod(p,100) == 0) then
            print *, p
        end if
        !Init force vectors to zero
        call initZero(force)

        !Update the bonded pairs list
        if (mod(p, bondedUpdateStep) == 1) then
            call findNonBondedPairs(p, bondedList, bondedLength, pos, dim,&
                                                          & "N", bins, numBins)
        end if

        !Adjust the velocities in the system such that the net velocity 
        !in each direction is zero. This prevents wandering ice-cube problem
        if (mod(p, zeroMomentTimeStep) == 0) then
            call zeroNetMomentum(vel)
        end if

        potentialEnergy = 0.0
        kineticEnergy = 0.0

        !Calculate forces at the current timestep
        call calcForces(p, potentialEnergy, pos, force, bins, numBins, dim, &
                                                &"N", bondedList, bondedLength)

        !Update positions based on forces
        call leapFrogAndShake(kineticEnergy, force, vel, pos, oldPos, oldVel)
        
        !Find the kinetic energy of the system and 
        !the total energy of the system
        totEnergy = potentialEnergy + kineticEnergy
        kineticEnergySum = kineticEnergySum + kineticEnergy

        if (mod(p, temperatureStep/10) == 1) then
            !Scale the temperature of the system down to the desired value
            call scaleTemp(kineticEnergy, vel)
            kineticEnergySum = 0.0
        end if

        !Write vel and position information every 100 timesteps
        !Output the energy at the current timestep
        call writeEnergy(p, kineticEnergy, potentialEnergy, totEnergy)
        call writeEquilTrajectory(p, dim(1), pos, vel)
    end do

    !Free all heap memory
    deallocate(pos)
    deallocate(vel)
    deallocate(oldPos)
    deallocate(oldVel)
    deallocate(force)

    call closeEquilFiles()

    CALL cpu_time(end_time)
    print *, "Time usage:", (end_time - start_time)/60.0, " minutes", &
                            &mod((end_time - start_time),60.0), " seconds"

end program equilSim
