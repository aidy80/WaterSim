!Subrouintes for md functions, force calculations, leapfrogs, initializations, etc

module mdSubroutines
    use Constants
    implicit none

contains
    !Helper subroutine to initialize a 2D array of dimension 
    !(arraySize, numDimensions) with all zeros. "array" is the 2D array passed 
    !and returned with all zeros
    subroutine initZero(array)
        implicit none
        
        real(dp), dimension(numMolecules,numAtomsPerMolecule, numDimensions), &
                                                        &intent(inout) :: array

        integer :: m,l,k

        do m = 1, numMolecules 
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    array(m,l,k) = 0.0
                end do
            end do
        end do
    end subroutine initZero

    subroutine openNVEFiles()
        integer :: numAtoms

        character :: nullChar

        !Open file containing position and velocity information 
        !about argon atoms in the cubic boundary from argon.gro
        open(unit=11, file='Equil_final.gro') 
        open(unit=90, file='NVE_temperature.dat')
        open(unit=91, file='NVE_totEnergy.dat')
        open(unit=92, file='NVE_potentialEnergy.dat')
        open(unit=93, file='NVE_kineticEnergy.dat')
        open(unit=94, file='NVE.gro')
        open(unit=95, file='NVE_final.gro')
        open(unit=96, file='NVE_TCFNorm.dat')
        open(unit=97, file='NVE_grOO.dat')
        open(unit=98, file='NVE_MSD.dat')
        open(unit=99, file='NVE_grHH.dat')
        open(unit=100, file='NVE_TCF.dat')
        open(unit=101, file='NVE_grOH.dat')
        open(unit=111, file='NVE_speedDist.dat') 
        open(unit=112, file='NVE_maxwellDist.dat')
        
        !Read in header information from the file
        read(11, *) nullChar
        read(11, 30) numAtoms
        if (numAtoms / numAtomsPerMolecule.NE.numMolecules) then
            print *, "Incorrect number of molecules read in"
        end if

        30 format(I5)
    end subroutine openNVEFiles

    subroutine openNVTFiles()
        character :: nullChar

        integer :: numAtoms

        !Open file containing position and velocity information 
        !about argon atoms in the cubic boundary from argon.gro
        open(unit=11, file='Equil_final.gro') 
        open(unit=90, file='NVT_temperature.dat')
        open(unit=91, file='NVT_totEnergy.dat')
        open(unit=92, file='NVT_potentialEnergy.dat')
        open(unit=93, file='NVT_kineticEnergy.dat')
        open(unit=94, file='NVT.gro')
        open(unit=95, file='NVT_final.gro')
        open(unit=96, file='NVT_TCFNorm.dat')
        open(unit=97, file='NVT_grOO.dat')
        open(unit=98, file='NVT_MSD.dat')
        open(unit=99, file='NVT_grHH.dat')
        open(unit=100, file='NVT_TCF.dat')
        open(unit=101, file='NVT_grOH.dat')
        open(unit=111, file='NVT_speedDist.dat') 
        open(unit=112, file='NVT_maxwellDist.dat')
        
        !Read in header information from the file
        read(11, *) nullChar
        print *, nullChar
        read(11, 30) numAtoms
        print *, numAtoms
        if (numAtoms / numAtomsPerMolecule.NE.numMolecules) then
            print *, "Incorrect number of molecules read in"
        end if

        30 format(I5)
    end subroutine openNVTFiles
    
    subroutine openEquilFiles()
        character :: nullChar

        integer :: numAtoms

        !Open file containing position and velocity information 
        !about argon atoms in the cubic boundary from argon.gro
        open(unit=11, file='initWater.gro') 
        open(unit=90, file='Equil_temperature.dat')
        open(unit=91, file='Equil_totEnergy.dat')
        open(unit=92, file='Equil_potentialEnergy.dat')
        open(unit=93, file='Equil_kineticEnergy.dat')
        open(unit=94, file='equil.gro')
        open(unit=95, file='Equil_final.gro')
    
        !Read in header information from the file
        read(11, *) nullChar
        print *, nullChar
        read(11, 30) numAtoms
        print *, numAtoms
        if (numAtoms / numAtomsPerMolecule.NE.numMolecules) then
            print *, "Incorrect number of molecules read in"
        end if

        30 format(I5)
    end subroutine openEquilFiles

    !Read in position, velocity and dimension information from the initial file
    subroutine readIn(pos, vel, dim)
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: pos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: vel
        real(dp), dimension(numDimensions), intent(inout) :: dim

        !Garbage variables
        character :: nullChar
        character :: nullChar1
        integer :: nullInt
        integer :: nullInt1

        !Iterators
        integer :: m,l,k

        !Read in position and velocity information from the .gro file
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                read(11, 10) nullInt, nullChar, nullChar1, nullInt1, pos(m,l,1), &
                        & pos(m,l,2), pos(m,l,3), vel(m,l,1), vel(m,l,2), vel(m,l,3)
            end do
        end do

        !Read in the dimensions of the box and convert to [m]
        read(11,*) dim(1),dim(2),dim(3)
        do m = 1, numDimensions
            dim(m) = dim(m) * 1.0E-9
        end do

        close (unit=11)

        !Convert read-in pos/velocity information from nm and nm/ps to m and m/s
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    pos(m,l,k) = pos(m,l,k) * 1.0E-9
                    vel(m,l,k) = vel(m,l,k) * 1.0E3
                end do
            end do
        end do

        30 format(I5)
        20 format(F10.5, F10.5, F10.5)
        10 format(i5,2a5,i5,3f8.3,3f8.3)
    end subroutine readIn

    !Zeros the net magnitude of the vectors passed in. Used to prevent the 
    !Wandering ice cube problem
    subroutine zeroNetMomentum(vel)
        implicit none
        
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                          &intent(inout) :: vel

        integer :: m, l, k

        real(dp), dimension(numDimensions) :: netVel
        real(dp), dimension(numDimensions) :: velScale

        !Find the current velocity of the system
        do k = 1, numDimensions
            netVel(k) = 0.0
        end do

        !Store the netVelcoity in each direction of the system
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    netVel(k) = netVel(k) + vel(m,l,k) 
                end do
            end do
        end do

        !Calc how much to scale each velocity by such that net velocity is zero
        do k = 1, numDimensions
            velScale(k) = netVel(k) / real(numMolecules * numAtomsPerMolecule)
        end do

        !Scale the velocity of each atom in the system
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    vel(m, l, k) = vel(m, l, k) - velScale(k)
                end do
            end do
        end do
    end subroutine zeroNetMomentum

    !Calculate the forces between each pair of atoms and record the net force
    !on each atom
    subroutine calcForces(p, potentialEnergy, pos, force, bins, numBins, dim, &
                                            &ensemble, bondedList, bondedLength)
        integer, intent(in) :: p
        real(dp), intent(inout) :: potentialEnergy
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                               &intent(in) :: pos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                          &intent(inout) :: force
        integer, dimension(numMolecules, numMolecules) :: bondedList
        integer, dimension(numMolecules) :: bondedLength
        integer, intent(in) :: numBins
        integer, dimension(numCombos, numBins), intent(inout) :: bins
        real(dp), dimension(numDimensions), intent(in) :: dim
        character, intent(in) :: ensemble

        real(dp) :: potentialSum 
        real(dp) :: distanceSq
        real(dp) :: distanceMag
        real(dp) :: currPotential
        real(dp) :: Fmag
        real(dp), dimension(numDimensions) :: distance

        integer :: pairIndex

        real(dp) :: sigmaDistTwo 
        real(dp) :: sigmaDistSix 
        real(dp) :: sigmaDistTwelve 

        integer :: currIndex
        integer :: i,j,m,l,k

        !Compute Force between each pair of atoms if their distance is below
        !the cutoff radius. Each pair is evaluated only once
        do i = 1, numMolecules 
            potentialSum = 0.0
            do m = 1, numAtomsPerMolecule
                do j = 1, bondedLength(i)
                    pairIndex = bondedList(i, j)
                    do l = 1, numAtomsPerMolecule
                    !Calculate the distance between the current atoms, applying 
                    !periodic boundary conditions to find the minimum possible 
                    !distance for each coordinate
                        do k = 1, numDimensions
                            distance(k) = pos(i, m, k) - pos(pairIndex, l, k)
                            distance(k) = distance(k) - dim(k) * &
                                                  &anint(distance(k) / dim(k))
                        end do

                        distanceSq = distance(1)**2 + distance(2)**2 + distance(3)**2

                        !If the distance between the two atoms is below the cutoff, 
                        !calculate the force exerted on each of them based on the 
                        !lennard jones potential
                        if (distanceSq < cutoffSq) then
                            distanceMag = sqrt(distanceSq)
                            Fmag = 0.0

                            !Calc LJ between oxygen atoms
                            if (l==1.AND.m==1) then
                                sigmaDistTwo = sigmaSq / distanceSq                
                                sigmaDistSix = sigmaDistTwo**3                     
                                sigmaDistTwelve = sigmaDistSix**2                  
                                                                                  
                                !Calc potential from lennard jones equation between the 
                                !current pair of atoms and add it to the current   
                                !potential energy sum                              
                                potentialSum = potentialSum + forceShift(distanceMag)&
                                            &* fourEps * (sigmaDistTwelve - sigmaDistSix) 
                                                                                  
                                !Calculate the resulting force on the current two atoms 
                                !based on the lennard jones potential between them. 
                                !Calculated using the negative gradient                       
                                Fmag = twentyFourEps * (2.0 * sigmaDistTwelve - &
                                                    &sigmaDistSix) / distanceMag
                            end if

                            !Calculate the electrostatic charges 
                            Fmag = Fmag + electricCoeff * molCharges(l) * molCharges(m) / &
                                         &distanceSq - electricCoeff * molCharges(l) * &
                                         &molCharges(m) / cutoffSq
                            !print *, "electric: ", Fmag
                            potentialSum = potentialSum + (forceShift(distanceMag) * &
                                        &electricCoeff * molCharges(l) * molCharges(m) / &
                                        &distanceMag) 

                            !Apply the electrostatic force to the net force of each atom
                            do k = 1, numDimensions
                                force(i, m, k) = force(i, m, k) + &
                                                    &Fmag * (distance(k) / distanceMag)
                                force(pairIndex, l, k) = force(pairIndex, l, k) - &
                                                    &Fmag * (distance(k) / distanceMag)
                            end do
                        end if

                    end do
                end do
            end do
            potentialEnergy = potentialEnergy + potentialSum
        end do
    end subroutine calcForces

    real function forceShift(radius)
        real(dp), intent(in) :: radius

        forceShift = (1 - (radius / cutoff))**2
        return 
    end function 

    !Update position and velocity using leapfrog algorithm and shake algorithm
    subroutine leapFrogAndShake(kineticEnergy, force, vel, pos, oldPos, oldVel)
        real(dp), intent(inout) :: kineticEnergy

        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: vel
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: pos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: oldPos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: oldVel
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: force

        real(dp) :: addedKE
        integer :: m, l, k

        !Use the leap-frog verlet algorithm to calculate new position and 
        !velocity vectors for all atoms based on the forces 
        !calculated between them.
        do m = 1, numMolecules
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions 
                    oldVel(m, l, k) = vel(m, l, k)
                    oldPos(m, l, k) = pos(m, l, k)
                    vel(m, l, k) = vel(m, l, k) + (force(m, l, k) / &
                                                        &molMasses(l)) * timeStep
                    pos(m, l, k) = pos(m, l, k) + vel(m, l, k) * timeStep
                end do
            end do

            call Shake(m, oldPos, pos)

            !Update velocity baseed on changed velocity and calc KE 
            !post-bondlength adjustments
            addedKE = 0.0
            do l = 1, numAtomsPerMolecule
                do k = 1, numDimensions
                    vel(m,l,k) = (pos(m,l,k) - oldPos(m,l,k)) / timeStep
                    addedKE = addedKE + 0.5 * (0.5*(oldVel(m,l,k) + vel(m,l,k)))**2 * molMasses(l)
                end do
            end do

            kineticEnergy = kineticEnergy + addedKE
        end do
    end subroutine leapFrogAndShake 

    !SHAKE Alg to fix bond angles
    subroutine Shake(m, oldPos, pos)
        integer, intent(in) :: m
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(in) :: oldPos
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout) :: pos

        logical :: hasConverge !Used in iterative while loop
        integer :: convergeCount

        real(dp), dimension(numDimensions) :: pastBondLength
        real(dp), dimension(numDimensions) :: currBondLength
        real(dp) :: currError
        real(dp) :: lambda
        real(dp) :: currBondLengthSq

        integer :: i,j,k
        
        hasConverge = .false.
        convergeCount = 0

        !Use iteration to solve for lambda, a coefficient adjusting the amount
        !that positions/velocities must be adjusted to prevent bond lengths from 
        !changing
        do while(.NOT.hasConverge)
            do i = 1, numAtomsPerMolecule - 1
                do j = i + 1, numAtomsPerMolecule
                    do k = 1, numDimensions
                        pastBondLength(k) = oldPos(m, i, k) - oldPos(m, j, k)
                        currBondLength(k) = pos(m, i, k) - pos(m, j, k)
                    end do
           
                    lambda = (dot(currBondLength, currBondLength) - bondLengthSq(i, j)) / &
                        &(2.0*(invMolMasses(i) + invMolMasses(j))*dot(currBondLength, pastBondLength))
            
                    !Use lambda to adjust the position and velocities of atoms in the simulation
                    do k = 1, numDimensions
                        pos(m, i, k) = pos(m, i, k) - pastBondLength(k) * invMolMasses(i) * lambda
                        pos(m, j, k) = pos(m, j, k) + pastBondLength(k) * invMolMasses(j) * lambda
                    end do
                end do
            end do

            do i = 1, numAtomsPerMolecule - 1
                do j = i + 1, numAtomsPerMolecule
                    currBondLengthSq = 0.0
                    do k = 1, numAtomsPerMolecule
                        currBondLengthSq = currBondLengthSq + (pos(m,i,k) - pos(m,j,k))**2
                    end do

                    !Find the error between the current bond 
                    !length and desired bond length
                    currError = ABS(currBondLengthSq - bondLengthSq(i, j)) / (bondLengthSq(i, j))

                    !Repeat the iteration until the bond is within 1E-18 of the desired 
                    !length or the iteration occurs 500 times
                    if(currError.LE.allowedError) then
                        hasConverge = .true.
                    else 
                        hasConverge = .false.
                        convergeCount = convergeCount + 1

                        if (convergeCount > 500) then
                            print *, "Did not converge"
                            hasConverge = .true.
                        end if
                    end if
                end do
            end do
        end do

        if (m==1) then
            call testBondLengths(pos)
        end if
    end subroutine Shake

    !Scale the velocities in the system such that the system has a 
    !temperature equal to desiredTemperature
    subroutine scaleTemp(kineticEnergy, vel)
        implicit none 
    
        real(dp), intent(in) :: kineticEnergy
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                            &intent(inout)::vel

        real(dp) :: currTemp 
        real(dp) :: tempScale

        integer :: m,l,k

        !Use equipartition thm to find the current temperature of the system
        currTemp = (kineticEnergy * 2.0D0) / &
                        &(real(degreesFreedom) * real(numMolecules) * Bolz)
        print *, "Scaling temp", currTemp
        tempScale = sqrt(temperature / currTemp)


        !Scale each velocity such that the kintic energy of the system can 
        !be applied to the equipartition theorem to find that the system has 
        !a temperature of "desiredTemperature"
        do m=1, numMolecules
            do l=1, numAtomsPerMolecule
                do k=1, numDimensions
                    vel(m,l,k) = vel(m,l,k) * tempScale
                end do
            end do
        end do
    end subroutine scaleTemp

    !Update the nonBonded pairs list
    subroutine findNonBondedPairs(p, bondedList, bondedLength, pos, dim, &
                                                    &ensemble, bins, numBins)
        integer, intent(in) :: p
        character, intent(in) :: ensemble
        integer, dimension(numMolecules, numMolecules), intent(inout) :: bondedList
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                                                                &intent(in) :: pos
        real(dp), dimension(numDimensions), intent(in) :: dim

        integer, intent(in) :: numBins
        integer, dimension(numCombos, numBins), intent(inout) :: bins

        integer, dimension(numMolecules) :: bondedLength
        real(dp), dimension(numDimensions) :: distance
        real(dp) :: distanceSq

        integer :: nextIndexi
        integer :: nextIndexj
        integer :: currIndex

        integer :: m,l,k
        integer :: i,j,q,r

        do m = 1, numMolecules
            bondedLength(m) = 0
        end do

        do i = 1, numMolecules - 1
            do q = 1, numAtomsPerMolecule
                do j = i + 1, numMolecules 
                    do r = 1, numAtomsPerMolecule
                        do k = 1, numDimensions
                            distance(k) = pos(i, q, k) - pos(j, r, k)
                            distance(k) = distance(k) - dim(k) * &
                                                  &anint(distance(k) / dim(k))
                        end do

                        distanceSq = distance(1)**2 + distance(2)**2 + distance(3)**2

                        if (q==1.AND.r==1.AND.distanceSq.LE.nonBondCutSq) then
                            bondedList(i, bondedLength(i) + 1) = j
                            bondedLength(i) = bondedLength(i) + 1
                        end if

                        if (p.GE.nonAnalSteps.AND.ensemble=="Y") then
                            currIndex = Int(sqrt(distanceSq)/delR) + 1
                            if (q==1.AND.r==1) then
                                bins(1, currIndex) = bins(1, currIndex) + 2
                            else if(q==1.OR.r==1) then
                                bins(3, currIndex) = bins(3, currIndex) + 2
                            else 
                                bins(2, currIndex) = bins(2, currIndex) + 2
                            end if
                        end if
                    end do
                end do
            end do
        end do
    end subroutine findNonBondedPairs

    subroutine measureGr(pos, bins, numBins, dim)
        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions) :: pos
        integer, intent(in) :: numBins
        integer, dimension(numCombos, numBins), intent(inout) :: bins
        real(dp), dimension(numDimensions), intent(in) :: dim

        integer :: currIndex
        real(dp), dimension(numDimensions) :: distance
        real(dp) :: distanceSq

        integer :: i, j, q, r, k

        do i = 1, numMolecules - 1
            do q = 1, numAtomsPerMolecule
                do j = i + 1, numMolecules
                    do r = 1, numAtomsPerMolecule
                        do k = 1, numDimensions
                            distance(k) = pos(i, q, k) - pos(j, r, k)
                            distance(k) = distance(k) - dim(k) * &
                                                  &anint(distance(k) / dim(k))
                        end do

                        distanceSq = distance(1)**2 + distance(2)**2 + distance(3)**2

                        currIndex = Int(sqrt(distanceSq)/delR) + 1
                        if (q==1.AND.r==1) then
                            bins(1, currIndex) = bins(1, currIndex) + 2
                        else if(q==1.OR.r==1) then
                            bins(3, currIndex) = bins(3, currIndex) + 2
                        else 
                            bins(2, currIndex) = bins(2, currIndex) + 2
                        end if
                    end do
                end do
            end do
        end do

    end subroutine measureGr

    !See if the system needs its temperature scaled
    logical function needsScaling(kineticEnergySum)
        real(dp) :: kineticEnergySum

        real(dp) :: avgKE
        real(dp) :: currTemp

        avgKE = kineticEnergySum / real(temperatureStep)
        currTemp = (avgKE * 2.0D0) / &
                        &(real(degreesFreedom) * real(numMolecules) * Bolz)
        print *, "Curr temp pre sacling", currTemp
        needsScaling = .false.

        if(ABS(currTemp - temperature).GE.temperatureTol) then
            needsScaling = .true. 
        end if

        return 
    end function

    !Close all NVE files
    subroutine closeAnalFiles()
        close (unit=90)
        close (unit=91)
        close (unit=92)
        close (unit=93)
        close (unit=94)
        close (unit=95)
        close (unit=96)
        close (unit=97)
        close (unit=98)
        close (unit=99)
        close (unit=100)
        close (unit=101)
        close (unit=111)
        close (unit=112)
    end subroutine closeAnalFiles

    !Close all files
    subroutine closeEquilFiles()
        close (unit=90)
        close (unit=91)
        close (unit=92)
        close (unit=93)
        close (unit=94)
        close (unit=95)
    end subroutine closeEquilFiles

    !Test that the Shake Algorithm successfully fixes bond length
    subroutine testBondLengths(pos)
        implicit none

        real(dp), dimension(numMolecules, numAtomsPerMolecule, numDimensions), &
                        &intent(in) :: pos

        integer :: m,l,k,i,j

        real(dp) :: currBondLengthSq


        do m = 1, numMolecules
            do i = 1, numAtomsPerMolecule
                do j = 1, numAtomsPerMolecule
                    currBondLengthSq = (pos(m, i, 1) - pos(m, j, 1))**2 +&
                                    & (pos(m, i, 2) - pos(m, j, 2))**2 +&
                                    & (pos(m, i, 3) - pos(m, j, 3))**2
                    if (ABS(sqrt(bondLengthSq(i,j)) - sqrt(currBondLengthSq)).GE.1.0E-12) then
                        print *, "desiredBondLength: ", sqrt(bondLengthSq(i,j)), &
                                &"currBondLength: ", sqrt(currBondLengthSq)
                    end if
                end do
            end do
        end do
    end subroutine testBondLengths

    !Compute the dot product of two vectors v1 and v2
    real function dot(v1, v2)
        implicit none

        real(dp), dimension(numDimensions),intent(in) :: v1
        real(dp), dimension(numDimensions),intent(in) :: v2

        dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)
        return
    end function

end module mdSubroutines
