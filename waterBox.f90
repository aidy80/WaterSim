!Written by Aidan Fike
!February 24, 2018
!
!This program outputs a gromacs-formatted-file filled with 500 evenly-space 
!O2 molecules at 320 Kelvin. These atoms exist in a cube  of dimensions 
!22.45x22.45x22.45 nanometers and are in random orientations

program waterBox 
    implicit none

    !Initial vrms of atoms in the simulation nm/ps
    real :: vrms = 0.908387

    !Number of atoms being simulated
    integer, parameter :: numMolecules = 901
    integer :: currNumMolecules = numMolecules
    integer, parameter :: atomsPerMolecule = 3
    integer, parameter :: numDimensions = 3
    integer :: atomicCount = 1
    real, parameter :: pi  = 3.1415926

    !Number of atoms per side of the cube
    integer :: xNum = 10
    integer :: yNum = 10
    integer :: zNum = 10

    !Dimension of the cube in nanometers
    real :: dim = 3.000

    !The distance between oxygen atoms
    real :: ohBondLength = 0.1000000
    real :: hhBondLength = 0.1631064
    real, parameter :: bondAngle = 109.28
    real :: rad = bondAngle * pi / 180.0
    real, dimension(3) :: orientation
    real, dimension(3) :: orientationAxis

    !Used to hold temporary variables for velocity values
    real :: a
    real :: b
    real :: c
    real :: velMag !Used to normalize the randomly-created normal vector
    real :: orientMag !Used to normalize the randomly-created normal vector
    real :: orientAxisMag !Used to normalize the randomly-created normal vector
    real :: oldOrientOne !Stores orient1

    !Vectors to hold velocity and position information about a given vector
    real, dimension(3) :: pos, vel

    !Counters to look through all possible atoms
    integer :: x, y, z, i, j, k

    !Relation between curr iterator (x,y,z) to dimension
    real :: dmdx
    real :: dmdy 
    real :: dmdz 

    !Create a random seed for random number generation. 
    !This seed will be modifed upon each iteration using scale, increment, max
    integer :: seed = 1010

    !Initialize the use of rand()
    call random_seed(size = seed)
    call srand(seed)

    !Calc distance between molecules
    dmdx = dim / real(xNum)
    dmdy = dim / real(yNum)
    dmdz = dim / real(zNum)

    !Open a file for writing argon atom information formatted for gromacs
    open(unit = 11, file='initWater.gro')

    !The file's header information
    write(11, *) "Initial box of 901 water molecules, t = 0.0ns"
    write(11, 30) numMolecules * atomsPerMolecule

    !Treating the box as a 3D grid, loop through it and create equally 
    !spaced atoms throughout.
    do x = 0, xNum - 1
        do y = 0, yNum - 1
            do z = 0, zNum - 1
                !Because this loop will go through more than numAtoms 
                !number of atoms, stop outputting atom information when you 
                !have outputted numAtoms of information
                if (currNumMolecules == 0) then
                    goto 100
                end if

                !Create a velocity vector in a random direction with 
                !magnitude vrms (nm/ps)
                do i = 1, 3
                    vel(i) = 2 * rand() - 1
                    orientationAxis(i) = rand()
                end do

                !orientation(1) = (rand() * 2) - 1
                !orientation(2) = (rand() * 2) - 1
                orientation(1) = 1.0
                orientation(2) = 0.0
                orientation(3) = 0.0

                !Normalize vectors
                orientMag = sqrt(orientation(1)**2 + orientation(2)**2 + orientation(3)**2)
                orientAxisMag = sqrt(orientationAxis(1)**2 + orientationAxis(2)**2 + orientationAxis(3)**2)
                velMag = sqrt(vel(1)**2 + vel(2)**2 + vel(3)**2)
                do i = 1,3
                    orientation(i) = (real(orientation(i)) / orientMag) * ohBondLength 
                    orientationAxis(i) = orientationAxis(i) / orientAxisMag
                    vel(1) = (real(vel(1)) / velMag) * vrms
                end do

                !Create a position vector depending on 
                !the current values of x,y, and z 
                !(the current coordinated of the atom in the grid)
                pos(1) = real(x) * dmdx
                pos(2) = real(y) * dmdy
                pos(3) = real(z) * dmdz

                !Output the position and velocity information for the 
                !most-recently-created atom
                write(11, 10) abs(numMolecules + 1 - currNumMolecules), "WATER", "O1", &
                             &atomicCount, pos(1), pos(2), pos(3), &
                             &vel(1), vel(2), vel(3)
                atomicCount = atomicCount + 1

                write(11, 10) abs(numMolecules + 1 - currNumMolecules), "WATER", "H1", &
                             &atomicCount, pos(1) + orientation(1), &
                             &pos(2) + orientation(2), pos(3) + orientation(3), &
                             &vel(1), vel(2), vel(3)
                atomicCount = atomicCount + 1

                oldOrientOne = orientation(1)
                orientation(1) = orientation(1) * cos(rad) - orientation(2) * sin(rad)
                orientation(2) = oldOrientOne * sin(rad) + orientation(2) * cos(rad)

                write(11, 10) abs(numMolecules + 1 - currNumMolecules), "WATER", "H2", &
                             &atomicCount, pos(1) + orientation(1), &
                             &pos(2) + orientation(2), pos(3) + orientation(3), &
                             &vel(1), vel(2), vel(3)
                atomicCount = atomicCount + 1

                !Decrement the number of atoms left to output
                currNumMolecules = currNumMolecules - 1 
            end do
        end do
    end do

    100 continue

    !Output the box's dimensions
    write(11, 20) dim, dim, dim

    !Format information for the outputted file
    30 format(I5)
    20 format(F10.5, F10.5, F10.5)
    10 format(i5,2a5,i5,3f8.3,3f8.3)

    close(unit = 11)

end program waterBox
