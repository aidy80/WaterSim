module Constants 
    integer, parameter :: dp = selected_real_kind(15, 307)!Give reals double precision

    integer, parameter :: numAtomsPerMolecule = 3 !Number of atoms in each Oxygen 
    integer, parameter :: numCombos = 3!Number of atomic combos (OO, HH, OH)
    integer, parameter :: numMolecules = 901 !Number of molecules in this simulation
    integer, parameter :: numDimensions = 3 !Dimension of position 
                                            !and velocity space
    integer, parameter :: CvvStep = 500 !Number of items in Cvv
    integer, parameter :: MSDStep = 5000 !Number of items in MSD

    real(dp), parameter :: desiredOHBondLength = 1.0E-10 ![m]Bond length betweem 
                                                         !oxygen atoms
    real(dp), parameter :: desiredOHBondLengthSq = desiredOHBondLength**2![m^2]
    real(dp), parameter :: desiredHHBondLength = 1.631E-10 ![m]Bond length betweem 
                                                         !oxygen atoms
    real(dp), parameter :: desiredHHBondLengthSq = desiredHHBondLength**2![m^2]
    real(dp), parameter :: allowedError = 1.0E-8 ![m]Allowed error between actual 
                                                  !bond length and desired

    !Constants
    real(dp), parameter :: e = 2.71828 !Mathematical constant
    real(dp), parameter :: Bolz = 1.38064852E-23 ![J/K] Boltzmann constant
    real(dp), parameter :: avagadrosNum = 6.0221409E23![molecules/mole]
    real(dp), parameter :: pi = 3.14159265359
    real(dp), parameter :: eCharge = 1.6021766208E-19![C] Charge of an electron

    real(dp), parameter :: diElec = 1.0
    real(dp), parameter :: electricConst = 8.854187817E-12 ![F/m] Permittivity free space
    real(dp), parameter :: electricCoeff = 1 / (4.0 * pi * diElec * electricConst)
    real(dp), parameter :: hCharge = 0.41 * eCharge ![C]
    real(dp), parameter :: oCharge = -2.0 * hCharge ![C]
    real(dp), parameter :: hMass = (1.008 / avagadrosNum) / 1000.0 ![kg] Mass 
                                                                !of an hydrogen atom
    real(dp), parameter :: oMass = (15.999 / avagadrosNum) / 1000.0 ![kg] Mass 
                                                                !of an oxygen atom
    real(dp), parameter :: timeStep = 2.0E-15 ![s] Time between calculations
    real(dp), parameter :: delR = 0.003E-9 ![m]Distance between shells used for g(r)
    real(dp), parameter :: nonBondCut = 14.0E-10 ![m]Cutoff for nonBonded pairs
    real(dp), parameter :: nonBondCutSq = nonBondCut**2 ![m]Cutoff for nonBonded pairs
    real(dp), parameter :: ljA = 0.37122E-9 ![(kJ/mol)^1/6]Parameter for LJ equation
    real(dp), parameter :: ljB = 0.3428E-9 ![(kJ/mol)^1/12]Parameter for LJ equation
    real(dp), parameter :: sigma = 3.5533E-10 / (2**(1.0/6.0))
    real(dp), parameter :: sigmaSq = sigma**2
    real(dp), parameter :: epsilon = 0.1553 * 4.184 * 1000.0 / avagadrosNum ![kJ/mol]
    real(dp), parameter :: fourEps = 4.0 * epsilon
    real(dp), parameter :: twentyFourEps = 24.0 * epsilon

    !Speed Distribution
    integer, parameter :: numVelBox = 200
    integer, parameter :: numMaxwell = 800
    real(dp), parameter :: boxScale = 0.1

    real(dp), parameter :: cutoff = 12.0E-10![m] Cutoff distance for short-range interaction
    real(dp), parameter :: cutoffSq  = cutoff**2![m^2] The cutoff radius squared 
                                                !used for optimization

    integer, parameter :: degreesFreedom = 6
    real(dp), parameter :: temperature = 298.0 ![K] initial const temperature 
                                                    !of the system
    real(dp), parameter :: temperatureTol = 5.0 ![K] Amount temperature is able to
                                                    !to fluctuate

    integer, parameter :: numEquilSteps = 50000 !Number of timesteps for equil sim
    integer, parameter :: nonAnalSteps = 50000 !Number of timesteps without analysis
    integer, parameter :: analSteps = 300000 !Number of timesteps in the program

    integer, parameter :: zeroMomentTimeStep = 100 !Number of timesteps 
                                                   !between momentum-zeroing
    integer, parameter :: numTrajSteps = 100 !Number of timesteps between 
                                             !trajectory outputs to .ttr file
    integer, parameter :: temperatureStep = 1000 !Number of timesteps between 
                                                 !temperature decrements
    integer, parameter :: bondedUpdateStep = 20 !How often you update the non-Bonded list


    !Charges of water molecule
    real(dp), dimension(numAtomsPerMolecule), parameter :: molCharges = &
                                                      &(/oCharge, hCharge, hCharge/)
    !Masses of water molecule
    real(dp), dimension(numAtomsPerMolecule), parameter :: molMasses = &
                                                            &(/oMass, hMass, hMass/)
    !Inverse masses of water molecule
    real(dp), dimension(numAtomsPerMolecule), parameter :: invMolMasses = &
                                                &(/(1.0/oMass), (1.0/hMass), (1.0/hMass)/)
    real(dp), dimension(numAtomsPerMolecule, numAtomsPerMolecule),parameter :: bondLengthSq = &
                        &reshape((/0.0D0, desiredOHBondLengthSq, desiredOHBondLengthSq,&
                        &desiredOHBondLengthSq, 0.0D0, desiredHHBondLengthSq,&
                        &desiredOHBondLengthSq, desiredHHBondLengthSq, 0.0D0/), (/3,3/))

end module Constants
