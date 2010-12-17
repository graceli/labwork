      PROGRAM avgE 

      IMPLICIT NONE

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 !                                                                                 !
 !                                                                                 !
 ! Reference: Park and Pande (2007) "Choosing weights for simulated tempering"     !
 !            Physical Review E. 76, 016703.                                       !
 !                                                                                 !
 ! Written: April 7, 2008 (Sarah Rauscher)                                         !
 !                                                                                 !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Necessary input files                                                           !
 !                                                                                 !
 ! temp_list : a file containing the list of temperatures at which data is         !
 !       generated                                                                 !
 ! energy***.xvg : an potential energy file for each T (***)                       !
 ! change variables as necessary (e.g. ntemps, nbins, delta_U, tolerance, nframes) !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! List of variables:                                                              !
 !                                                                                 !
 ! Integers:                                                                       !
 ! i,k,m : counters, i counts over frames, k counts over temperatures, m counts    !
 !     over energy bins                                                            !
 ! temp_index(k) : an array containing the values of the temperatures              !
 ! T(k) : a variable counting temperatures from 1 to 66                            !
 ! nframes : the number of energy values from each temperature                     !
 ! ntemps : the number of temperatures sampled                                     !
 !                                                                                 !
 ! Reals:                                                                          !
 ! time : the time, which is read in the energy***.xvg files (***=temp)            !
 ! potentialE(i,T) : the potential energies read in from the energy***.xvg files,  !
 !     kept in columns for each temperature                                        !
 ! k_b : Boltzmann's constant, 8.3142773x10-3 kJ/mol                               !
 ! beta(k) : the inverse temperature corresponding to T(k)                         !
 ! A(k) : the Helmholtz free energy at temperature T(k), A(k) = f(k)/beta(k)       !  
 ! sum_energy(k), avg_energy(k) :: the sum of all the energy values read in for    !
 !     temperature k, and the average                                              !
 ! first_A_value : the first A value in the set of A values                        !
 !                                                                                 !
 ! Other:                                                                          !
 ! blah : is a character variable which is used in the reading in of the xvg files !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      INTEGER :: i,k,m,nframes,ntemps
      INTEGER, DIMENSION(:), ALLOCATABLE :: T

      INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(p=8)

      REAL (KIND=dbl) :: time 
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: potentialE
      REAL (KIND=dbl) :: k_b
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: beta
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: A

      CHARACTER(LEN=100) :: blah     
      CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: temp_index

      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: sum_energy,avg_energy
      INTEGER :: prev
      REAL (KIND=dbl) :: first_A_value

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Setting the values of the variables to zero or there correct initial values...  !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   i = 0
   k = 0
   m = 0

   nframes = 0
   ntemps = 0
   first_A_value = 0.0_dbl

 ! Reading in the necessary user-specified parameters

      OPEN(1,file='avgE_parameters.prn',form='formatted')

          READ(1,*)blah
          READ(1,*)ntemps,nframes,first_A_value

      CLOSE(1)

 ! The units of the Boltzmann constant are in kJ/mol*K

   k_b = 0.0083142773_dbl

 ! Allocate the necessary arrays

   ALLOCATE (T(ntemps)) 
   ALLOCATE (potentialE(nframes,ntemps))
   ALLOCATE (beta(ntemps))
   ALLOCATE (A(ntemps))
   ALLOCATE (temp_index(ntemps))
   ALLOCATE (sum_energy(ntemps))
   ALLOCATE (avg_energy(ntemps))

 ! Fill in the arrays with the set of temperatures simulated

      OPEN(1,file='temp_list',form='formatted')

          DO k=1,ntemps
             READ(1,*)temp_index(k)
          END DO

      CLOSE(1)

      OPEN(1,file='temp_list',form='formatted')
 
          DO k=1,ntemps
             READ(1,*)T(k)
          END DO

      CLOSE(1)

 ! Calculating the inverse temperatures...

   DO k=1,ntemps
     beta(k) = 1.0_dbl/(k_b*(DBLE(T(k))))
   END DO

 ! Set the histograms and potential energies to zero

    DO i=1,nframes 
       DO k=1,ntemps
          potentialE(i,k) = 0.0_dbl 
       END DO
    END DO

 ! Set the Helmholtz free energies to zero    

   DO k=1,ntemps
      A(k) = 0.0_dbl
   END DO

   DO i=1,ntemps
      sum_energy(i) = 0.0_dbl
      avg_energy(i) = 0.0_dbl
   END DO

   prev = 0

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Generating the potential energy histograms from the simulation data             !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 ! Reading the potential energies from the .xvg files...

    DO k=1,ntemps

      OPEN(1,file='energy'//temp_index(k)//'.xvg',form='formatted') 

          DO i=1,21
             READ(1,*)blah
          END DO

          DO i=1,nframes
             READ(1,*)time,potentialE(i,k)
          END DO

      CLOSE(1)

   END DO

   DO k=1,ntemps
      DO i=1,nframes
         sum_energy(k) = sum_energy(k) + potentialE(i,k)
      END DO
   END DO

   DO k=1,ntemps
      avg_energy(k) = sum_energy(k)/DBLE(nframes)
   END DO

   1024 FORMAT(I3,F20.6)

      OPEN(88,file='AVGE_values.prn',form='formatted')

         DO k=1,ntemps
            WRITE(88,1024)T(k),avg_energy(k)
         END DO

      CLOSE(88)

 ! Calculating the new A values

   A(1) = first_A_value 

   DO k=2,ntemps
      prev = k-1
      A(k) = ((beta(prev)*A(prev)) + (((beta(k)-beta(prev))/2.0_dbl)*(avg_energy(prev) + avg_energy(k))))/beta(k)
   END DO

      OPEN(89,file='new_A_values_park_pande.prn',form='formatted')

         DO k=1,ntemps
            WRITE(89,1024)T(k),A(k)
         END DO

      CLOSE(89)

      END PROGRAM avgE 
