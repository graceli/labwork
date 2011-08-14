      PROGRAM dssp

      IMPLICIT NONE

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 !                                                                                 !
 !                                        DSSP                                     !
 ! Reference : Kabsch, W. and C. Sander. (1983). Dictionary of Protein Secondary   !
 !             Structure: Pattern recognition of Hydrogen-Bonded and Geometrical   !
 !             Features. Biopolymers. 22. 2577-2637.                               !
 !                                                                                 !
 ! Written: Oct. 1, 2005 (Sarah Mansour)                                           !
 ! Commented: Jan. 17, 2008 (Sarah Rauscher)                                       !
 ! Adding allocatable arrays & other cool things: Jan. 31, 2008 (SR)               !
 !                                                                                 !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! List of variables:                                                              !
 !                                                                                 !
 ! Integers:                                                                       !
 ! i,j,t : these are both counters (i,j for the residue # and t for the time frame)!
 ! nres : the total number of residues                                             !
 ! resid : the number of the residue                                               !
 ! nframes : the number of time frames being analysed                              !
 !                                                                                 !
 ! Integer Arrays:                                                                 !
 ! hb_array(i,j,t): an array describing the existence of HBonds between residues i !
 !         and j at frame t, with a 1 if there is and HBond, and 0 if not          !
 ! bend(i,t): an array of integers classifying whether there exists a bend (1) or  !
 !         not (0) at each residue i at each time frame t [bends are only possible !
 !         for i=3,nres-2]                                                         !
 ! chirality(i,t): is the chirality of residue i at time t; the chirality is either!
 !         +1 if the dihedral is between 0 and 180 or -1 if the dihedral is between!
 !         -180 and 0                                                              !
 ! turn_2(i,t): if there is a 2-turn (gamma) at residue i and i+2                  !
 ! turn_3(i,t): if there is a 3-turn (beta) at residue i and i+3                   !
 ! turn_4(i,t): if there is a 4-turn at residue i and i+4                          !
 ! turn_5(i,t): if there is a 5-turn at residue i and i+5                          !
 ! turn_6(i,t): if there is a 6-turn at residue i and i+6                          !
 ! turn_7(i,t): if there is a 7-turn at residue i and i+7                          !
 ! helix_3(i,t): if there is a 3-helix at residue i (3-10 helix)                   !
 ! helix_4(i,t): if there is a 4-helix at residue i (alpha helix)                  !
 ! helix_5(i,t): if there is a 5-helix at residue i (pi helix)                     !
 ! parallel_bridge(i,t): if there is a parallel beta bridge at residue i           !
 ! antiparallel_bridge(i,t): if there is an antiparallel beta bridge               !
 ! parallel_ladder(i,t): if the bridge at i participates in a parallel ladder      !
 ! antiparallel_ladder(i,t): if the bridge participates in an antiparallel ladder  !
 ! all the variables starting with "final" are the corresponding values after      !
 ! including the hierarchy of structures***                                        !
 ! final_bend(i,t): bends after assigning only one structure to each residue       !
 ! final_turn_2(i,t): 2-turns after assigning only one structure to each residue   !
 ! final_turn_3(i,t): 3-turns after assigning only one structure to each residue   !
 ! final_turn_4(i,t): 4-turns after assigning only one structure to each residue   !
 ! final_turn_5(i,t): 5-turns after assigning only one structure to each residue   !
 ! final_turn_6(i,t): 6-turns after assigning only one structure to each residue   !
 ! final_turn_7(i,t): 7-turns after assigning only one structure to each residue   !
 ! final_helix_3(i,t): 3-10 helix after assigning only one structure to each res   !
 ! final_helix_4(i,t): alpha helix after assigning only one structure to each res  !
 ! final_helix_5(i,t): pi helix after assigning only one structure to each res     !
 ! final_parallel_bridge(i,t): ...                                                 !
 ! final_antiparallel_bridge(i,t): ...                                             !
 ! final_parallel_ladder(i,t): ...                                                 !
 ! final_antiparallel_ladder(i,t): ...                                             !
 !                                                                                 !
 ! Real Numbers:                                                                   !
 !                                                                                 !
 ! Real Arrays:                                                                    !
 ! Cx(i),Cy(i),Cz(i): the coordinates of the C atoms (C=O C) for each residue i    !
 ! Ox(i),Oy(i),Oz(i): the coordinates of the O atoms (C=O O) for each residue i    !
 ! Nx(i),Ny(i),Nz(i): the coordinates of the N atoms (N-H N) for each residue i    !
 ! Hx(i),Hy(i),Hz(i): the coordinates of the H atoms (N-H H) for each residue i    !
 ! Cax(i),Cay(i),Caz(i): the coordinates of the C alpha atoms for each residue i   !
 ! dON(i,j): the distance between the O of residue i and the N of residue j        !
 ! dCH(i,j): the distance between the C of residue i and the H of residue j        !
 ! dOH(i,j): the distance between the O of residue i and the H of residue j        !
 ! dCN(i,j): the distance between the C of residue i and the N of residue j        !
 ! E_hb(i,j): the energy of the HB between residue i and residue j                 !
 ! vector_a(i,3): the vector from Ca(i-2) to Ca(i) for each residue i              !
 ! vector_b(i,3): the vector from Ca(i) to Ca(i+2) for each residue i              !
 ! magnitude_a(i): the magnitude of vector_a for residue i                         !
 ! magnitude_b(i): the magnitude of vector_b for residue i                         !
 ! a_dot_b(i): the dot product of vectors a and b for residue i                    !
 ! cos_theta(i): the cosine of the angle between vectors a and b                   !
 ! normal_1(i,3): is the normal to the plane of C_alpha's i-1, i and i+1           !
 ! normal_2(i,3): is the normal to the plane of C_alpha's i, i+1, and i+2          !
 ! normal_1_dot_normal_2(i): is the dot product of normal_1 and normal_2           !
 ! magnitude_normal_1(i): is the magnitude of normal_1                             !
 ! magnitude_normal_2(i): is the magnitude of normal_2                             !
 ! cos_chi(i): is the cosine of the angle between the planes (the dihedral)        !
 ! chi(i): is the dihedral angle (in radians)                                      !
 ! chi_deg(i): is the dihedral angle (in degrees)                                  !
 !                                                                                 !
 ! Character Variables:                                                            !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Defining variables and arrays                                                   !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! Defining all integer variables and arrays

      INTEGER :: i,j,t,b,d,k,l
      INTEGER :: nres
      INTEGER :: resid
      INTEGER :: nframes

 ! The size of the arrays is now dynamically allocated 

      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: hb_array
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: bend
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: chirality
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: turn_1,turn_2,turn_3,turn_4,turn_5,turn_6,turn_7
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: turn_1_II,turn_2_II,turn_3_II                                              
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: turn_4_II,turn_5_II,turn_6_II,turn_7_II    
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: helix_3,helix_4,helix_5
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: parallel_bridge,antiparallel_bridge
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: parallel_ladder,antiparallel_ladder 
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: final_bend
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: final_turn_1,final_turn_2,final_turn_3
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: final_turn_4,final_turn_5,final_turn_6,final_turn_7
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: final_turn_1_II,final_turn_2_II,final_turn_3_II
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: final_turn_4_II,final_turn_5_II,final_turn_6_II,final_turn_7_II
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: final_helix_3,final_helix_4,final_helix_5
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: final_parallel_bridge,final_antiparallel_bridge
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: final_parallel_ladder,final_antiparallel_ladder 

 ! These variables keep track of the total amount of each structure over all the snapshots

      INTEGER :: total_hbonds
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_hbonds_at_t
      INTEGER :: total_bends
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_bends_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_bends_at_i

      INTEGER :: total_turn_1
      INTEGER :: total_turn_2
      INTEGER :: total_turn_3
      INTEGER :: total_turn_4
      INTEGER :: total_turn_5
      INTEGER :: total_turn_6
      INTEGER :: total_turn_7
      INTEGER :: total_helix_3
      INTEGER :: total_helix_4
      INTEGER :: total_helix_5
      INTEGER :: total_parallel_bridge
      INTEGER :: total_antiparallel_bridge
      INTEGER :: total_parallel_ladder
      INTEGER :: total_antiparallel_ladder
      INTEGER :: total_turn_1_II
      INTEGER :: total_turn_2_II
      INTEGER :: total_turn_3_II
      INTEGER :: total_turn_4_II
      INTEGER :: total_turn_5_II
      INTEGER :: total_turn_6_II
      INTEGER :: total_turn_7_II

      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_1_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_2_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_3_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_4_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_5_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_6_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_7_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_1_II_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_2_II_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_3_II_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_4_II_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_5_II_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_6_II_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_7_II_at_t

      INTEGER, DIMENSION(:), ALLOCATABLE :: total_helix_3_at_t 
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_helix_4_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_helix_5_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_parallel_bridge_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_antiparallel_bridge_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_parallel_ladder_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_antiparallel_ladder_at_t

      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_1_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_2_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_3_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_4_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_5_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_6_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_7_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_1_II_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_2_II_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_3_II_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_4_II_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_5_II_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_6_II_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_turn_7_II_at_i

      INTEGER, DIMENSION(:), ALLOCATABLE :: total_helix_3_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_helix_4_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_helix_5_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_parallel_bridge_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_antiparallel_bridge_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_parallel_ladder_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_antiparallel_ladder_at_i

 ! Incorporating PPII
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PP2
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PP2_only1s

 ! This variable tells if there is PPII after accounting for the hierarchy of structure
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PP2_final

 ! These count the total PPII, for each frame, for each residue, and over all snapshots
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_PP2_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_PP2_at_i
      INTEGER :: total_PP2

 ! This variable is used to define double precision

      INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(p=8)

 ! These variables contain the coordinates of the backbone atoms for each residue 
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: Cx,Cy,Cz
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: Ox,Oy,Oz
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: Nx,Ny,Nz
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: Hx,Hy,Hz
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: Cax,Cay,Caz

 ! These variables are used to determine the presence of peptide-peptide HBs
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dON
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dCH
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dOH
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dCN
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dHN
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dCO
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: q
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: costheta1
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: E_hb
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: vector_a,vector_b
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: magnitude_a,magnitude_b
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: a_dot_b,cos_theta
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: normal_1,normal_2
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: normal_1_dot_normal_2
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: magnitude_normal_1
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: magnitude_normal_2
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: cos_chi,chi,chi_deg 

 ! These variables are used to calculate phi and psi at each snapshot

      REAL :: torsion

      REAL, DIMENSION(:,:), ALLOCATABLE :: phi
      REAL, DIMENSION(:,:), ALLOCATABLE :: psi

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: phi_int,psi_int

      REAL :: Cim1(3)
      REAL :: Ni(3)
      REAL :: Cai(3)
      REAL :: Ci(3)
      REAL :: Nip1(3)

 ! These variables are needed to write out the contact map

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: hb_array_total
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: hb_array_total_normalized
      INTEGER, DIMENSION(:), ALLOCATABLE :: total_hbonds_at_res_i

 ! These variables are needed to keep track of the end to end distance

      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: eed_at_t
      INTEGER, DIMENSION(:), ALLOCATABLE :: nint_eed_at_t
      INTEGER :: eed_dist(400)
      REAL (KIND=dbl) :: eed_normalized(400) 
      REAL (KIND=dbl) :: eed_list,eed_avg,eed_sum

 ! These variables keep track of the Ca3, Ca4 and Ca5 distances

      REAL, DIMENSION(:,:), ALLOCATABLE :: dCa3,dCa4,dCa5
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: nint_dCa3,nint_dCa4,nint_dCa5
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: dCa3_dist,dCa4_dist,dCa5_dist
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dCa3_normalized
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dCa4_normalized
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: dCa5_normalized

 ! These variables are used to output normalized values of the average properties

      REAL (KIND=dbl) :: total_hbonds_normalized,total_bends_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_bends_at_i_normalized 
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_hbonds_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_PP2_at_i_normalized
      REAL (KIND=dbl) :: total_PP2_normalized

      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_2_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_3_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_4_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_5_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_6_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_7_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_2_II_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_3_II_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_4_II_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_5_II_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_6_II_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_turn_7_II_at_i_normalized

      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_helix_3_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_helix_4_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_helix_5_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_parallel_bridge_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_antiparallel_bridge_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_parallel_ladder_at_i_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: total_antiparallel_ladder_at_i_normalized

      REAL (KIND=dbl) :: total_turn_2_normalized
      REAL (KIND=dbl) :: total_turn_3_normalized
      REAL (KIND=dbl) :: total_turn_4_normalized
      REAL (KIND=dbl) :: total_turn_5_normalized
      REAL (KIND=dbl) :: total_turn_6_normalized
      REAL (KIND=dbl) :: total_turn_7_normalized
      REAL (KIND=dbl) :: total_turn_2_II_normalized
      REAL (KIND=dbl) :: total_turn_3_II_normalized
      REAL (KIND=dbl) :: total_turn_4_II_normalized
      REAL (KIND=dbl) :: total_turn_5_II_normalized
      REAL (KIND=dbl) :: total_turn_6_II_normalized
      REAL (KIND=dbl) :: total_turn_7_II_normalized

      REAL (KIND=dbl) :: total_helix_3_normalized
      REAL (KIND=dbl) :: total_helix_4_normalized
      REAL (KIND=dbl) :: total_helix_5_normalized
      REAL (KIND=dbl) :: total_parallel_bridge_normalized
      REAL (KIND=dbl) :: total_antiparallel_bridge_normalized
      REAL (KIND=dbl) :: total_parallel_ladder_normalized
      REAL (KIND=dbl) :: total_antiparallel_ladder_normalized

 ! These are just variables created in order to read the .pdb files

      REAL :: m,n

      CHARACTER(LEN=100) :: blah 
      CHARACTER(LEN=4) :: atom
      CHARACTER(LEN=1) :: c
      CHARACTER(LEN=3) :: res

 ! Keeping track of the number of prolines

      INTEGER :: nprolines
      INTEGER, DIMENSION(:), ALLOCATABLE :: proline_resnumbers

 ! These variables are for water

      INTEGER :: nwat

 ! These arrays contain the coordinates of all the waters at each time frame

      REAL, DIMENSION(:,:,:), ALLOCATABLE :: Owat,H1wat,H2wat

 ! These variables are used in the calculation of hydration
 
      REAL, DIMENSION(:,:), ALLOCATABLE :: dN_Owat,dO_Owat
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: HB_N,HB_CO
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: HB_N_total_at_i,HB_CO_total_at_i,HB_total_at_i
      INTEGER, DIMENSION(:), ALLOCATABLE :: HB_total_at_t
      INTEGER :: HBtotal
      INTEGER :: inow,jnow,tnow

      REAL :: dOwat_H
      REAL :: dHN_wat
      REAL :: costheta
      REAL :: dH1wat_H
      REAL :: dH1wat_N
      REAL :: dH2wat_H
      REAL :: dH2wat_N
      REAL :: E_hb_H1
      REAL :: E_hb_H2

      REAL :: dO_H1wat
      REAL :: dO_H2wat
      REAL :: dOwat_H1wat
      REAL :: dOwat_H2wat
      REAL :: dC_H1wat
      REAL :: dC_H2wat
      REAL :: dC_Owat
      REAL :: E_hb_wat

      INTEGER :: one,two
      INTEGER :: two_x_nres

      INTEGER, DIMENSION(:), ALLOCATABLE :: HB_total_at_resi
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: HB_both
      INTEGER, DIMENSION(:), ALLOCATABLE :: HB_both_sum
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: HB_bridging
      INTEGER, DIMENSION(:), ALLOCATABLE :: nbound_t,nbridging_t
 
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: HB_bridging_sum
      INTEGER :: nbound_sum,nbridging_sum
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: HB_bridging_avg 
      REAL (KIND=dbl) :: nbound_avg,nbridging_avg

      INTEGER :: nbound_dist(0:40),nbridging_dist(0:40)
      REAL (KIND=dbl) :: nbound_dist_normalized(0:40),nbridging_dist_normalized(0:40)   

      ! Variables for counting bridging waters with both N-H and C=O                

      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: HB_bridging_NH_CO
      INTEGER :: int_1,int_2
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: HB_bridging_NH_CO_sum
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: HB_bridging_NH_CO_avg

      ! Variables for writing normalized values of the important hydration parameters

      REAL (KIND=dbl) :: HBtotal_normalized
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: HB_total_at_resi_normalized

 ! These variables keep track of the desired output files
 
      INTEGER :: hierarchy,time_evolution,phi_psi,contact_map,eed,turn_locations
      INTEGER :: inter_Ca_dist,water

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Initialize all variables to zero or their starting values                       !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      i = 0
      j = 0
      t = 0 
      b = 0
      d = 0
      k = 0
      l = 0
      resid = 0
      nres = 0
      nframes = 0

      hierarchy = 0
      time_evolution = 0
      phi_psi = 0
      contact_map = 0
      eed = 0
      turn_locations = 0
      inter_Ca_dist = 0
      water = 0

      OPEN(1,file='desired_output_files.prn',form='formatted')

          READ(1,*)blah
          READ(1,*)hierarchy,time_evolution,phi_psi,contact_map,eed,turn_locations,&
                   &inter_Ca_dist,water

      CLOSE(1)

      OPEN(1,file='nres_nframes_nprolines_nwat.prn',form='formatted')

         READ(1,*)nres,nframes,nprolines

      CLOSE(1)

      IF (water.eq.1) THEN

      OPEN(1,file='nres_nframes_nprolines_nwat.prn',form='formatted')

         READ(1,*)nres,nframes,nprolines,nwat

      CLOSE(1)

      END IF

      two_x_nres = 2*nres

      IF (nprolines > 0) THEN

         ALLOCATE (proline_resnumbers(nprolines))

         OPEN(1,file='proline_resnumbers.prn',form='formatted')

             DO i=1,nprolines
                READ(1,*)proline_resnumbers(i)
             END DO

         CLOSE(1) 
      
      END IF 

      ALLOCATE (Cx(nres,nframes))
      ALLOCATE (Cy(nres,nframes))
      ALLOCATE (Cz(nres,nframes))
      ALLOCATE (Ox(nres,nframes))
      ALLOCATE (Oy(nres,nframes))
      ALLOCATE (Oz(nres,nframes))
      ALLOCATE (Nx(nres,nframes))
      ALLOCATE (Ny(nres,nframes))
      ALLOCATE (Nz(nres,nframes))
      ALLOCATE (Hx(nres,nframes))
      ALLOCATE (Hy(nres,nframes))
      ALLOCATE (Hz(nres,nframes))  
      ALLOCATE (Cax(nres,nframes))
      ALLOCATE (Cay(nres,nframes))
      ALLOCATE (Caz(nres,nframes))  

      DO i=1,nres
      DO t=1,nframes
         Cx(i,t) = 0.0_dbl
         Cy(i,t) = 0.0_dbl
         Cz(i,t) = 0.0_dbl
         Ox(i,t) = 0.0_dbl
         Oy(i,t) = 0.0_dbl
         Oz(i,t) = 0.0_dbl
         Nx(i,t) = 0.0_dbl
         Ny(i,t) = 0.0_dbl
         Nz(i,t) = 0.0_dbl
         Hx(i,t) = 0.0_dbl
         Hy(i,t) = 0.0_dbl
         Hz(i,t) = 0.0_dbl
         Cax(i,t) = 0.0_dbl
         Cay(i,t) = 0.0_dbl
         Caz(i,t) = 0.0_dbl
      END DO
      END DO
 
      ALLOCATE (dON(nres,nres))
      ALLOCATE (dCH(nres,nres))
      ALLOCATE (dOH(nres,nres))
      ALLOCATE (dCN(nres,nres))
      ALLOCATE (dHN(nres,nres))
      ALLOCATE (dCO(nres,nres))
      ALLOCATE (E_hb(nres,nres))
      ALLOCATE (hb_array(nres,nres,nframes))
      ALLOCATE (q(nres,nres))
      ALLOCATE (costheta1(nres,nres))

      DO i=1,nres
         DO j=1,nres
                 dON(i,j) = 0.0_dbl
                 dCH(i,j) = 0.0_dbl
                 dOH(i,j) = 0.0_dbl
                 dCN(i,j) = 0.0_dbl
                 E_hb(i,j) = 0.0_dbl
                 dHN(i,j) = 0.0_dbl
                 dCO(i,j) = 0.0_dbl
                 q(i,j) = 0.0_dbl
                 costheta1(i,j) = 0.0_dbl
         END DO
      END DO
 
      DO i=1,nres
         DO j=1,nres
            DO t=1,nframes
               hb_array(i,j,t) = 0
            END DO
         END DO
      END DO

      ALLOCATE (vector_a(nres,3))
      ALLOCATE (vector_b(nres,3))
      ALLOCATE (magnitude_a(nres))
      ALLOCATE (magnitude_b(nres))
      ALLOCATE (a_dot_b(nres))
      ALLOCATE (cos_theta(nres))
      ALLOCATE (normal_1(nres,3))
      ALLOCATE (normal_2(nres,3))
      ALLOCATE (normal_1_dot_normal_2(nres))
      ALLOCATE (magnitude_normal_1(nres))
      ALLOCATE (magnitude_normal_2(nres))
      ALLOCATE (cos_chi(nres))
      ALLOCATE (chi(nres))
      ALLOCATE (chi_deg(nres))
      ALLOCATE (bend(nres,nframes))
      ALLOCATE (chirality(nres,nframes))

      DO i=1,nres
         vector_a(i,1) = 0.0_dbl
         vector_a(i,2) = 0.0_dbl
         vector_a(i,3) = 0.0_dbl
         vector_b(i,1) = 0.0_dbl
         vector_b(i,2) = 0.0_dbl
         vector_b(i,3) = 0.0_dbl
         magnitude_a(i) = 0.0_dbl
         magnitude_b(i) = 0.0_dbl
         a_dot_b(i) = 0.0_dbl
         cos_theta(i) = 0.0_dbl
         normal_1(i,1) = 0.0_dbl
         normal_1(i,2) = 0.0_dbl
         normal_1(i,3) = 0.0_dbl
         normal_2(i,1) = 0.0_dbl
         normal_2(i,2) = 0.0_dbl
         normal_2(i,3) = 0.0_dbl
         normal_1_dot_normal_2(i) = 0.0_dbl
         magnitude_normal_1(i) = 0.0_dbl
         magnitude_normal_2(i) = 0.0_dbl
         cos_chi(i) = 0.0_dbl
         chi(i) = 0.0_dbl
         chi_deg(i) = 0.0_dbl
       END DO

      DO t=1,nframes
         DO i=1,nres
             bend(i,t) = 0
             chirality(i,t) = 0
         END DO
      END DO

      ALLOCATE (phi(nres,nframes))
      ALLOCATE (psi(nres,nframes))
      ALLOCATE (phi_int(nres,nframes))
      ALLOCATE (psi_int(nres,nframes))

      DO i=1,nres
         DO t=1,nframes
            phi(i,t) = 0.0
            psi(i,t) = 0.0
            phi_int(i,t) = 0
            psi_int(i,t) = 0
         END DO
      END DO

      DO i=1,3
         Cim1(i) = 0.0
         Ni(i) = 0.0
         Cai(i) = 0.0
         Ci(i) = 0.0
         Nip1(i) = 0.0 
      END DO

      ALLOCATE (eed_at_t(nframes))

      DO t=1,nframes
         eed_at_t(t) = 0.0_dbl
      END DO

      IF (inter_Ca_dist.eq.1) THEN

      ALLOCATE (dCa3(nres-2,nframes))
      ALLOCATE (dCa4(nres-3,nframes))
      ALLOCATE (dCa5(nres-4,nframes))

      DO i=1,nres-2
         DO t=1,nframes
            dCa3(i,t) = 0.00000
         END DO
      END DO

      DO i=1,nres-3
         DO t=1,nframes
            dCa4(i,t) = 0.00000
         END DO
      END DO

      DO i=1,nres-4
         DO t=1,nframes
            dCa5(i,t) = 0.00000
         END DO
      END DO

      END IF

      total_hbonds_normalized = 0.0_dbl
      total_bends_normalized = 0.0_dbl

      ALLOCATE (total_bends_at_i_normalized(nres))
      ALLOCATE (total_hbonds_at_i_normalized(nres))
      ALLOCATE (total_PP2_at_i_normalized(nres))

      DO i=1,nres
         total_bends_at_i_normalized(i) = 0.0_dbl
         total_hbonds_at_i_normalized(i) = 0.0_dbl
         total_PP2_at_i_normalized(i) = 0.0_dbl
      END DO

      total_PP2_normalized = 0.0_dbl

      total_turn_2_normalized = 0.0_dbl
      total_turn_3_normalized = 0.0_dbl
      total_turn_4_normalized = 0.0_dbl
      total_turn_5_normalized = 0.0_dbl
      total_turn_6_normalized = 0.0_dbl
      total_turn_7_normalized = 0.0_dbl
      total_turn_2_II_normalized = 0.0_dbl
      total_turn_3_II_normalized = 0.0_dbl
      total_turn_4_II_normalized = 0.0_dbl
      total_turn_5_II_normalized = 0.0_dbl
      total_turn_6_II_normalized = 0.0_dbl
      total_turn_7_II_normalized = 0.0_dbl

      total_helix_3_normalized = 0.0_dbl
      total_helix_4_normalized = 0.0_dbl
      total_helix_5_normalized = 0.0_dbl
      total_parallel_bridge_normalized = 0.0_dbl
      total_antiparallel_bridge_normalized = 0.0_dbl
      total_parallel_ladder_normalized = 0.0_dbl
      total_antiparallel_ladder_normalized = 0.0_dbl

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Read in the C,O,N and H coordinates and calculate the relevant distances for    !
 ! HBonds for each time frame                                                      !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   OPEN(1,file='outc.pdb',form='formatted')

        DO t=1,nframes
  
           DO i=1,5
              READ(1,*)blah
           END DO
                 DO i=1,nres
                      READ(1,*)atom,b,c,res,d,Cx(i,t),Cy(i,t),Cz(i,t),m,n
                 END DO
           DO i=1,2
              READ(1,*)blah        
           END DO      

        END DO

   CLOSE(1)

   OPEN(2,file='outo.pdb',form='formatted')

        DO t=1,nframes
 
           DO i=1,5
              READ(2,*)blah        
           END DO      
                 DO i=1,nres
                      READ(2,*)atom,b,c,res,d,Ox(i,t),Oy(i,t),Oz(i,t),m,n
                 END DO
           DO i=1,2
              READ(2,*)blah
           END DO

        END DO

   CLOSE(2)

   OPEN(3,file='outn.pdb',form='formatted')

        DO t=1,nframes

           DO i=1,5
              READ(3,*)blah
           END DO
                 DO i=1,nres
                      READ(3,*)atom,b,c,res,d,Nx(i,t),Ny(i,t),Nz(i,t),m,n
                 END DO
           DO i=1,2
              READ(3,*)blah
           END DO

        END DO

   CLOSE(3)

   OPEN(4,file='outh.pdb',form='formatted')

        DO t=1,nframes

           DO i=1,5
              READ(4,*)blah
           END DO
                 DO i=1,nres
                      READ(4,*)atom,b,c,res,d,Hx(i,t),Hy(i,t),Hz(i,t),m,n
                 END DO
           DO i=1,2
              READ(4,*)blah
           END DO

        END DO

   CLOSE(4)

   OPEN(5,file='outca.pdb',form='formatted')

        DO t=1,nframes

           DO i=1,5
              READ(5,*)blah
           END DO
                 DO i=1,nres
                      READ(5,*)atom,b,c,res,d,Cax(i,t),Cay(i,t),Caz(i,t),m,n
                 END DO
           DO i=1,2
              READ(5,*)blah
           END DO

        END DO

   CLOSE(5)

 ! Calculate the relevant parameters for each time frame

 DO t=1,nframes

 ! Calculating phi and psi

         DO i=2,nres-1
            Cim1(1)=Cx(i-1,t)
            Cim1(2)=Cy(i-1,t)
            Cim1(3)=Cz(i-1,t)
            Ni(1)=Nx(i,t)
            Ni(2)=Ny(i,t)
            Ni(3)=Nz(i,t)
            Cai(1)=Cax(i,t)
            Cai(2)=Cay(i,t)
            Cai(3)=Caz(i,t)
            Ci(1)=Cx(i,t)
            Ci(2)=Cy(i,t)
            Ci(3)=Cz(i,t)
            Nip1(1)=Nx(i+1,t)
            Nip1(2)=Ny(i+1,t)
            Nip1(3)=Nz(i+1,t)
            phi(i,t) = torsion ( Cim1, Ni, Cai, Ci )
            psi(i,t) = torsion ( Ni, Cai, Ci, Nip1 )
         END DO

 ! Calculate the end to end distance

   eed_at_t(t) = SQRT(((Cax(nres,t)-Cax(1,t))**2) + ((Cay(nres,t)-Cay(1,t))**2) + ((Caz(nres,t)-Caz(1,t))**2))

 IF (inter_Ca_dist.eq.1) THEN

 ! Calculate the inter Ca distances

   DO i=1,nres-2
      j=i+2
      dCa3(i,t) = SQRT(((Cax(i,t)-Cax(j,t))**2) + ((Cay(i,t)-Cay(j,t))**2) + ((Caz(i,t)-Caz(j,t))**2))
   END DO

   DO i=1,nres-3
      j=i+3
      dCa4(i,t) = SQRT(((Cax(i,t)-Cax(j,t))**2) + ((Cay(i,t)-Cay(j,t))**2) + ((Caz(i,t)-Caz(j,t))**2))
   END DO

   DO i=1,nres-4
      j=i+4
      dCa5(i,t) = SQRT(((Cax(i,t)-Cax(j,t))**2) + ((Cay(i,t)-Cay(j,t))**2) + ((Caz(i,t)-Caz(j,t))**2))
   END DO

 END IF
 
 ! Calculate the relevant distances to determine if there are HBs

                 DO i=1,nres
                      DO j=1,nres
                          dON(i,j) = DBLE(SQRT(((Nx(j,t)-Ox(i,t))**2) + &
                                       & ((Ny(j,t)-Oy(i,t))**2) + ((Nz(j,t)-Oz(i,t))**2)))
                          dCH(i,j) = DBLE(SQRT(((Hx(j,t)-Cx(i,t))**2) + &
                                       & ((Hy(j,t)-Cy(i,t))**2) + ((Hz(j,t)-Cz(i,t))**2)))
                          dOH(i,j) = DBLE(SQRT(((Hx(j,t)-Ox(i,t))**2) + &
                                       & ((Hy(j,t)-Oy(i,t))**2) + ((Hz(j,t)-Oz(i,t))**2)))
                          dCN(i,j) = DBLE(SQRT(((Nx(j,t)-Cx(i,t))**2) + &
                                       & ((Ny(j,t)-Cy(i,t))**2) + ((Nz(j,t)-Cz(i,t))**2)))
                          dHN(i,j) = DBLE(SQRT(((Nx(j,t)-Hx(j,t))**2) + &
                                       & ((Ny(j,t)-Hy(j,t))**2) + ((Nz(j,t)-Hz(j,t))**2)))
                          dCO(i,j) = DBLE(SQRT(((Cx(i,t)-Ox(i,t))**2) + &
                                       & ((Cy(i,t)-Oy(i,t))**2) + ((Cz(i,t)-Oz(i,t))**2)))
                      END DO
                 END DO

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Creating the array of energies of potential H-Bonds and the H-Bond array        !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! Calculate the energy of the HB

         DO i=1,nres
            DO j=1,nres
               E_hb(i,j) = (2.78746712_dbl/dON(i,j)) + (2.78746712_dbl/dCH(i,j)) &
                                 & - (2.78746712_dbl/dOH(i,j)) - (2.78746712_dbl/dCN(i,j))
            END DO
         END DO

 !  Calculate the angle criteria for HB's:

         DO i=1,nres
            DO j=1,nres
               q(i,j) = DBLE(((dON(i,j)**2) - (dOH(i,j)**2) - (dHN(i,j)**2))/(2*dHN(i,j)))
               costheta1(i,j) = DBLE(q(i,j)/dOH(i,j))
            END DO
         END DO

 ! Using a more strict definition of a hydrogren bond...
 ! 1. The energy must be less than E=-0.5 kcal/mol
 ! 2. The distance between O and H must be less than 2.5 angstroms (Torshin et al 2002)
 ! 3. The distance between O and N must be less than 3.5 angstroms
 ! 4. The angle theta1 must be less than 60 degrees

         DO i=1,nres
            DO j=1,nres
               IF (E_hb(i,j).lt.-0.05) THEN
                  IF (dON(i,j).lt.3.50) THEN
                     IF (dOH(i,j).lt.2.50) THEN
                        IF (costheta1(i,j).gt.0.5) THEN
                                  hb_array(i,j,t) = 1
                        END IF
                     END IF
                  END IF
               END IF
            END DO
         END DO

         DO i=1,nres
            DO j=1,nres
               IF (dOH(i,j).gt.2.50) THEN
                  hb_array(i,j,t) = 0
               END IF
               IF (dON(i,j).gt.3.50) THEN
                  hb_array(i,j,t) = 0
               END IF
               IF (costheta1(i,j).lt.0.5) THEN
                  hb_array(i,j,t) = 0
               END IF
            END DO
         END DO


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Read in the C_alpha coordinates and calculate the vectors to get cos_theta(i),  !
 ! which will be used to classify if there is a bend at each residue i then use    !
 ! the C_alpha coordinates to determine the chirality at position i                !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate cos_theta first to see if there is a "bend" at residue i

                 DO i=3,nres-2
                      vector_a(i,1) = DBLE(Cax(i,t) - Cax(i-2,t))
                      vector_a(i,2) = DBLE(Cay(i,t) - Cay(i-2,t))
                      vector_a(i,3) = DBLE(Caz(i,t) - Caz(i-2,t))
                      vector_b(i,1) = DBLE(Cax(i+2,t) - Cax(i,t))
                      vector_b(i,2) = DBLE(Cay(i+2,t) - Cay(i,t))
                      vector_b(i,3) = DBLE(Caz(i+2,t) - Caz(i,t))
                 END DO

                 DO i=3,nres-2
                      magnitude_a(i) = DBLE(SQRT((vector_a(i,1)**2) + &
                                & (vector_a(i,2)**2) + (vector_a(i,3)**2)))
                      magnitude_b(i) = DBLE(SQRT((vector_b(i,1)**2) + &
                                & (vector_b(i,2)**2) + (vector_b(i,3)**2)))
                 END DO

                 DO i=3,nres-2
                      a_dot_b(i) = DBLE((vector_a(i,1)*vector_b(i,1)) + &
                                & (vector_a(i,2)*vector_b(i,2)) + &
                                & (vector_a(i,3)*vector_b(i,3)))
                      cos_theta(i) = DBLE(a_dot_b(i)/(magnitude_a(i)*magnitude_b(i)))
                 END DO

                 DO i=3,nres-2
                      IF (cos_theta(i).lt.0.342_dbl) THEN
                             bend(i,t) = 1
                      ELSE
                             bend(i,t) = 0
                      END IF
                 END DO

! Calculate the dihedral angle between Ca's to determine "chirality"

                 DO i=2,nres-2
                      normal_1(i,1) = ((Cay(i,t) - Cay(i-1,t))*(Caz(i+1,t) - Caz(i,t))) - &
                                & ((Caz(i,t) - Caz(i-1,t))*(Cay(i+1,t) - Cay(i,t)))
                      normal_1(i,2) = - ((Cax(i,t) - Cax(i-1,t))*(Caz(i+1,t) - Caz(i,t))) + &
                                & ((Caz(i,t) - Caz(i-1,t))*(Cax(i+1,t) - Cax(i,t)))
                      normal_1(i,3) = ((Cax(i,t) - Cax(i-1,t))*(Cay(i+1,t) - Cay(i,t))) - &
                                & ((Cay(i,t) - Cay(i-1,t))*(Cax(i+1,t) - Cax(i,t)))
                      normal_2(i,1) = ((Cay(i+1,t) - Cay(i,t))*(Caz(i+2,t) - Caz(i+1,t))) - &
                                & ((Caz(i+1,t) - Caz(i,t))*(Cay(i+2,t) - Cay(i+1,t)))
                      normal_2(i,2) = - ((Cax(i+1,t) - Cax(i,t))*(Caz(i+2,t) - Caz(i+1,t))) + &
                                & ((Caz(i+1,t) - Caz(i,t))*(Cax(i+2,t) - Cax(i+1,t)))
                      normal_2(i,3) = ((Cax(i+1,t) - Cax(i,t))*(Cay(i+2,t) - Cay(i+1,t))) - &
                                & ((Cay(i+1,t) - Cay(i,t))*(Cax(i+2,t) - Cax(i+1,t)))
                 END DO

                 DO i=2,nres-2
                      normal_1_dot_normal_2(i) = (normal_1(i,1)*normal_2(i,1)) + &
                                & (normal_1(i,2)*normal_2(i,2)) + &
                                & (normal_1(i,3)*normal_2(i,3))
                      magnitude_normal_1(i) = SQRT((normal_1(i,1)**2) + &
                                & (normal_1(i,2)**2) + (normal_1(i,3)**2))
                      magnitude_normal_2(i) = SQRT((normal_2(i,1)**2) + &
                                & (normal_2(i,2)**2) + (normal_2(i,3)**2))
                 END DO

                 DO i=2,nres-2
                      cos_chi(i) = normal_1_dot_normal_2(i)/(magnitude_normal_1(i)*&
                                &magnitude_normal_2(i))
                      chi(i) = ACOS(cos_chi(i))
                      chi_deg(i) = chi(i)*57.2957795
                 END DO

                 DO i=2,nres-2
                      IF (chi_deg(i).lt.180.0) THEN
                         IF (chi_deg(i).gt.0.0) THEN
                            chirality(i,t) = 1
                         ELSE IF (chi_deg(i).lt.0.0) THEN
                            IF (chi_deg(i).gt.-180.0) THEN
                               chirality(i,t) = -1
                            END IF
                         END IF
                      END IF
                 END DO

 ! Resetting the necessary variables to zero

      DO i=1,nres
         DO j=1,nres
                 dON(i,j) = 0.0_dbl
                 dCH(i,j) = 0.0_dbl
                 dOH(i,j) = 0.0_dbl
                 dCN(i,j) = 0.0_dbl
                 E_hb(i,j) = 0.0_dbl
                 dHN(i,j) = 0.0_dbl
                 dCO(i,j) = 0.0_dbl
                 q(i,j) = 0.0_dbl
                 costheta1(i,j) = 0.0_dbl
         END DO
      END DO

      DO i=1,nres
         vector_a(i,1) = 0.0_dbl
         vector_a(i,2) = 0.0_dbl
         vector_a(i,3) = 0.0_dbl
         vector_b(i,1) = 0.0_dbl
         vector_b(i,2) = 0.0_dbl
         vector_b(i,3) = 0.0_dbl
         magnitude_a(i) = 0.0_dbl
         magnitude_b(i) = 0.0_dbl
         a_dot_b(i) = 0.0_dbl
         cos_theta(i) = 0.0_dbl
         normal_1(i,1) = 0.0_dbl
         normal_1(i,2) = 0.0_dbl
         normal_1(i,3) = 0.0_dbl
         normal_2(i,1) = 0.0_dbl
         normal_2(i,2) = 0.0_dbl
         normal_2(i,3) = 0.0_dbl
         normal_1_dot_normal_2(i) = 0.0_dbl
         magnitude_normal_1(i) = 0.0_dbl
         magnitude_normal_2(i) = 0.0_dbl
         cos_chi(i) = 0.0_dbl
         chi(i) = 0.0_dbl
         chi_deg(i) = 0.0_dbl
       END DO

       DO i=1,3
            Cim1(i)=0.0
            Ni(i)=0.0
            Cai(i)=0.0
            Ci(i)=0.0
            Nip1(i)=0.0
       END DO

      END DO

 ! This ends the cycle over all the frames

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Make sure there are no miscounted HBs within a peptide group, or with prolines  ! 
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! Remove the possibility of an HB within a peptide group

     DO t=1,nframes
         DO i=1,nres
            DO j=1,nres
               IF (i.eq.j) THEN
                  hb_array(i,j,t) = 0
               END IF
            END DO
         END DO
     END DO

 ! prolines (if there is any, make sure that they don't have any HBs with an N-H) 

      IF ( nprolines > 0 ) THEN

         DO i=1,nprolines

         DO t=1,nframes
            j=proline_resnumbers(i)
               DO k=1,nres
                  hb_array(k,j,t) = 0
               END DO
         END DO
   
         END DO

      END IF

      DEALLOCATE (dON)
      DEALLOCATE (dCH)
      DEALLOCATE (dOH)
      DEALLOCATE (dCN)
      DEALLOCATE (dHN)
      DEALLOCATE (dCO)
      DEALLOCATE (E_hb)
      DEALLOCATE (q)
      DEALLOCATE (costheta1)

      DEALLOCATE (vector_a)
      DEALLOCATE (vector_b)
      DEALLOCATE (magnitude_a)
      DEALLOCATE (magnitude_b)
      DEALLOCATE (a_dot_b)
      DEALLOCATE (cos_theta)
      DEALLOCATE (normal_1)
      DEALLOCATE (normal_2)
      DEALLOCATE (normal_1_dot_normal_2)
      DEALLOCATE (magnitude_normal_1)
      DEALLOCATE (magnitude_normal_2)
      DEALLOCATE (cos_chi)
      DEALLOCATE (chi)
      DEALLOCATE (chi_deg)
      DEALLOCATE (chirality)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Here we analyze hydration                                                       !
 ! VARIABLE : water                                                                !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! Here is where the analysis of hydration will occur
 ! For each time frame, the water coordinates will need to be read in
 ! Then calculate the presence of hydrogen bonds to water
 ! as well as the presence of bridging waters (these are waters with
 ! hydrogen bonds to more than one residue)

   IF (water.eq.1) THEN      

   ALLOCATE (Owat(3,nwat,nframes))               
   ALLOCATE (H1wat(3,nwat,nframes))
   ALLOCATE (H2wat(3,nwat,nframes))

   DO i=1,3                  
      DO j=1,nwat
         DO t=1,nframes
            Owat(i,j,t) = 0.000000
            H1wat(i,j,t) = 0.000000
            H2wat(i,j,t) = 0.000000
         END DO
      END DO
   END DO

      ALLOCATE (HB_total_at_resi(nres))
      ALLOCATE (HB_both(nres,nwat))
      ALLOCATE (HB_both_sum(nwat))
      ALLOCATE (HB_bridging(nres,nres,nframes))
      ALLOCATE (nbound_t(nframes))
      ALLOCATE (nbridging_t(nframes))
      ALLOCATE (HB_bridging_sum(nres,nres))
      ALLOCATE (HB_bridging_avg(nres,nres))
      ALLOCATE (HB_bridging_NH_CO(two_x_nres,two_x_nres,nframes))
      ALLOCATE (HB_bridging_NH_CO_sum(two_x_nres,two_x_nres))
      ALLOCATE (HB_bridging_NH_CO_avg(two_x_nres,two_x_nres))

      nbound_sum = 0
      nbridging_sum = 0
      nbound_avg = 0.0_dbl
      nbridging_avg = 0.0_dbl

      DO i=0,40
         nbound_dist(i) = 0
         nbridging_dist(i) = 0
         nbound_dist_normalized(i) = 0.0_dbl
         nbridging_dist_normalized(i) = 0.0_dbl
      END DO

      DO i=1,nres
         HB_total_at_resi(i) = 0
      END DO

      DO i=1,nres
         DO j=1,nwat
            HB_both(i,j) = 0
         END DO
      END DO

      DO j=1,nwat
         HB_both_sum(j) = 0
      END DO

      DO t=1,nframes
         DO i=1,nres
            DO j=1,nres
               HB_bridging(i,j,t) = 0
            END DO
         END DO
      END DO

      DO t=1,nframes
         nbound_t(t) = 0
         nbridging_t(t) = 0
      END DO

      DO i=1,nres
         DO j=1,nres
            HB_bridging_sum(i,j) = 0
            HB_bridging_avg(i,j) = 0.0_dbl
         END DO
      END DO

      DO t=1,nframes
         DO i=1,two_x_nres
            DO j=1,two_x_nres
               HB_bridging_NH_CO(i,j,t) = 0
            END DO
         END DO
      END DO

      DO i=1,two_x_nres
         DO j=1,two_x_nres
            HB_bridging_NH_CO_sum(i,j) = 0.0_dbl
            HB_bridging_NH_CO_avg(i,j) = 0.0_dbl
         END DO
      END DO

     ALLOCATE (dN_Owat(nres,nwat))
     ALLOCATE (dO_Owat(nres,nwat))
     ALLOCATE (HB_N(nres,nwat))
     ALLOCATE (HB_CO(nres,nwat))
     ALLOCATE (HB_N_total_at_i(nres,nframes))
     ALLOCATE (HB_CO_total_at_i(nres,nframes))
     ALLOCATE (HB_total_at_i(nres,nframes))
     ALLOCATE (HB_total_at_t(nframes))

     DO i=1,nres
        DO j=1,nwat
           dN_Owat(i,j) = 0.000000
           dO_Owat(i,j) = 0.000000
           HB_N(i,j) = 0
           HB_CO(i,j) = 0
        END DO
     END DO

     DO i=1,nres
        DO t=1,nframes
           HB_N_total_at_i(i,t) = 0
           HB_CO_total_at_i(i,t) = 0
           HB_total_at_i(i,t) = 0
        END DO
     END DO

     DO t=1,nframes
        HB_total_at_t(t) = 0
     END DO

     HBtotal = 0

     HBtotal_normalized = 0.0_dbl
  
     ALLOCATE (HB_total_at_resi_normalized(nres))

     DO i=1,nres
        HB_total_at_resi_normalized(i) = 0.0_dbl
     END DO

 ! Here we read in the water coordinates from the PDB file

      OPEN(6,file='waters.pdb',form='formatted')

        DO t=1,nframes

           DO i=1,5
              READ(6,*)blah
           END DO
                 DO i=1,nwat
                      READ(6,*)atom,b,c,res,d,Owat(1,i,t),Owat(2,i,t),Owat(3,i,t),m,n
                      READ(6,*)atom,b,c,res,d,H1wat(1,i,t),H1wat(2,i,t),H1wat(3,i,t),m,n
                      READ(6,*)atom,b,c,res,d,H2wat(1,i,t),H2wat(2,i,t),H2wat(3,i,t),m,n
                 END DO
           DO i=1,2
              READ(6,*)blah
           END DO

        END DO

    CLOSE(6)

   ! Calculate the relevant quantities for each time frame

 DO t=1,nframes

 ! HB's with backbone NH's

       DO i=1,nres
          DO j=1,nwat
             dN_Owat(i,j) = SQRT(((Nx(i,t)-Owat(1,j,t))**2) + &
                             & ((Ny(i,t)-Owat(2,j,t))**2) + ((Nz(i,t)-Owat(3,j,t))**2))
          END DO
       END DO

       DO i=1,nres
          DO j=1,nwat
             IF (dN_Owat(i,j).lt.3.5000) THEN
                inow = i
                jnow = j

            dOwat_H = SQRT(((Hx(i,t)-Owat(1,j,t))**2) + &
                                       & ((Hy(i,t)-Owat(2,j,t))**2) + ((Hz(i,t)-Owat(3,j,t))**2))

            IF (dOwat_H.lt.2.5000) THEN
               dHN_wat =  SQRT(((Hx(i,t)-Nx(i,t))**2) + &
                                       & ((Hy(i,t)-Ny(i,t))**2) + ((Hz(i,t)-Nz(i,t))**2))

               costheta = ((dN_Owat(i,j)**2) - (dOwat_H**2) - (dHN_wat**2))/(2*(dHN_wat*dOwat_H))
               dH1wat_H = SQRT(((Hx(i,t)-H1wat(1,j,t))**2) + &
                                       & ((Hy(i,t)-H1wat(2,j,t))**2) + ((Hz(i,t)-H1wat(3,j,t))**2))
               dH1wat_N = SQRT(((Nx(i,t)-H1wat(1,j,t))**2) + &
                                       & ((Ny(i,t)-H1wat(2,j,t))**2) + ((Nz(i,t)-H1wat(3,j,t))**2))
               dH2wat_H = SQRT(((Hx(i,t)-H2wat(1,j,t))**2) + &
                                       & ((Hy(i,t)-H2wat(2,j,t))**2) + ((Hz(i,t)-H2wat(3,j,t))**2))
               dH2wat_N = SQRT(((Nx(i,t)-H2wat(1,j,t))**2) + &
                                       & ((Ny(i,t)-H2wat(2,j,t))**2) + ((Nz(i,t)-H2wat(3,j,t))**2))
               E_hb_H1 =  (2.78746712/dN_Owat(i,j)) + (2.78746712/dH1wat_H) &
                                 & - (2.78746712/dOwat_H) - (2.78746712/dH1wat_N)
               E_hb_H2 =  (2.78746712/dN_Owat(i,j)) + (2.78746712/dH2wat_H) &
                                 & - (2.78746712/dOwat_H) - (2.78746712/dH2wat_N)
                  IF (costheta.gt.0.5000) THEN
                     IF (E_hb_H1.lt.-0.05) THEN
                        HB_N(i,j) = 1
                     END IF
                     IF (E_hb_H2.lt.-0.05) THEN
                        HB_N(i,j) = 1
                     END IF
                  END IF
             END IF

             END IF
          END DO
       END DO

 ! prolines (if there is any, make sure that they don't have any HBs with an N-H)

      IF ( nprolines > 0 ) THEN

         DO i=1,nprolines
            j=proline_resnumbers(i)
            DO k=1,nwat
               HB_N(j,k) = 0
            END DO
         END DO

      END IF

       DO i=1,nres
          DO j=1,nwat
             HB_N_total_at_i(i,t) = HB_N_total_at_i(i,t) + HB_N(i,j)
          END DO
       END DO

 ! HB's with backbone CO's

       DO i=1,nres
          DO j=1,nwat
             dO_Owat(i,j) = SQRT(((Ox(i,t)-Owat(1,j,t))**2) + &
                             & ((Oy(i,t)-Owat(2,j,t))**2) + ((Oz(i,t)-Owat(3,j,t))**2))
          END DO
       END DO

       DO i=1,nres
          DO j=1,nwat
             IF (dO_Owat(i,j).lt.3.50000) THEN
                inow = i
                jnow = j

             dO_H1wat = SQRT(((Ox(i,t)-H1wat(1,j,t))**2) + &
                                       & ((Oy(i,t)-H1wat(2,j,t))**2) + ((Oz(i,t)-H1wat(3,j,t))**2))

             dO_H2wat = SQRT(((Ox(i,t)-H2wat(1,j,t))**2) + &
                                       & ((Oy(i,t)-H2wat(2,j,t))**2) + ((Oz(i,t)-H2wat(3,j,t))**2))

        ! If H1 is the H involved in the Hbond

             IF (dO_H1wat.lt.dO_H2wat) THEN
                IF (dO_H1wat.lt.2.50000) THEN
                   dOwat_H1wat = SQRT(((Owat(1,j,t)-H1wat(1,j,t))**2) + &
                                       & ((Owat(2,j,t)-H1wat(2,j,t))**2) + ((Owat(3,j,t)-H1wat(3,j,t))**2))
                   costheta = ((dO_Owat(i,j)**2) - (dO_H1wat**2) - (dOwat_H1wat**2))/(2*(dOwat_H1wat*dO_H1wat))
                   dC_H1wat = SQRT(((Cx(i,t)-H1wat(1,j,t))**2) + &
                                       & ((Cy(i,t)-H1wat(2,j,t))**2) + ((Cz(i,t)-H1wat(3,j,t))**2))
                   dC_Owat = SQRT(((Cx(i,t)-Owat(1,j,t))**2) + &
                                       & ((Cy(i,t)-Owat(2,j,t))**2) + ((Cz(i,t)-Owat(3,j,t))**2))
                   E_hb_wat =  (2.78746712/dO_Owat(i,j)) + (2.78746712/dC_H1wat) &
                                 & - (2.78746712/dO_H1wat) - (2.78746712/dC_Owat)
                   IF (costheta.gt.0.50000) THEN
                      IF (E_hb_wat.lt.-0.05) THEN
                         HB_CO(i,j) = 1
                      END IF
                   END IF
                END IF
              END IF

        ! If H2 is the H involved in the Hbond

             IF (dO_H2wat.lt.dO_H1wat) THEN
                IF (dO_H2wat.lt.2.500000) THEN
                   dOwat_H2wat = SQRT(((Owat(1,j,t)-H2wat(1,j,t))**2) + &
                                       & ((Owat(2,j,t)-H2wat(2,j,t))**2) + ((Owat(3,j,t)-H2wat(3,j,t))**2))
                   costheta = ((dO_Owat(i,j)**2) - (dO_H2wat**2) - (dOwat_H2wat**2))/(2*(dOwat_H2wat*dO_H2wat))
                   dC_H2wat = SQRT(((Cx(i,t)-H2wat(1,j,t))**2) + &
                                       & ((Cy(i,t)-H2wat(2,j,t))**2) + ((Cz(i,t)-H2wat(3,j,t))**2))
                   dC_Owat = SQRT(((Cx(i,t)-Owat(1,j,t))**2) + &
                                       & ((Cy(i,t)-Owat(2,j,t))**2) + ((Cz(i,t)-Owat(3,j,t))**2))
                   E_hb_wat =  (2.78746712/dO_Owat(i,j)) + (2.78746712/dC_H2wat) &
                                 & - (2.78746712/dO_H2wat) - (2.78746712/dC_Owat)
                   IF (costheta.gt.0.50000) THEN
                      IF (E_hb_wat.lt.-0.05) THEN
                         HB_CO(i,j) = 1
                      END IF
                   END IF
                 END IF
              END IF

             END IF
          END DO
       END DO

       DO i=1,nres
          DO j=1,nwat
             HB_CO_total_at_i(i,t) = HB_CO_total_at_i(i,t) + HB_CO(i,j)
          END DO
       END DO


 ! Total HB's at i

                 DO i=1,nres
                    HB_total_at_i(i,t) = HB_N_total_at_i(i,t) + HB_CO_total_at_i(i,t)
                 END DO

     DO i=1,nres
        DO j=1,nwat
           IF (HB_CO(i,j).gt.0) THEN
              HB_both(i,j) = 1
           END IF
           IF (HB_N(i,j).gt.0) THEN
              HB_both(i,j) = 1
           END IF
        END DO
     END DO

     DO j=1,nwat
        DO i=1,nres
           IF (HB_both(i,j).gt.0) THEN
               HB_both_sum(j) = HB_both_sum(j) + 1
           END IF
        END DO
     END DO

     ! Counting the number of bound waters
     ! By this definition, bound waters include bridging waters

     DO j=1,nwat
        IF (HB_both_sum(j).gt.0) THEN
           nbound_t(t) = nbound_t(t) + 1
        END IF
     END DO

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 !  Analysis of bridging waters                                                    !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! Counting the number of bridging waters

     DO j=1,nwat
        IF (HB_both_sum(j).gt.1) THEN
           nbridging_t(t) = nbridging_t(t) + 1
        END IF
     END DO

     ! If bridging, which two residues were bridging?

  3100 FORMAT(' ',I6)

  DO j=1,nwat
     IF (HB_both_sum(j).gt.1) THEN
        OPEN(2202,file='bridged_residues.prn',form='formatted')
            DO i=1,nres
               IF (HB_both(i,j).gt.0) THEN
                  WRITE(2202,3100)i
               END IF
            END DO
         CLOSE(2202)
        OPEN(2202,file='bridged_residues.prn',form='formatted')
        READ(2202,*)one
        READ(2202,*)two
        CLOSE(2202)
        HB_bridging(one,two,t) = 1

        ! Possibility #1 N-H, N-H
        IF (HB_N(one,j).gt.0) THEN
           IF (HB_N(two,j).gt.0) THEN
              HB_bridging_NH_CO(one,two,t) = 1
           END IF
        END IF

        ! Possibility #2 C=O, C=O
        IF (HB_CO(one,j).gt.0) THEN
           IF (HB_CO(two,j).gt.0) THEN
              int_1 = one + 8
              int_2 = two + 8
              HB_bridging_NH_CO(int_1,int_2,t) = 1
           END IF
        END IF

        ! Possibility #3 N-H, C=O
        IF (HB_N(one,j).gt.0) THEN
           IF (HB_CO(two,j).gt.0) THEN
              int_2 = two + 8
              HB_bridging_NH_CO(one,int_2,t) = 1
           END IF
        END IF

        ! Possibility #4 C=O, N-H
        IF (HB_CO(one,j).gt.0) THEN
           IF (HB_N(two,j).gt.0) THEN
              int_1 = one + 8
              HB_bridging_NH_CO(int_1,two,t) = 1
            END IF
        END IF

     END IF
  END DO

 ! Here we are looking for bridging within a residue 
 ! This means that one water has HBs to both the C=O and N-H

  DO j=1,nwat
     DO i=1,nres

     IF (HB_N(i,j).gt.0) THEN
        IF (HB_CO(i,j).gt.0) THEN
        OPEN(2202,file='bridged_residues.prn',form='formatted')
            WRITE(2202,3100)i
        CLOSE(2202)
        OPEN(2202,file='bridged_residues.prn',form='formatted')
            READ(2202,*)one
        CLOSE(2202)

        int_1 = one + 8

        HB_bridging_NH_CO(one,int_1,t) = 1
        HB_bridging_NH_CO(int_1,one,t) = 1

        ! Add bridging within a residue to the count of bridging waters
        nbridging_t(t) = nbridging_t(t) + 1

        END IF
     END IF

     END DO
  END DO

 ! Resetting the necessary variables to zero before the next time frame

       DO i=1,nres
          DO j=1,nwat
             dN_Owat(i,j) = 0.000000
             dO_Owat(i,j) = 0.000000
             HB_N(i,j) = 0
             HB_CO(i,j) = 0
             HB_both(i,j) = 0
          END DO
       END DO

      dOwat_H = 0.000000
      costheta = 0.000000 
      dH1wat_H = 0.000000
      dH1wat_N = 0.000000
      dH2wat_H = 0.000000
      dH2wat_N = 0.000000
      E_hb_H1 = 0.000000
      E_hb_H2 = 0.000000
      dO_H1wat = 0.000000
      dO_H2wat = 0.000000
      dOwat_H1wat = 0.000000 
      dC_H1wat = 0.000000
      dC_Owat = 0.000000
      E_hb_wat = 0.000000
      dOwat_H2wat = 0.000000 
      dC_H2wat = 0.000000

     DO j=1,nwat
        HB_both_sum(j) = 0
     END DO

 END DO

    DEALLOCATE (Owat)
    DEALLOCATE (H1wat)
    DEALLOCATE (H2wat)

    DEALLOCATE (dN_Owat)
    DEALLOCATE (dO_Owat)
    DEALLOCATE (HB_N)
    DEALLOCATE (HB_CO)
    DEALLOCATE (HB_both)
    DEALLOCATE (HB_both_sum)

 ! Calculating and writing out the average hydration properties

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Analysis of Relevant Details about location of HBonds, total # of hbonds etc    !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO t=1,nframes
         DO i=1,nres
            HB_total_at_t(t) = HB_total_at_t(t) +  HB_total_at_i(i,t)
         END DO
      END DO

      DO t=1,nframes
         HBtotal = HBtotal + HB_total_at_t(t)
      END DO  

      DO i=1,nres                
         DO t=1,nframes            
            HB_total_at_resi(i) = HB_total_at_resi(i) +  HB_total_at_i(i,t)
         END DO              
      END DO

 1201 FORMAT(' ',22I6)

      IF (time_evolution.eq.1) THEN

      OPEN(1,file='HB_total_at_i_at_t.prn',form='formatted')
      OPEN(2,file='HB_total_at_t.prn',form='formatted')
      OPEN(3,file='HB_CO_total_at_i_at_t.prn',form='formatted')
      OPEN(4,file='HB_N_total_at_i_at_t.prn',form='formatted')

      DO t=1,nframes             

      WRITE(1,1201)t,( HB_total_at_i(i,t), i=1,nres )
      WRITE(2,1201)t,HB_total_at_t(t)
      WRITE(3,1201)t,( HB_CO_total_at_i(i,t), i=1,nres )
      WRITE(4,1201)t,( HB_N_total_at_i(i,t), i=1,nres )

      END DO

      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)

      END IF

 ! This normalized the total number of HBs to water by the number of frames and
 ! the number of hydrogen bonding groups
 
      HBtotal_normalized = DBLE(HBtotal)/DBLE(nframes*((2*nres)-nprolines))

 ! This normalizes by the number of groups and the number of frames

      DO i=1,nres
         HB_total_at_resi_normalized(i) = DBLE(HB_total_at_resi(i))/DBLE(2*nframes)
      END DO 

      ! normalized by only one group for the prolines 

      IF ( nprolines > 0 ) THEN

         DO i=1,nprolines
            j=proline_resnumbers(i)
            HB_total_at_resi_normalized(j) = DBLE(HB_total_at_resi(j))/DBLE(nframes) 
         END DO

      END IF

      OPEN(1,file='HBtotal.prn',form='formatted')
      OPEN(2,file='HB_total_at_i.prn',form='formatted')

   1211 FORMAT(' ',I10,F10.6)
   1212 FORMAT(' ',2I10,F10.6)

      WRITE(1,1211)HBtotal,HBtotal_normalized

      DO i=1,nres
         WRITE(2,1212)i,HB_total_at_resi(i),HB_total_at_resi_normalized(i)
      END DO

      CLOSE(1)
      CLOSE(2)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Analysis of bound and bridging waters                                           !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO t=1,nframes
         nbound_sum = nbound_sum + nbound_t(t)
         nbridging_sum = nbridging_sum + nbridging_t(t)
      END DO

 ! Average number of bound and bridging waters

      nbound_avg = DBLE(nbound_sum)/DBLE(nframes)
      nbridging_avg = DBLE(nbridging_sum)/DBLE(nframes)

 ! Average of how many bridging waters at each position

      DO t=1,nframes
         DO i=1,nres
            DO j=1,nres
               HB_bridging_sum(i,j) = HB_bridging_sum(i,j) + HB_bridging(i,j,t)
            END DO
         END DO
      END DO

      DO i=1,nres
         DO j=1,nres
            HB_bridging_avg(i,j) = DBLE(HB_bridging_sum(i,j))/DBLE(nframes)
         END DO
      END DO

 ! Writing out the relevant quantities

   IF (time_evolution.eq.1) THEN

   OPEN(1,file='nboundwaters_at_t.prn',form='formatted')

   DO t=1,nframes
      WRITE(1,1007)t,nbound_t(t)
   END DO

   CLOSE(1)

   OPEN(2,file='nbridgingwaters_at_t.prn',form='formatted')

   DO t=1,nframes
      WRITE(2,1007)t,nbridging_t(t)
   END DO

   CLOSE(2)

   END IF

   1209 FORMAT(' ',2F10.6)

   OPEN(1,file='nbound_nbridging_avg.prn',form='formatted')
 
   WRITE(1,1209)nbound_avg,nbridging_avg

   CLOSE(1)

   1010 FORMAT(' ',10F9.6)

   OPEN(1,file='nbridging_i_j.prn',form='formatted')

   WRITE(1,1010)REAL(0),( REAL(j), j=1,nres )
   DO i=1,nres
   WRITE(1,1010)REAL(i),( HB_bridging_avg(i,j), j=1,nres )
   END DO

   CLOSE(1)

   ! Distributions of nbound and nbridging

   DO t=1,nframes
      nbound_dist(nbound_t(t)) = nbound_dist(nbound_t(t)) + 1
      nbridging_dist(nbridging_t(t)) = nbridging_dist(nbridging_t(t)) + 1
   END DO

   DO i=0,40
     nbound_dist_normalized(i) = DBLE(nbound_dist(i))/DBLE(nframes)
     nbridging_dist_normalized(i) = DBLE(nbridging_dist(i))/DBLE(nframes)
   END DO

   OPEN(16,file='nbound_distribution.prn',form='formatted')

   DO i=0,40
      WRITE(16,1211)i,nbound_dist_normalized(i)
   END DO

   CLOSE(16)

   OPEN(17,file='nbridging_distribution.prn',form='formatted')

   DO i=0,40
      WRITE(17,1211)i,nbridging_dist_normalized(i)
   END DO

   CLOSE(17)

 ! Average of how many bridging waters at each position

      DO t=1,nframes
         DO i=1,two_x_nres
            DO j=1,two_x_nres
               HB_bridging_NH_CO_sum(i,j) = HB_bridging_NH_CO_sum(i,j) + HB_bridging_NH_CO(i,j,t)
            END DO
         END DO
      END DO

      DO i=1,two_x_nres
         DO j=1,two_x_nres
            HB_bridging_NH_CO_avg(i,j) = DBLE(HB_bridging_NH_CO_sum(i,j))/DBLE(nframes)
         END DO
      END DO

   1012 FORMAT(' ',16F9.6)

   OPEN(18,file='nbridging_i_j_NH_CO.prn',form='formatted')

   DO i=1,two_x_nres
   WRITE(18,1012)( HB_bridging_NH_CO_avg(i,j), j=1,two_x_nres )
   END DO

   CLOSE(18)

   1210 FORMAT(' ',I10,12F10.6)

   OPEN(1,file='bridging_waters_plot1_NH_NH.prn',form='formatted')

   WRITE(1,1210)0,( REAL(j), j=1,nres )
   DO i=1,nres
      WRITE(1,1210)i,( HB_bridging_NH_CO_avg(i,j), j=1,nres )
   END DO

   CLOSE(1)

   OPEN(2,file='bridging_waters_plot2_CO_CO.prn',form='formatted')

   WRITE(2,1210)0,( REAL(j), j=nres+1,two_x_nres )
   DO i=nres+1,two_x_nres
      WRITE(2,1210)i,( HB_bridging_NH_CO_avg(i,j), j=nres+1,two_x_nres )
   END DO

   CLOSE(2)

   OPEN(3,file='bridging_waters_plot3_NH_CO.prn',form='formatted')

   WRITE(3,1210)0,( REAL(j), j=nres+1,two_x_nres )
   DO i=1,nres
      WRITE(3,1210)i,( HB_bridging_NH_CO_avg(i,j), j=nres+1,two_x_nres )            
   END DO

   CLOSE(3)

   OPEN(4,file='bridging_waters_plot4_CO_NH.prn',form='formatted')

   WRITE(4,1210)0,( REAL(j), j=1,nres )
   DO i=nres+1,two_x_nres           
      WRITE(4,1210)i,( HB_bridging_NH_CO_avg(i,j), j=1,nres )
   END DO

   CLOSE(4)

   OPEN(19,file='nbridging_i_j_NH_CO_plot1.prn',form='formatted')

   WRITE(19,1010)HB_bridging_NH_CO_avg(1,9),HB_bridging_NH_CO_avg(1,10),HB_bridging_NH_CO_avg(1,11),&
                &HB_bridging_NH_CO_avg(1,12),HB_bridging_NH_CO_avg(1,13),HB_bridging_NH_CO_avg(1,14),&
                &HB_bridging_NH_CO_avg(1,15),HB_bridging_NH_CO_avg(1,16)
   WRITE(19,1010)HB_bridging_NH_CO_avg(9,2),HB_bridging_NH_CO_avg(2,10),HB_bridging_NH_CO_avg(2,11),&
                &HB_bridging_NH_CO_avg(2,12),HB_bridging_NH_CO_avg(2,13),HB_bridging_NH_CO_avg(2,14),&
                &HB_bridging_NH_CO_avg(2,15),HB_bridging_NH_CO_avg(2,16)
   WRITE(19,1010)HB_bridging_NH_CO_avg(9,3),HB_bridging_NH_CO_avg(10,3),HB_bridging_NH_CO_avg(3,11),&
                &HB_bridging_NH_CO_avg(3,12),HB_bridging_NH_CO_avg(3,13),HB_bridging_NH_CO_avg(3,14),&
                &HB_bridging_NH_CO_avg(3,15),HB_bridging_NH_CO_avg(3,16)
   WRITE(19,1010)HB_bridging_NH_CO_avg(9,4),HB_bridging_NH_CO_avg(10,4),HB_bridging_NH_CO_avg(11,4),&
                &HB_bridging_NH_CO_avg(4,12),HB_bridging_NH_CO_avg(4,13),HB_bridging_NH_CO_avg(4,14),&
                &HB_bridging_NH_CO_avg(4,15),HB_bridging_NH_CO_avg(4,16)
   WRITE(19,1010)HB_bridging_NH_CO_avg(9,5),HB_bridging_NH_CO_avg(10,5),HB_bridging_NH_CO_avg(11,5),&
                &HB_bridging_NH_CO_avg(12,5),HB_bridging_NH_CO_avg(5,13),HB_bridging_NH_CO_avg(5,14),&
                &HB_bridging_NH_CO_avg(5,15),HB_bridging_NH_CO_avg(5,16)
   WRITE(19,1010)HB_bridging_NH_CO_avg(9,6),HB_bridging_NH_CO_avg(10,6),HB_bridging_NH_CO_avg(11,6),&
                &HB_bridging_NH_CO_avg(12,6),HB_bridging_NH_CO_avg(13,6),HB_bridging_NH_CO_avg(6,14),&
                &HB_bridging_NH_CO_avg(6,15),HB_bridging_NH_CO_avg(6,16)
   WRITE(19,1010)HB_bridging_NH_CO_avg(9,7),HB_bridging_NH_CO_avg(10,7),HB_bridging_NH_CO_avg(11,7),&
                &HB_bridging_NH_CO_avg(12,7),HB_bridging_NH_CO_avg(13,7),HB_bridging_NH_CO_avg(14,7),&
                &HB_bridging_NH_CO_avg(7,15),HB_bridging_NH_CO_avg(7,16)
   WRITE(19,1010)HB_bridging_NH_CO_avg(9,8),HB_bridging_NH_CO_avg(10,8),HB_bridging_NH_CO_avg(11,8),&
                &HB_bridging_NH_CO_avg(12,8),HB_bridging_NH_CO_avg(13,8),HB_bridging_NH_CO_avg(14,8),&
                &HB_bridging_NH_CO_avg(15,8),HB_bridging_NH_CO_avg(8,16)
   CLOSE(19)

   OPEN(20,file='nbridging_i_j_NH_NH_CO_CO_plot2.prn',form='formatted')

   WRITE(20,1010)HB_bridging_NH_CO_avg(1,1),HB_bridging_NH_CO_avg(1,2),HB_bridging_NH_CO_avg(1,3),&
                &HB_bridging_NH_CO_avg(1,4),HB_bridging_NH_CO_avg(1,5),HB_bridging_NH_CO_avg(1,6),&
                &HB_bridging_NH_CO_avg(1,7),HB_bridging_NH_CO_avg(1,8)
   WRITE(20,1010)HB_bridging_NH_CO_avg(9,10),HB_bridging_NH_CO_avg(2,2),HB_bridging_NH_CO_avg(2,3),&
                &HB_bridging_NH_CO_avg(2,4),HB_bridging_NH_CO_avg(2,5),HB_bridging_NH_CO_avg(2,6),&
                &HB_bridging_NH_CO_avg(2,7),HB_bridging_NH_CO_avg(2,8)
    WRITE(20,1010)HB_bridging_NH_CO_avg(9,11),HB_bridging_NH_CO_avg(10,11),HB_bridging_NH_CO_avg(3,3),&
                &HB_bridging_NH_CO_avg(3,4),HB_bridging_NH_CO_avg(3,5),HB_bridging_NH_CO_avg(3,6),&
                &HB_bridging_NH_CO_avg(3,7),HB_bridging_NH_CO_avg(3,8)
    WRITE(20,1010)HB_bridging_NH_CO_avg(9,12),HB_bridging_NH_CO_avg(10,12),HB_bridging_NH_CO_avg(11,12),&
                &HB_bridging_NH_CO_avg(4,4),HB_bridging_NH_CO_avg(4,5),HB_bridging_NH_CO_avg(4,6),&
                &HB_bridging_NH_CO_avg(4,7),HB_bridging_NH_CO_avg(4,8)
    WRITE(20,1010)HB_bridging_NH_CO_avg(9,13),HB_bridging_NH_CO_avg(10,13),HB_bridging_NH_CO_avg(11,13),&
                &HB_bridging_NH_CO_avg(12,13),HB_bridging_NH_CO_avg(5,5),HB_bridging_NH_CO_avg(5,6),&
                &HB_bridging_NH_CO_avg(5,7),HB_bridging_NH_CO_avg(5,8)
    WRITE(20,1010)HB_bridging_NH_CO_avg(9,14),HB_bridging_NH_CO_avg(10,14),HB_bridging_NH_CO_avg(11,14),&
                &HB_bridging_NH_CO_avg(12,14),HB_bridging_NH_CO_avg(13,14),HB_bridging_NH_CO_avg(6,6),&
                &HB_bridging_NH_CO_avg(6,7),HB_bridging_NH_CO_avg(6,8)
    WRITE(20,1010)HB_bridging_NH_CO_avg(9,15),HB_bridging_NH_CO_avg(10,15),HB_bridging_NH_CO_avg(11,15),&
                &HB_bridging_NH_CO_avg(12,15),HB_bridging_NH_CO_avg(13,15),HB_bridging_NH_CO_avg(14,15),&
                &HB_bridging_NH_CO_avg(7,7),HB_bridging_NH_CO_avg(7,8)
    WRITE(20,1010)HB_bridging_NH_CO_avg(9,16),HB_bridging_NH_CO_avg(10,16),HB_bridging_NH_CO_avg(11,16),&
                &HB_bridging_NH_CO_avg(12,16),HB_bridging_NH_CO_avg(13,16),HB_bridging_NH_CO_avg(14,16),&
                &HB_bridging_NH_CO_avg(15,16),HB_bridging_NH_CO_avg(8,8)
   CLOSE(20)

 ! Deallocating the hydration variables after the results have been written

   DEALLOCATE (HB_N_total_at_i)
   DEALLOCATE (HB_CO_total_at_i)
   DEALLOCATE (HB_total_at_i)
   DEALLOCATE (HB_total_at_t)
   DEALLOCATE (HB_total_at_resi)
   DEALLOCATE (HB_bridging)
   DEALLOCATE (nbound_t)
   DEALLOCATE (nbridging_t)
   DEALLOCATE (HB_bridging_sum)
   DEALLOCATE (HB_bridging_avg)
   DEALLOCATE (HB_bridging_NH_CO)
   DEALLOCATE (HB_bridging_NH_CO_sum)
   DEALLOCATE (HB_bridging_NH_CO_avg)

   END IF

      DEALLOCATE (Cx)
      DEALLOCATE (Cy)
      DEALLOCATE (Cz)
      DEALLOCATE (Ox)
      DEALLOCATE (Oy)
      DEALLOCATE (Oz)
      DEALLOCATE (Nx)
      DEALLOCATE (Ny)
      DEALLOCATE (Nz)
      DEALLOCATE (Hx)
      DEALLOCATE (Hy)
      DEALLOCATE (Hz)
      DEALLOCATE (Cax)
      DEALLOCATE (Cay)
      DEALLOCATE (Caz)

      IF ( nprolines > 0 ) THEN

      DEALLOCATE (proline_resnumbers)

      END IF


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Initialize all variables to zero or their starting values                       !       
 ! These are all the variables that keep track of the secondary structure          !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE (turn_1(nres,nframes))
      ALLOCATE (turn_2(nres,nframes))
      ALLOCATE (turn_3(nres,nframes))
      ALLOCATE (turn_4(nres,nframes))
      ALLOCATE (turn_5(nres,nframes))
      ALLOCATE (turn_6(nres,nframes))
      ALLOCATE (turn_7(nres,nframes)) 
      ALLOCATE (turn_1_II(nres,nframes))
      ALLOCATE (turn_2_II(nres,nframes))
      ALLOCATE (turn_3_II(nres,nframes))
      ALLOCATE (turn_4_II(nres,nframes))
      ALLOCATE (turn_5_II(nres,nframes))
      ALLOCATE (turn_6_II(nres,nframes))
      ALLOCATE (turn_7_II(nres,nframes))
      ALLOCATE (helix_3(nres,nframes))
      ALLOCATE (helix_4(nres,nframes))
      ALLOCATE (helix_5(nres,nframes))
      ALLOCATE (parallel_bridge(nres,nframes))
      ALLOCATE (antiparallel_bridge(nres,nframes))
      ALLOCATE (parallel_ladder(nres,nframes))
      ALLOCATE (antiparallel_ladder(nres,nframes))

      ALLOCATE (final_bend(nres,nframes))

      ALLOCATE (final_turn_1(nres,nframes))
      ALLOCATE (final_turn_2(nres,nframes))
      ALLOCATE (final_turn_3(nres,nframes))
      ALLOCATE (final_turn_4(nres,nframes))
      ALLOCATE (final_turn_5(nres,nframes))
      ALLOCATE (final_turn_6(nres,nframes))
      ALLOCATE (final_turn_7(nres,nframes))
      ALLOCATE (final_turn_1_II(nres,nframes))
      ALLOCATE (final_turn_2_II(nres,nframes))
      ALLOCATE (final_turn_3_II(nres,nframes))
      ALLOCATE (final_turn_4_II(nres,nframes))
      ALLOCATE (final_turn_5_II(nres,nframes))
      ALLOCATE (final_turn_6_II(nres,nframes))
      ALLOCATE (final_turn_7_II(nres,nframes))
      ALLOCATE (final_helix_3(nres,nframes))
      ALLOCATE (final_helix_4(nres,nframes))
      ALLOCATE (final_helix_5(nres,nframes))
      ALLOCATE (final_parallel_bridge(nres,nframes))
      ALLOCATE (final_antiparallel_bridge(nres,nframes))
      ALLOCATE (final_parallel_ladder(nres,nframes))
      ALLOCATE (final_antiparallel_ladder(nres,nframes))

      DO t=1,nframes
         DO i=1,nres
             turn_1(i,t) = 0
             turn_2(i,t) = 0
             turn_3(i,t) = 0
             turn_4(i,t) = 0
             turn_5(i,t) = 0
             turn_6(i,t) = 0
             turn_7(i,t) = 0
             turn_1_II(i,t) = 0
             turn_2_II(i,t) = 0
             turn_3_II(i,t) = 0
             turn_4_II(i,t) = 0
             turn_5_II(i,t) = 0
             turn_6_II(i,t) = 0
             turn_7_II(i,t) = 0
             helix_3(i,t) = 0
             helix_4(i,t) = 0
             helix_5(i,t) = 0
             parallel_bridge(i,t) = 0
             antiparallel_bridge(i,t) = 0
             parallel_ladder(i,t) = 0
             antiparallel_ladder(i,t) = 0
             final_bend(i,t) = 0
             final_turn_1(i,t) = 0
             final_turn_2(i,t) = 0
             final_turn_3(i,t) = 0
             final_turn_4(i,t) = 0
             final_turn_5(i,t) = 0
             final_turn_6(i,t) = 0
             final_turn_7(i,t) = 0
             final_turn_1_II(i,t) = 0
             final_turn_2_II(i,t) = 0
             final_turn_3_II(i,t) = 0
             final_turn_4_II(i,t) = 0
             final_turn_5_II(i,t) = 0
             final_turn_6_II(i,t) = 0
             final_turn_7_II(i,t) = 0
             final_helix_3(i,t) = 0
             final_helix_4(i,t) = 0
             final_helix_5(i,t) = 0
             final_parallel_bridge(i,t) = 0
             final_antiparallel_bridge(i,t) = 0
             final_parallel_ladder(i,t) = 0
             final_antiparallel_ladder(i,t) = 0
         END DO
      END DO

      ALLOCATE (total_hbonds_at_t(nframes))
      ALLOCATE (total_bends_at_t(nframes))
      ALLOCATE (total_bends_at_i(nres))

     total_hbonds = 0
     total_bends = 0

     DO t=1,nframes
        total_hbonds_at_t(t) = 0
        total_bends_at_t(t) = 0
     END DO

     DO i=1,nres
        total_bends_at_i(i) = 0
     END DO

      ALLOCATE (total_turn_1_at_t(nframes))
      ALLOCATE (total_turn_1_II_at_t(nframes))
      ALLOCATE (total_turn_2_at_t(nframes))
      ALLOCATE (total_turn_2_II_at_t(nframes))
      ALLOCATE (total_turn_3_at_t(nframes))
      ALLOCATE (total_turn_3_II_at_t(nframes))
      ALLOCATE (total_turn_4_at_t(nframes))
      ALLOCATE (total_turn_4_II_at_t(nframes))
      ALLOCATE (total_turn_5_at_t(nframes))
      ALLOCATE (total_turn_5_II_at_t(nframes))
      ALLOCATE (total_turn_6_at_t(nframes))
      ALLOCATE (total_turn_6_II_at_t(nframes))
      ALLOCATE (total_turn_7_at_t(nframes))
      ALLOCATE (total_turn_7_II_at_t(nframes))
      ALLOCATE (total_helix_3_at_t(nframes))
      ALLOCATE (total_helix_4_at_t(nframes)) 
      ALLOCATE (total_helix_5_at_t(nframes))
      ALLOCATE (total_parallel_bridge_at_t(nframes))
      ALLOCATE (total_antiparallel_bridge_at_t(nframes))
      ALLOCATE (total_parallel_ladder_at_t(nframes))
      ALLOCATE (total_antiparallel_ladder_at_t(nframes))

      DO t=1,nframes
         total_turn_1_at_t(t) = 0
         total_turn_1_II_at_t(t) = 0
         total_turn_2_at_t(t) = 0
         total_turn_2_II_at_t(t) = 0
         total_turn_3_at_t(t) = 0
         total_turn_3_II_at_t(t) = 0
         total_turn_4_at_t(t) = 0
         total_turn_4_II_at_t(t) = 0
         total_turn_5_at_t(t) = 0
         total_turn_5_II_at_t(t) = 0
         total_turn_6_at_t(t) = 0
         total_turn_6_II_at_t(t) = 0
         total_turn_7_at_t(t) = 0
         total_turn_7_II_at_t(t) = 0
         total_helix_3_at_t(t) = 0
         total_helix_4_at_t(t) = 0
         total_helix_5_at_t(t) = 0
         total_parallel_bridge_at_t(t) = 0
         total_antiparallel_bridge_at_t(t) = 0
         total_parallel_ladder_at_t(t) = 0
         total_antiparallel_ladder_at_t(t) = 0
      END DO

      ALLOCATE (total_turn_1_at_i(nres))
      ALLOCATE (total_turn_1_II_at_i(nres))
      ALLOCATE (total_turn_2_at_i(nres))
      ALLOCATE (total_turn_2_II_at_i(nres))
      ALLOCATE (total_turn_3_at_i(nres))
      ALLOCATE (total_turn_3_II_at_i(nres))
      ALLOCATE (total_turn_4_at_i(nres))
      ALLOCATE (total_turn_4_II_at_i(nres))
      ALLOCATE (total_turn_5_at_i(nres))
      ALLOCATE (total_turn_5_II_at_i(nres))
      ALLOCATE (total_turn_6_at_i(nres))
      ALLOCATE (total_turn_6_II_at_i(nres))
      ALLOCATE (total_turn_7_at_i(nres))
      ALLOCATE (total_turn_7_II_at_i(nres))
      ALLOCATE (total_helix_3_at_i(nres))  
      ALLOCATE (total_helix_4_at_i(nres)) 
      ALLOCATE (total_helix_5_at_i(nres))
      ALLOCATE (total_parallel_bridge_at_i(nres))
      ALLOCATE (total_antiparallel_bridge_at_i(nres))
      ALLOCATE (total_parallel_ladder_at_i(nres))
      ALLOCATE (total_antiparallel_ladder_at_i(nres))

      DO i=1,nres
         total_turn_1_at_i(i) = 0
         total_turn_1_II_at_i(i) = 0
         total_turn_2_at_i(i) = 0
         total_turn_2_II_at_i(i) = 0
         total_turn_3_at_i(i) = 0
         total_turn_3_II_at_i(i) = 0
         total_turn_4_at_i(i) = 0
         total_turn_4_II_at_i(i) = 0
         total_turn_5_at_i(i) = 0
         total_turn_5_II_at_i(i) = 0
         total_turn_6_at_i(i) = 0
         total_turn_6_II_at_i(i) = 0
         total_turn_7_at_i(i) = 0
         total_turn_7_II_at_i(i) = 0
         total_helix_3_at_i(i) = 0
         total_helix_4_at_i(i) = 0
         total_helix_5_at_i(i) = 0
         total_parallel_bridge_at_i(i) = 0
         total_antiparallel_bridge_at_i(i) = 0
         total_parallel_ladder_at_i(i) = 0
         total_antiparallel_ladder_at_i(i) = 0
      END DO

      ALLOCATE (PP2(nres,nframes))
      ALLOCATE (PP2_only1s(nres,nframes))
      ALLOCATE (PP2_final(nres,nframes))

      ALLOCATE (total_PP2_at_t(nframes))
      ALLOCATE (total_PP2_at_i(nres))

      DO t=1,nframes
         DO i=1,nres
            PP2(i,t) = 0
            PP2_only1s(i,t) = 0
            PP2_final(i,t) = 0
         END DO
      END DO

      DO t=1,nframes
         total_PP2_at_t(t) = 0
      END DO

      DO i=1,nres
         total_PP2_at_i(i) = 0
      END DO

      total_PP2 = 0

      total_turn_1 = 0
      total_turn_2 = 0
      total_turn_3 = 0
      total_turn_4 = 0
      total_turn_5 = 0
      total_turn_6 = 0
      total_turn_7 = 0
      total_turn_1_II = 0
      total_turn_2_II = 0
      total_turn_3_II = 0
      total_turn_4_II = 0
      total_turn_5_II = 0
      total_turn_6_II = 0
      total_turn_7_II = 0
      total_helix_3 = 0
      total_helix_4 = 0
      total_helix_5 = 0
      total_parallel_bridge = 0
      total_antiparallel_bridge = 0
      total_parallel_ladder = 0
      total_antiparallel_ladder = 0

 ! Calculate whether each residue is in a PPII conformation or not

 ! For the first and last residues, set the phi and psi values to +180

   DO t=1,nframes
      phi(1,t) = 180.0
      psi(1,t) = 180.0
      phi(nres,t) = 180.0
      psi(nres,t) = 180.0
   END DO

 ! Define the PP2 array

 DO t=1,nframes
    DO i=2,nres-1
       IF (phi(i,t).gt.-95.0) THEN
          IF (phi(i,t).lt.-55.0) THEN
             IF (psi(i,t).gt.125.0) THEN
                IF (psi(i,t).lt.165.0) THEN
                   PP2(i,t) = 1
                END IF
            END IF
          END IF
       END IF
    END DO
 END DO

 DO t=1,nframes
    DO i=1,nres
       total_PP2_at_t(t) = total_PP2_at_t(t) + PP2(i,t)
    END DO
 END DO

 DO i=1,nres
    DO t=1,nframes
       total_PP2_at_i(i) = total_PP2_at_i(i) + PP2(i,t)
    END DO
 END DO

 DO i=1,nres
    total_PP2 = total_PP2 + total_PP2_at_i(i)
 END DO

  DO i=1,nres
     total_PP2_at_i_normalized(i) = (DBLE(total_PP2_at_i(i)))/(DBLE(nframes))
  END DO

  total_PP2_normalized = (DBLE(total_PP2))/(DBLE((nres-2)*nframes))

 OPEN(1,file='PP2total.prn',form='formatted')
 OPEN(2,file='PP2total_at_res.prn',form='formatted')

 1004 FORMAT(' ',I15,F15.6)
 1003 FORMAT(' ',I6,I15)

      WRITE(1,1004)total_PP2,total_PP2_normalized

      IF (time_evolution.eq.1) THEN
      OPEN(3,file='PP2total_te.prn',form='formatted')
         DO t=1,nframes
            WRITE(3,1003)t,total_PP2_at_t(t)
         END DO
      CLOSE(3)
      END IF

 1009 FORMAT(' ',2I15,F15.6)

      DO i=1,nres
         WRITE(2,1009)i,total_PP2_at_i(i),total_PP2_at_i_normalized(i)
      END DO

 CLOSE(1)
 CLOSE(2)

 ! PP2_only1s

     DO t=1,nframes
        DO i=1,nres
           IF (PP2(i,t).gt.0) THEN
              PP2_only1s(i,t) = 1
           END IF
        END DO
     END DO

 ! Writing out the integer values of the phi and psi for each residue
 ! Only if the user has specified this as desired output

 IF (phi_psi.eq.1) THEN

    DO i=1,nres
       DO t=1,nframes
          phi_int(i,t) = NINT(phi(i,t))
          psi_int(i,t) = NINT(psi(i,t))
       END DO
    END DO

    OPEN(1,file='phi_of_residues.prn',form='formatted')

    WRITE(1,1201)( i, i=2,nres-1 )
    
    DO t=1,nframes
        WRITE(1,1201)( phi_int(i,t), i=2,nres-1 )
    END DO

    CLOSE(1)

    OPEN(2,file='psi_of_residues.prn',form='formatted')

    WRITE(2,1201)( i, i=2,nres-1 )
    
    DO t=1,nframes
        WRITE(2,1201)( psi_int(i,t), i=2,nres-1 )
    END DO

    CLOSE(2)

 END IF

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Searching for n-turns (2-turns,3-turns,...,7-turns)                             !
 ! Also searching for non-classical turns with HBonds from N-H to C=O              !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1-turns (aka delta turns)

      DO t=1,nframes
         DO i=1,nres-1
            j=i+1
            IF (hb_array(j,i,t).gt.0) THEN
               turn_1(i,t) = 1
               turn_1(j,t) = 1
            END IF
         END DO
      END DO
 
! 2-turns

      DO t=1,nframes
         DO i=1,nres-2
            j=i+2
            IF (hb_array(i,j,t).gt.0) THEN
               turn_2(i,t) = 1
               turn_2(j,t) = 1
            END IF
         END DO
      END DO

! 2-turns with the opposite direction

      DO t=1,nframes
         DO i=3,nres
            j=i-2
            IF (hb_array(i,j,t).gt.0) THEN
               turn_2_II(i,t) = 1
               turn_2_II(j,t) = 1
            END IF
         END DO
      END DO

! 3-turns

      DO t=1,nframes
         DO i=1,nres-3
            j=i+3
            IF (hb_array(i,j,t).gt.0) THEN
               turn_3(i,t) = 1
               turn_3(j,t) = 1
            END IF
         END DO
      END DO

! 3-turns with the opposite direction

      DO t=1,nframes
         DO i=4,nres
            j=i-3
            IF (hb_array(i,j,t).gt.0) THEN
               turn_3_II(i,t) = 1
               turn_3_II(j,t) = 1
            END IF
         END DO
      END DO


! 4-turns

      DO t=1,nframes
         DO i=1,nres-4
            j=i+4
            IF (hb_array(i,j,t).gt.0) THEN
               turn_4(i,t) = 1
               turn_4(j,t) = 1
            END IF
         END DO
      END DO

! 4-turns with the opposite direction

      DO t=1,nframes
         DO i=5,nres
            j=i-4
            IF (hb_array(i,j,t).gt.0) THEN
               turn_4_II(i,t) = 1
               turn_4_II(j,t) = 1
            END IF
         END DO
      END DO


! 5-turns

      DO t=1,nframes
         DO i=1,nres-5
            j=i+5
            IF (hb_array(i,j,t).gt.0) THEN
               turn_5(i,t) = 1
               turn_5(j,t) = 1
            END IF
         END DO
      END DO

! 5-turns with the opposite direction

      DO t=1,nframes
         DO i=6,nres
            j=i-5
            IF (hb_array(i,j,t).gt.0) THEN
               turn_5_II(i,t) = 1
               turn_5_II(j,t) = 1
            END IF
         END DO
      END DO


! 6-turns

      DO t=1,nframes
         DO i=1,nres-6
            j=i+6
            IF (hb_array(i,j,t).gt.0) THEN
               turn_6(i,t) = 1
               turn_6(j,t) = 1
            END IF
         END DO
      END DO

! 6-turns with the opposite direction

      DO t=1,nframes
         DO i=7,nres
            j=i-6
            IF (hb_array(i,j,t).gt.0) THEN
               turn_6_II(i,t) = 1
               turn_6_II(j,t) = 1
            END IF
         END DO
      END DO


! 7-turns

      DO t=1,nframes
         DO i=1,nres-7
            j=i+7
            IF (hb_array(i,j,t).gt.0) THEN
               turn_7(i,t) = 1
               turn_7(j,t) = 1
            END IF
         END DO
      END DO

! 7-turns with the opposite direction

      DO t=1,nframes
         DO i=8,nres
            j=i-7
            IF (hb_array(i,j,t).gt.0) THEN
               turn_7_II(i,t) = 1
               turn_7_II(j,t) = 1
            END IF
         END DO
      END DO

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Searching for n-helices (3-helices,4-helices,5-helices)                         !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 3-helices, aka 3-10 helices

      DO t=1,nframes
         DO i=2,nres-3
            j=i-1
            IF (turn_3(j,t).gt.0) THEN
               IF (turn_3(i,t).gt.0) THEN
                  helix_3(i,t) = 1
                  helix_3(i+1,t) = 1
                  helix_3(i+2,t) = 1
               END IF
            END IF
         END DO
      END DO

! 4-helices, aka alpha helices

      DO t=1,nframes
         DO i=2,nres-4
            j=i-1
            IF (turn_4(j,t).gt.0) THEN
               IF (turn_4(i,t).gt.0) THEN
                  helix_4(i,t) = 1
                  helix_4(i+1,t) = 1
                  helix_4(i+2,t) = 1
                  helix_4(i+3,t) = 1
               END IF
            END IF
         END DO
      END DO

! 5-helices, aka pi helices

      DO t=1,nframes
         DO i=2,nres-5
            j=i-1
            IF (turn_5(j,t).gt.0) THEN
               IF (turn_5(i,t).gt.0) THEN
                  helix_5(i,t) = 1
                  helix_5(i+1,t) = 1
                  helix_5(i+2,t) = 1
                  helix_5(i+3,t) = 1
                  helix_5(i+4,t) = 1
               END IF
            END IF
         END DO
      END DO

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Searching for parallel bridges and antiparallel bridges                         !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Parallel bridges

      DO t=1,nframes
         DO i=2,nres-1
            DO j=2,nres-1
               IF (hb_array(i-1,j,t).gt.0) THEN
                  IF (hb_array(j,i+1,t).gt.0) THEN
                     parallel_bridge(i,t) = 1
                     parallel_bridge(j,t) = 1
                  END IF
               ELSE IF (hb_array(i-1,j+1,t).gt.0) THEN
                  IF (hb_array(j-1,i+1,t).gt.0) THEN
                     parallel_bridge(i,t) = 1
                     parallel_bridge(j,t) = 1
                  END IF
               END IF   
            END DO
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
            DO j=1,nres
               k = i - j
               l = ABS(k)
               IF (l.lt.5) THEN
                  parallel_bridge(i,t) = 0
               END IF
            END DO
         END DO
      END DO

! Antiparallel bridges

      DO t=1,nframes
         DO i=2,nres-1
            DO j=2,nres-1
               IF (hb_array(i,j,t).gt.0) THEN
                  IF (hb_array(j,i,t).gt.0) THEN
                     antiparallel_bridge(i,t) = 1
                     antiparallel_bridge(j,t) = 1
                  END IF
               ELSE IF (hb_array(i-1,j+1,t).gt.0) THEN
                  IF (hb_array(j-1,i+1,t).gt.0) THEN
                     antiparallel_bridge(i,t) = 1
                     antiparallel_bridge(j,t) = 1
                  END IF
               ELSE IF (hb_array(j+1,i-1,t).gt.0) THEN
                  IF (hb_array(i+1,j-1,t).gt.0) THEN
                     antiparallel_bridge(i,t) = 1
                     antiparallel_bridge(j,t) = 1
                  END IF
               END IF 
            END DO
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
            DO j=1,nres
               k = i - j
               l = ABS(k)
               IF (l.lt.5) THEN
                  antiparallel_bridge(i,t) = 0
               END IF
            END DO
         END DO
      END DO


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Searching for parallel ladders and antiparallel ladders                         !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Parallel ladders

      DO t=1,nframes
         DO i=2,nres-1
            DO j=2,nres-1
               IF (parallel_bridge(i,t).gt.0) THEN
                  IF (parallel_bridge(j,t).gt.0) THEN
                     IF (parallel_bridge(i+1,t).gt.0) THEN
                        IF (parallel_bridge(j+1,t).gt.0) THEN
                           parallel_ladder(i,t) = 1
                           parallel_ladder(j,t) = 1
                           parallel_ladder(i+1,t) = 1
                           parallel_ladder(j+1,t) = 1
                        END IF
                     END IF
                  END IF
               END IF
            END DO
         END DO
      END DO

! Antiparallel ladders
    
      DO t=1,nframes
         DO i=2,nres-1
            DO j=2,nres-1
               IF (antiparallel_bridge(i,t).gt.0) THEN
                  IF (antiparallel_bridge(j,t).gt.0) THEN
                     IF (antiparallel_bridge(i+1,t).gt.0) THEN
                        IF (antiparallel_bridge(j-1,t).gt.0) THEN
                           antiparallel_ladder(i,t) = 1
                           antiparallel_ladder(j,t) = 1
                           antiparallel_ladder(i+1,t) = 1
                           antiparallel_ladder(j-1,t) = 1
                        END IF
                     END IF
                  END IF
               END IF
            END DO
         END DO
      END DO

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Hierarchy of secondary structure - assigning only one type to each residue      !
 ! 1. 4-helix (alpha helix)                                                        !
 ! 2. beta ladder                                                                  !
 ! 3. beta bridge                                                                  !
 ! 4. 3-helix (3-10 helix)                                                         !
 ! 5. 5-helix (pi helix)                                                           !
 ! 6. turns (h-bonded turns)                                                       !
 ! 7. bends                                                                        !
 ! **** order switched 2 & 3 from the DSSP paper                                   !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO t=1,nframes
         DO i=1,nres
             final_helix_4(i,t) = helix_4(i,t)
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_helix_3(i,t) = helix_3(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_helix_3(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_helix_5(i,t) = helix_5(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_helix_5(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_helix_5(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_2(i,t) = turn_2(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_2(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_2(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_2(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_3(i,t) = turn_3(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_3(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_3(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_3(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_3(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_4(i,t) = turn_4(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_4(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_4(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_4(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_4(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_4(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_5(i,t) = turn_5(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_5(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_5(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_5(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_5(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_5(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_turn_5(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_parallel_ladder(i,t) = parallel_ladder(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_parallel_ladder(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_parallel_ladder(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_parallel_ladder(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_parallel_ladder(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_parallel_ladder(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_parallel_ladder(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_parallel_ladder(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_antiparallel_ladder(i,t) = antiparallel_ladder(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_antiparallel_ladder(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_antiparallel_ladder(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_antiparallel_ladder(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_antiparallel_ladder(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_antiparallel_ladder(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_antiparallel_ladder(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_antiparallel_ladder(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_antiparallel_ladder(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_parallel_bridge(i,t) = parallel_bridge(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_parallel_bridge(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_parallel_bridge(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_parallel_bridge(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_parallel_bridge(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_parallel_bridge(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_parallel_bridge(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_parallel_bridge(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_parallel_bridge(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_parallel_bridge(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_antiparallel_bridge(i,t) = antiparallel_bridge(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_antiparallel_bridge(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_antiparallel_bridge(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_antiparallel_bridge(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_antiparallel_bridge(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_antiparallel_bridge(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_antiparallel_bridge(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_antiparallel_bridge(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_antiparallel_bridge(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_antiparallel_bridge(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_antiparallel_bridge(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_6(i,t) = turn_6(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  final_turn_6(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_7(i,t) = turn_7(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
               IF (final_turn_6(i,t).gt.0) THEN
                  final_turn_7(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_2_II(i,t) = turn_2_II(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_turn_6(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
               IF (final_turn_7(i,t).gt.0) THEN
                  final_turn_2_II(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_3_II(i,t) = turn_3_II(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_turn_6(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_turn_7(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
               IF (final_turn_2_II(i,t).gt.0) THEN
                  final_turn_3_II(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_4_II(i,t) = turn_4_II(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_turn_6(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_turn_7(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_turn_2_II(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
               IF (final_turn_3_II(i,t).gt.0) THEN
                  final_turn_4_II(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_5_II(i,t) = turn_5_II(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_turn_6(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_turn_7(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_turn_2_II(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_turn_3_II(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
               IF (final_turn_4_II(i,t).gt.0) THEN
                  final_turn_5_II(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_6_II(i,t) = turn_6_II(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_turn_6(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_turn_7(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_turn_2_II(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_turn_3_II(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_turn_4_II(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
               IF (final_turn_5_II(i,t).gt.0) THEN
                  final_turn_6_II(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_7_II(i,t) = turn_7_II(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_6(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_7(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_2_II(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_3_II(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_4_II(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_5_II(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
               IF (final_turn_6_II(i,t).gt.0) THEN
                  final_turn_7_II(i,t) = 0
               END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
             final_turn_1(i,t) = turn_1(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF 
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_6(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_7(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_2_II(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_3_II(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_4_II(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_5_II(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_6_II(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
               IF (final_turn_7_II(i,t).gt.0) THEN
                  final_turn_1(i,t) = 0
               END IF
         END DO
      END DO 

      DO t=1,nframes
         DO i=1,nres
             final_bend(i,t) = bend(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_1(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_6(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_7(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_2_II(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_3_II(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_4_II(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_5_II(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_6_II(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
               IF (final_turn_7_II(i,t).gt.0) THEN
                  final_bend(i,t) = 0
               END IF
         END DO
      END DO

! PP2_only1s(i,t) - hierarchy for PPII

      DO t=1,nframes
         DO i=1,nres
             PP2_final(i,t) = PP2_only1s(i,t)
               IF (final_helix_4(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_antiparallel_ladder(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_parallel_ladder(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_antiparallel_bridge(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_parallel_bridge(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_helix_3(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_helix_5(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_1(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_2(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_3(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_4(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_5(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_6(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_7(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_2_II(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_3_II(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_4_II(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_5_II(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_6_II(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
               IF (final_turn_7_II(i,t).gt.0) THEN
                  PP2_final(i,t) = 0
               END IF
         END DO
      END DO

      IF (time_evolution.eq.1) THEN

      OPEN(1,file='alphahelix_te.prn',form='formatted')
      OPEN(2,file='parallel_beta_ladder_te.prn',form='formatted')
      OPEN(3,file='antiparallel_beta_ladder_te.prn',form='formatted')
      OPEN(4,file='parallel_beta_bridge_te.prn',form='formatted')
      OPEN(5,file='antiparallel_beta_bridge_te.prn',form='formatted')
      OPEN(6,file='3_10_helix_te.prn',form='formatted')
      OPEN(7,file='pi_helix_te.prn',form='formatted')
      OPEN(8,file='2_turns_te.prn',form='formatted')
      OPEN(9,file='3_turns_te.prn',form='formatted')
      OPEN(10,file='4_turns_te.prn',form='formatted')
      OPEN(11,file='5_turns_te.prn',form='formatted')
      OPEN(12,file='6_turns_te.prn',form='formatted')
      OPEN(13,file='7_turns_te.prn',form='formatted')
      OPEN(14,file='2_turns_rev_te.prn',form='formatted')
      OPEN(15,file='3_turns_rev_te.prn',form='formatted')
      OPEN(16,file='4_turns_rev_te.prn',form='formatted')
      OPEN(17,file='5_turns_rev_te.prn',form='formatted')
      OPEN(18,file='6_turns_rev_te.prn',form='formatted')
      OPEN(19,file='7_turns_rev_te.prn',form='formatted')
      OPEN(20,file='bends_te.prn',form='formatted')
      OPEN(21,file='pp2_final_hierarchy_te.prn',form='formatted')

      WRITE(1,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(1,1201)t,( helix_4(i,t), i=1,nres )
      END DO
      CLOSE(1)

      WRITE(2,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(2,1201)t,( parallel_ladder(i,t), i=1,nres )
      END DO
      CLOSE(2)

      WRITE(3,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(3,1201)t,( antiparallel_ladder(i,t), i=1,nres )
      END DO
      CLOSE(3)

      WRITE(4,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(4,1201)t,( parallel_bridge(i,t), i=1,nres )
      END DO
      CLOSE(4)

      WRITE(5,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(5,1201)t,( antiparallel_bridge(i,t), i=1,nres )
      END DO
      CLOSE(5)

      WRITE(6,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(6,1201)t,( helix_3(i,t), i=1,nres )
      END DO
      CLOSE(6)

      WRITE(7,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(7,1201)t,( helix_5(i,t), i=1,nres )
      END DO       
      CLOSE(7)    

      WRITE(8,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(8,1201)t,( turn_2(i,t), i=1,nres )
      END DO
      CLOSE(8)

      WRITE(9,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(9,1201)t,( turn_3(i,t), i=1,nres )
      END DO
      CLOSE(9)

      WRITE(10,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(10,1201)t,( turn_4(i,t), i=1,nres )
      END DO
      CLOSE(10)

      WRITE(11,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(11,1201)t,( turn_5(i,t), i=1,nres )
      END DO
      CLOSE(11)

      WRITE(12,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(12,1201)t,( turn_6(i,t), i=1,nres )
      END DO
      CLOSE(12)

      WRITE(13,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(13,1201)t,( turn_7(i,t), i=1,nres )
      END DO
      CLOSE(13)

      WRITE(14,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(14,1201)t,( turn_2_II(i,t), i=1,nres )
      END DO
      CLOSE(14)

      WRITE(15,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(15,1201)t,( turn_3_II(i,t), i=1,nres )
      END DO
      CLOSE(15)

      WRITE(16,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(16,1201)t,( turn_4_II(i,t), i=1,nres )
      END DO
      CLOSE(16)

      WRITE(17,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(17,1201)t,( turn_5_II(i,t), i=1,nres )
      END DO
      CLOSE(17)

      WRITE(18,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(18,1201)t,( turn_6_II(i,t), i=1,nres )
      END DO
      CLOSE(18)

      WRITE(19,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(19,1201)t,( turn_7_II(i,t), i=1,nres )
      END DO
      CLOSE(19)

      WRITE(20,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(20,1201)t,( bend(i,t), i=1,nres )
      END DO
      CLOSE(20)

      WRITE(21,1201)( i, i=1,nres )
      DO t=1,nframes
         WRITE(21,1201)t,( PP2_final(i,t), i=1,nres )
      END DO
      CLOSE(21)

      END IF

 ! PPII statistics with hierarchy

 ! First resetting the counting variables to zero

  DO t=1,nframes
     total_PP2_at_t(t) = 0
  END DO
 
  DO i=1,nres
     total_PP2_at_i(i) = 0
     total_PP2_at_i_normalized(i) = 0.0_dbl
  END DO

  total_PP2 = 0
  total_PP2_normalized = 0.0_dbl 

  DO t=1,nframes
     DO i=1,nres
        total_PP2_at_t(t) = total_PP2_at_t(t) + PP2_final(i,t)
     END DO
  END DO

  DO i=1,nres
     DO t=1,nframes
        total_PP2_at_i(i) = total_PP2_at_i(i) + PP2_final(i,t)
     END DO
  END DO

  DO i=1,nres
     total_PP2 = total_PP2 + total_PP2_at_i(i)
  END DO

  DO i=1,nres
     total_PP2_at_i_normalized(i) = (DBLE(total_PP2_at_i(i)))/(DBLE(nframes))
  END DO

  total_PP2_normalized = (DBLE(total_PP2))/(DBLE((nres-2)*nframes))

  OPEN(1,file='pp2total_hierarchy.prn',form='formatted')
  OPEN(2,file='pp2total_at_res_hierarchy.prn',form='formatted')

  1114 FORMAT(' ',I6,I15)

       WRITE(1,1004)total_PP2,total_PP2_normalized

       IF (time_evolution.eq.1) THEN
       OPEN(3,file='pp2total_te_hierarchy.prn',form='formatted')
          DO t=1,nframes
             WRITE(3,1114)t,total_PP2_at_t(t)
          END DO
       CLOSE(3)
       END IF

       DO i=1,nres
          WRITE(2,1009)i,total_PP2_at_i(i),total_PP2_at_i_normalized(i)
       END DO

  CLOSE(1)
  CLOSE(2)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 !                             Analysis of the HBonds                              !
 !                                                                                 !
 ! 1. Determine the total number of hbonds occuring during the simulation          !
 ! 2. Determine the time evolution of the number of hbonds during the simulation   !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Part 1: total number of hbonds occuring during the simulation.

      DO t=1,nframes
         DO i=1,nres
            DO j=1,nres
               total_hbonds = total_hbonds + hb_array(i,j,t)
            END DO
         END DO
      END DO

! Part 2: time evolution of hbonds during the simulation.

      DO t=1,nframes
         total_hbonds_at_t(t) = 0
         DO i=1,nres
            DO j=1,nres
               total_hbonds_at_t(t) = total_hbonds_at_t(t) + hb_array(i,j,t)
            END DO
         END DO
      END DO

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 !                             Analysis of the Bends                               !
 !                                                                                 !
 ! 1. Determine the total number of bends occuring during the simulation           !
 ! 2. Determine the time evolution of the number of bends during the simulation    !
 !    (the number of residues with bends for each timeframe)                       !
 ! 3. Determine the total number of bends occuring for each residue (at each       !
 !    position); the total number of bends for each residue position in the repeat !
 !    ; the total number of bends for each type of residue (e.g. G,V etc)          !
 !    is there a trend visible for all glycines, or for all glycines at a certain  !
 !    position for example?                                                        !
 ! 4. Determine the time evolution of each of the quantities in number 3.          !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Part 1: total number of bends during the simulation.

      DO t=1,nframes
         DO i=1,nres
               total_bends = total_bends + bend(i,t)
         END DO
      END DO

! Part 2: time evolution of the number of bends during the simulation.

      DO t=1,nframes
         total_bends_at_t(t) = 0
         DO i=1,nres
               total_bends_at_t(t) = total_bends_at_t(t) + bend(i,t)
         END DO
      END DO

! Part 3: total number of bends occurring for each residue position; for each
!         position in the repeat; and for each type of residue.

      DO t=1,nframes
         DO i=1,nres
               total_bends_at_i(i) = total_bends_at_i(i) + bend(i,t)
         END DO
      END DO
   

! Writing the information about h-bonds and bends to files.

      total_hbonds_normalized = (DBLE(total_hbonds))/(DBLE(nres*nframes))
      total_bends_normalized = (DBLE(total_bends))/(DBLE((nres-4)*nframes))

      DO i=1,nres
         total_bends_at_i_normalized(i) = (DBLE(total_bends_at_i(i)))/(DBLE(nframes))
      END DO

      OPEN(1,file='hbonds_and_bendstotal.prn',form='formatted')
      OPEN(2,file='bends_at_residue_i.prn',form='formatted')

      1008 FORMAT(' ',A15,I15,F10.6)
      1007 FORMAT(' ',2I15)

      WRITE(1,1008)'total_hbonds',total_hbonds,total_hbonds_normalized
      WRITE(1,1008)'total_bends',total_bends,total_bends_normalized 

      DO i=1,nres
         WRITE(2,1009)i,total_bends_at_i(i),total_bends_at_i_normalized(i)
      END DO

      CLOSE(1)
      CLOSE(2) 

      IF (time_evolution.eq.1) THEN

      OPEN(3,file='total_hbonds_te.prn',form='formatted')
      OPEN(4,file='total_bends_te.prn',form='formatted')

      DO t=1,nframes
      WRITE(3,1007)t,total_hbonds_at_t(t)
      WRITE(4,1007)t,total_bends_at_t(t)
      END DO

      CLOSE(3)
      CLOSE(4)

      END IF

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 !             Analysis of other Secondary Structure Elements                      !
 !                                                                                 !
 ! 1. Determine the total number of each occuring during the simulation            !
 ! 2. Determine the time evolution of the number of each during the simulation     !
 ! 3. Determine the total number of each occuring for each residue (at each        !
 !    position); the total number of each for each type of residue (e.g. G,V etc)  !
 !    is there a trend visible for all glycines, or for all glycines at a certain  !
 !    position for example?                                                        !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! Part 1: The total number of each secondary structure element occuring during the
 !         simulation

      DO t=1,nframes
         DO i=1,nres
               total_turn_1 = total_turn_1 + turn_1(i,t)
               total_turn_2 = total_turn_2 + turn_2(i,t)
               total_turn_2_II = total_turn_2_II + turn_2_II(i,t)
               total_turn_3 = total_turn_3 + turn_3(i,t)
               total_turn_3_II = total_turn_3_II + turn_3_II(i,t)
               total_turn_4 = total_turn_4 + turn_4(i,t)
               total_turn_4_II = total_turn_4_II + turn_4_II(i,t)
               total_turn_5 = total_turn_5 + turn_5(i,t)
               total_turn_5_II = total_turn_5_II + turn_5_II(i,t)
               total_turn_6 = total_turn_6 + turn_6(i,t)
               total_turn_6_II = total_turn_6_II + turn_6_II(i,t)
               total_turn_7 = total_turn_7 + turn_7(i,t)
               total_turn_7_II = total_turn_7_II + turn_7_II(i,t)
               total_helix_3 = total_helix_3 + helix_3(i,t)
               total_helix_4 = total_helix_4 + helix_4(i,t)
               total_helix_5 = total_helix_5 + helix_5(i,t)
               total_parallel_bridge = total_parallel_bridge + parallel_bridge(i,t) 
               total_antiparallel_bridge = total_antiparallel_bridge & 
                        & + antiparallel_bridge(i,t)
               total_parallel_ladder = total_parallel_ladder & 
                        & + parallel_ladder(i,t)
               total_antiparallel_ladder = total_antiparallel_ladder &
                        & + antiparallel_ladder(i,t)
         END DO
      END DO

 ! Part 2: The time evolution of each secondary structure element

      DO t=1,nframes
         DO i=1,nres
               total_turn_1_at_t(t) = total_turn_1_at_t(t) + turn_1(i,t)
               total_turn_2_at_t(t) = total_turn_2_at_t(t) + turn_2(i,t)
               total_turn_2_II_at_t(t) = total_turn_2_II_at_t(t) + turn_2_II(i,t)
               total_turn_3_at_t(t) = total_turn_3_at_t(t) + turn_3(i,t)
               total_turn_3_II_at_t(t) = total_turn_3_II_at_t(t) + turn_3_II(i,t)
               total_turn_4_at_t(t) = total_turn_4_at_t(t) + turn_4(i,t)
               total_turn_4_II_at_t(t) = total_turn_4_II_at_t(t) + turn_4_II(i,t)
               total_turn_5_at_t(t) = total_turn_5_at_t(t) + turn_5(i,t)
               total_turn_5_II_at_t(t) = total_turn_5_II_at_t(t) + turn_5_II(i,t)
               total_turn_6_at_t(t) = total_turn_6_at_t(t) + turn_6(i,t)
               total_turn_6_II_at_t(t) = total_turn_6_II_at_t(t) + turn_6_II(i,t)
               total_turn_7_at_t(t) = total_turn_7_at_t(t) + turn_7(i,t)
               total_turn_7_II_at_t(t) = total_turn_7_II_at_t(t) + turn_7_II(i,t)
               total_helix_3_at_t(t) = total_helix_3_at_t(t) + helix_3(i,t)
               total_helix_4_at_t(t) = total_helix_4_at_t(t) + helix_4(i,t)
               total_helix_5_at_t(t) = total_helix_5_at_t(t) + helix_5(i,t)
               total_parallel_bridge_at_t(t) = total_parallel_bridge_at_t(t) &
                       & + parallel_bridge(i,t)
               total_antiparallel_bridge_at_t(t) = total_antiparallel_bridge_at_t(t) & 
                       & + antiparallel_bridge(i,t)
               total_parallel_ladder_at_t(t) = total_parallel_ladder_at_t(t) & 
                       & + parallel_ladder(i,t)
               total_antiparallel_ladder_at_t(t) = total_antiparallel_ladder_at_t(t) &
                       & + antiparallel_ladder(i,t)
         END DO
      END DO

! Part 3: The total number of each secondary structure at each position i

      DO t=1,nframes
         DO i=1,nres
               total_turn_1_at_i(i) = total_turn_1_at_i(i) + turn_1(i,t)
               total_turn_2_at_i(i) = total_turn_2_at_i(i) + turn_2(i,t)
               total_turn_2_II_at_i(i) = total_turn_2_II_at_i(i) + turn_2_II(i,t)
               total_turn_3_at_i(i) = total_turn_3_at_i(i) + turn_3(i,t)
               total_turn_3_II_at_i(i) = total_turn_3_II_at_i(i) + turn_3_II(i,t)
               total_turn_4_at_i(i) = total_turn_4_at_i(i) + turn_4(i,t)
               total_turn_4_II_at_i(i) = total_turn_4_II_at_i(i) + turn_4_II(i,t)
               total_turn_5_at_i(i) = total_turn_5_at_i(i) + turn_5(i,t)
               total_turn_5_II_at_i(i) = total_turn_5_II_at_i(i) + turn_5_II(i,t)
               total_turn_6_at_i(i) = total_turn_6_at_i(i) + turn_6(i,t)
               total_turn_6_II_at_i(i) = total_turn_6_II_at_i(i) + turn_6_II(i,t)
               total_turn_7_at_i(i) = total_turn_7_at_i(i) + turn_7(i,t)
               total_turn_7_II_at_i(i) = total_turn_7_II_at_i(i) + turn_7_II(i,t)
               total_helix_3_at_i(i) = total_helix_3_at_i(i) + helix_3(i,t)
               total_helix_4_at_i(i) = total_helix_4_at_i(i) + helix_4(i,t)
               total_helix_5_at_i(i) = total_helix_5_at_i(i) + helix_5(i,t)
               total_parallel_bridge_at_i(i) = total_parallel_bridge_at_i(i) + &
                                      & parallel_bridge(i,t)
               total_antiparallel_bridge_at_i(i) = total_antiparallel_bridge_at_i(i) + &
                                      & antiparallel_bridge(i,t)
               total_parallel_ladder_at_i(i) = total_parallel_ladder_at_i(i) + &
                                      & parallel_ladder(i,t)
         END DO
      END DO

 ! Writing the information on each secondary structure type

      OPEN(1,file='total_sec_struct_types.prn',form='formatted')
      OPEN(2,file='total_sec_struct_typesII.prn',form='formatted')
      OPEN(3,file='secondary_struc_at_each_res.prn',form='formatted')

      1011 FORMAT(' ',15I7)

      WRITE(1,1011)total_turn_1,total_turn_2,total_turn_3,total_turn_4,&
                    &total_turn_5,total_turn_6,&
                    &total_turn_7,total_helix_3,total_helix_4,total_helix_5,&
                    &total_parallel_bridge,total_antiparallel_bridge,&
                    &total_parallel_ladder,total_antiparallel_ladder

      WRITE(2,1011)total_turn_2_II,total_turn_3_II,total_turn_4_II,&
                    &total_turn_5_II,total_turn_6_II,total_turn_7_II

     1147 FORMAT(' ',A20,15I8)

      WRITE(3,1147)'resnum',( i, i=1,nres )
      WRITE(3,1147)'helix_3',( total_helix_3_at_i(i), i=1,nres )
      WRITE(3,1147)'helix_4',( total_helix_4_at_i(i), i=1,nres )
      WRITE(3,1147)'helix_5',( total_helix_5_at_i(i), i=1,nres )
      WRITE(3,1147)'parallel_bridge',( total_parallel_bridge_at_i(i), i=1,nres )
      WRITE(3,1147)'antiparallel_bridge',( total_antiparallel_bridge_at_i(i), i=1,nres )
      WRITE(3,1147)'parallel_ladder',( total_parallel_ladder_at_i(i), i=1,nres )
      WRITE(3,1147)'antiparallel_ladder',( total_antiparallel_ladder_at_i(i), i=1,nres )

      CLOSE(3)
      CLOSE(1)
      CLOSE(2)

 ! Writing out the normalized values of each secondary structure type

      ALLOCATE (total_helix_3_at_i_normalized(nres))
      ALLOCATE (total_helix_4_at_i_normalized(nres))
      ALLOCATE (total_helix_5_at_i_normalized(nres))     
      ALLOCATE (total_parallel_bridge_at_i_normalized(nres))
      ALLOCATE (total_antiparallel_bridge_at_i_normalized(nres))
      ALLOCATE (total_parallel_ladder_at_i_normalized(nres))
      ALLOCATE (total_antiparallel_ladder_at_i_normalized(nres))

      DO i=1,nres
         total_helix_3_at_i_normalized(i) = 0.0_dbl 
         total_helix_4_at_i_normalized(i) = 0.0_dbl
         total_helix_5_at_i_normalized(i) = 0.0_dbl
         total_parallel_bridge_at_i_normalized(i) = 0.0_dbl
         total_antiparallel_bridge_at_i_normalized(i) = 0.0_dbl
         total_parallel_ladder_at_i_normalized(i) = 0.0_dbl
         total_antiparallel_ladder_at_i_normalized(i) = 0.0_dbl
      END DO

      DO i=1,nres
         total_helix_3_at_i_normalized(i) = (DBLE(total_helix_3_at_i(i)))/&
                                           &(DBLE(nframes))
         total_helix_4_at_i_normalized(i) = (DBLE(total_helix_4_at_i(i)))/&
                                           &(DBLE(nframes))
         total_helix_5_at_i_normalized(i) = (DBLE(total_helix_5_at_i(i)))/&
                                           &(DBLE(nframes))
         total_parallel_bridge_at_i_normalized(i) = (DBLE(total_parallel_bridge_at_i(i)))/&
                                           &(DBLE(nframes))
         total_antiparallel_bridge_at_i_normalized(i) = (DBLE(total_antiparallel_bridge_at_i(i)))/&
                                           &(DBLE(nframes))
         total_parallel_ladder_at_i_normalized(i) = (DBLE(total_parallel_ladder_at_i(i)))/&
                                           &(DBLE(nframes))
         total_antiparallel_ladder_at_i_normalized(i) = (DBLE(total_antiparallel_ladder_at_i(i)))/&
                                           &(DBLE(nframes))
      END DO

      OPEN(4,file='secondary_struc_at_each_res_normalized.prn',form='formatted')

      1148 FORMAT(' ',A20,13F10.6)

      WRITE(4,1148)'resnum',( REAL(i), i=1,nres )
      WRITE(4,1148)'helix_3',( total_helix_3_at_i_normalized(i), i=1,nres )
      WRITE(4,1148)'helix_4',( total_helix_4_at_i_normalized(i), i=1,nres )
      WRITE(4,1148)'helix_5',( total_helix_5_at_i_normalized(i), i=1,nres )
      WRITE(4,1148)'parallel_bridge',( total_parallel_bridge_at_i_normalized(i), i=1,nres )
      WRITE(4,1148)'antiparallel_bridge',( total_antiparallel_bridge_at_i_normalized(i), i=1,nres )
      WRITE(4,1148)'parallel_ladder',( total_parallel_ladder_at_i_normalized(i), i=1,nres )
      WRITE(4,1148)'antiparallel_ladder',( total_antiparallel_ladder_at_i_normalized(i), i=1,nres )

      CLOSE(4)

      DO i=1,nres
         total_helix_3_normalized = total_helix_3_normalized + total_helix_3_at_i_normalized(i)  
         total_helix_4_normalized = total_helix_4_normalized + total_helix_4_at_i_normalized(i)  
         total_helix_5_normalized = total_helix_5_normalized + total_helix_5_at_i_normalized(i)  
         total_parallel_bridge_normalized = total_parallel_bridge_normalized +&
                                            & total_parallel_bridge_at_i_normalized(i)
         total_antiparallel_bridge_normalized = total_antiparallel_bridge_normalized +&
                                            & total_antiparallel_bridge_at_i_normalized(i)
         total_parallel_ladder_normalized = total_parallel_ladder_normalized +&
                                            & total_parallel_ladder_at_i_normalized(i)
         total_antiparallel_ladder_normalized = total_antiparallel_ladder_normalized +&
                                            & total_antiparallel_ladder_at_i_normalized(i)

      END DO

      OPEN(5,file='total_sec_struct_types_normalized.prn',form='formatted')

      WRITE(5,1148)'helix_3',total_helix_3_normalized
      WRITE(5,1148)'helix_4',total_helix_4_normalized
      WRITE(5,1148)'helix_5',total_helix_5_normalized
      WRITE(5,1148)'parallel_bridge',total_parallel_bridge_normalized
      WRITE(5,1148)'antiparallel_bridge',total_antiparallel_bridge_normalized
      WRITE(5,1148)'parallel_ladder',total_parallel_ladder_normalized
      WRITE(5,1148)'antiparallel_ladder',total_antiparallel_ladder_normalized

      CLOSE(5)
      
      DEALLOCATE (total_helix_3_at_i_normalized)
      DEALLOCATE (total_helix_4_at_i_normalized)
      DEALLOCATE (total_helix_5_at_i_normalized)
      DEALLOCATE (total_parallel_bridge_at_i_normalized)
      DEALLOCATE (total_antiparallel_bridge_at_i_normalized)
      DEALLOCATE (total_parallel_ladder_at_i_normalized)
      DEALLOCATE (total_antiparallel_ladder_at_i_normalized)

      IF (time_evolution.eq.1) THEN

      OPEN(1,file='total_turn_2_te.prn',form='formatted')
      OPEN(2,file='total_turn_2_II_te.prn',form='formatted')
      OPEN(3,file='total_turn_3_te.prn',form='formatted')
      OPEN(4,file='total_turn_3_II_te.prn',form='formatted')
      OPEN(5,file='total_turn_4_te.prn',form='formatted')
      OPEN(6,file='total_turn_4_II_te.prn',form='formatted')
      OPEN(7,file='total_turn_5_te.prn',form='formatted')
      OPEN(8,file='total_turn_5_II_te.prn',form='formatted')
      OPEN(9,file='total_turn_6_te.prn',form='formatted')
      OPEN(10,file='total_turn_6_II_te.prn',form='formatted')
      OPEN(11,file='total_turn_7_te.prn',form='formatted')
      OPEN(12,file='total_turn_7_II_te.prn',form='formatted')
      OPEN(13,file='total_helix_3_te.prn',form='formatted')
      OPEN(14,file='total_helix_4_te.prn',form='formatted')
      OPEN(15,file='total_helix_5_te.prn',form='formatted')
      OPEN(16,file='total_parallel_bridge_te.prn',form='formatted')
      OPEN(17,file='total_antiparallel_bridge_te.prn',form='formatted')
      OPEN(18,file='total_parallel_ladder_te.prn',form='formatted')
      OPEN(19,file='total_antiparallel_ladder_te.prn',form='formatted')

      DO t=1,nframes 
      WRITE(1,1007)t,total_turn_2_at_t(t)
      WRITE(2,1007)t,total_turn_2_II_at_t(t)
      WRITE(3,1007)t,total_turn_3_at_t(t)
      WRITE(4,1007)t,total_turn_3_II_at_t(t)
      WRITE(5,1007)t,total_turn_4_at_t(t)
      WRITE(6,1007)t,total_turn_4_II_at_t(t)
      WRITE(7,1007)t,total_turn_5_at_t(t)
      WRITE(8,1007)t,total_turn_5_II_at_t(t)
      WRITE(9,1007)t,total_turn_6_at_t(t)
      WRITE(10,1007)t,total_turn_6_II_at_t(t)
      WRITE(11,1007)t,total_turn_7_at_t(t)
      WRITE(12,1007)t,total_turn_7_II_at_t(t)
      WRITE(13,1007)t,total_helix_3_at_t(t)
      WRITE(14,1007)t,total_helix_4_at_t(t)
      WRITE(15,1007)t,total_helix_5_at_t(t)
      WRITE(16,1007)t,total_parallel_bridge_at_t(t)
      WRITE(17,1007)t,total_antiparallel_bridge_at_t(t)
      WRITE(18,1007)t,total_parallel_ladder_at_t(t)
      WRITE(19,1007)t,total_antiparallel_ladder_at_t(t)
      END DO

      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      CLOSE(5)
      CLOSE(6)
      CLOSE(7)
      CLOSE(8)
      CLOSE(9)
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      CLOSE(19)

      END IF

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Writing out the contact map                                                     !
 ! VARIABLE : contact_map                                                          !                                   
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IF (contact_map.eq.1) THEN

     ALLOCATE (hb_array_total(nres,nres))
     ALLOCATE (hb_array_total_normalized(nres,nres))
     ALLOCATE (total_hbonds_at_res_i(nres))

     DO i=1,nres
        DO j=1,nres
           hb_array_total(i,j) = 0
           hb_array_total_normalized(i,j) = 0.0_dbl
        END DO
     END DO

     DO i=1,nres
        total_hbonds_at_res_i(i) = 0
     END DO

 ! Calculating the hydrogen bonded contact map

       DO i=1,nres
          DO j=1,nres
             DO t=1,nframes
                hb_array_total(i,j) = hb_array_total(i,j) + hb_array(i,j,t)
             END DO
          END DO
       END DO

 ! Calculating the **normalized** hydrogen bonded contact map

       DO i=1,nres
          DO j=1,nres
             hb_array_total_normalized(i,j) = DBLE(hb_array_total(i,j))/DBLE(nframes)
          END DO
       END DO

 ! Writing out the contact map as a list and as an array

      1022 FORMAT(' ',2I10,F10.6)

      OPEN(1,file='contactmap_mathematica.prn',form='formatted')

      DO i=1,nres
         DO j=1,nres
            WRITE(1,1022)i,j,hb_array_total_normalized(i,j)
         END DO
      END DO

      CLOSE(1)

      1023 FORMAT(' ',16F8.5)

      OPEN(2,file='contactmap_mathematica_squares.prn',form='formatted')

      DO j=1,nres
            WRITE(2,1023)( hb_array_total_normalized(i,j), i=1,nres )
      END DO

      CLOSE(2)

     DO i=1,nres
        DO j=1,nres
        total_hbonds_at_res_i(i) = total_hbonds_at_res_i(i) + hb_array_total(i,j)
        total_hbonds_at_res_i(i) = total_hbonds_at_res_i(i) + hb_array_total(j,i)
        END DO
     END DO

     DO i=1,nres
        total_hbonds_at_i_normalized(i) = (DBLE(total_hbonds_at_res_i(i)))/(DBLE(nframes))
     END DO

     OPEN(3,file='hbonds_at_residue_i.prn',form='formatted')
     DO i=1,nres
        WRITE(3,1009)i,total_hbonds_at_res_i(i),total_hbonds_at_i_normalized(i) 
     END DO
     CLOSE(3)

     DEALLOCATE (hb_array_total)
     DEALLOCATE (hb_array_total_normalized)
     DEALLOCATE (total_hbonds_at_res_i)

 END IF

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Writing out the end to end distance distribution                                !
 ! VARIABLE : eed                                                                  !        
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IF (eed.eq.1) THEN

   ALLOCATE (nint_eed_at_t(nframes))

   DO t=1,nframes
      nint_eed_at_t(t) = 0
   END DO

   DO i=1,400
      eed_dist(i) = 0
      eed_normalized(i) = 0.0_dbl
   END DO

   eed_list = 0.0_dbl
   eed_avg = 0.0_dbl
   eed_sum = 0.0_dbl

 ! Calculate the average eed

   DO t=1,nframes
      eed_sum = eed_sum + eed_at_t(t)
   END DO

      eed_avg = DBLE(eed_sum/DBLE(nframes))

   OPEN(1,file='eed_at_t.prn',form='formatted')

      3006 FORMAT(' ',I6,F10.4)

      DO t=1,nframes
         WRITE(1,3006)t,eed_at_t(t)
      END DO

         WRITE(1,3006)nframes,eed_avg

    CLOSE(1)

 ! Calculate the eed distribution

      DO t=1,nframes
         nint_eed_at_t(t) = NINT(eed_at_t(t)*2.0)
      END DO

      DO t=1,nframes
         eed_dist(nint_eed_at_t(t)) = eed_dist(nint_eed_at_t(t)) + 1
      END DO

      DO t=1,400
         eed_normalized(t) = DBLE(eed_dist(t))/DBLE(nframes)
      END DO

 ! Write out the eed distribution

    OPEN(2,file='eed_distribution.prn',form='formatted')

       3007 FORMAT(' ',2F10.5)

      DO t=1,60
         eed_list = DBLE(t)/2.0_dbl
         WRITE(2,3007)eed_list,eed_normalized(t)
      END DO

    CLOSE(2)

    DEALLOCATE (nint_eed_at_t)
    DEALLOCATE (eed_at_t)

 END IF

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! This part of the program writes out where the turns are                         !
 ! For each type of turn, it writes out the count of each type of turn with the    !
 ! first residue of the turn at residue i                                          !
 ! This first will involve recalculating the turn variables to only count the first!
 ! residue                                                                         !
 ! VARIABLE : turn_locations                                                       !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IF (turn_locations.eq.1) THEN

 ! First resetting these variables to zero

      DO t=1,nframes
        DO i=1,nres
            turn_2(i,t) = 0
            turn_3(i,t) = 0
            turn_4(i,t) = 0
            turn_5(i,t) = 0
            turn_6(i,t) = 0
            turn_7(i,t) = 0
            turn_2_II(i,t) = 0
            turn_3_II(i,t) = 0
            turn_4_II(i,t) = 0
            turn_5_II(i,t) = 0
            turn_6_II(i,t) = 0
            turn_7_II(i,t) = 0
         END DO
      END DO
      
      DO i=1,nres
         total_turn_2_at_i(i) = 0
         total_turn_3_at_i(i) = 0
         total_turn_4_at_i(i) = 0
         total_turn_5_at_i(i) = 0
         total_turn_6_at_i(i) = 0
         total_turn_7_at_i(i) = 0
         total_turn_2_II_at_i(i) = 0
         total_turn_3_II_at_i(i) = 0
         total_turn_4_II_at_i(i) = 0
         total_turn_5_II_at_i(i) = 0
         total_turn_6_II_at_i(i) = 0
         total_turn_7_II_at_i(i) = 0
      END DO

 ! Now only counting the first residue in the turn

! 2-turns
      DO t=1,nframes
         DO i=1,nres-2
            j=i+2
            IF (hb_array(i,j,t).gt.0) THEN
               turn_2(i,t) = 1
 !              turn_2(j,t) = 1
            END IF
         END DO
      END DO
! 2-turns with the opposite direction
      DO t=1,nframes
         DO i=3,nres
            j=i-2
            IF (hb_array(i,j,t).gt.0) THEN
 !              turn_2_II(i,t) = 1
               turn_2_II(j,t) = 1
            END IF
         END DO
      END DO
! 3-turns
      DO t=1,nframes
         DO i=1,nres-3
            j=i+3
            IF (hb_array(i,j,t).gt.0) THEN
               turn_3(i,t) = 1
 !              turn_3(j,t) = 1
            END IF
         END DO
      END DO
! 3-turns with the opposite direction
      DO t=1,nframes
         DO i=4,nres
            j=i-3
            IF (hb_array(i,j,t).gt.0) THEN
 !            turn_3_II(i,t) = 1
               turn_3_II(j,t) = 1
            END IF
         END DO
      END DO
! 4-turns
      DO t=1,nframes
         DO i=1,nres-4
            j=i+4
            IF (hb_array(i,j,t).gt.0) THEN
               turn_4(i,t) = 1
 !              turn_4(j,t) = 1
            END IF
         END DO
      END DO
! 4-turns with the opposite direction
      DO t=1,nframes
         DO i=5,nres
            j=i-4
            IF (hb_array(i,j,t).gt.0) THEN
 !              turn_4_II(i,t) = 1
               turn_4_II(j,t) = 1
            END IF
         END DO
      END DO
! 5-turns
      DO t=1,nframes
         DO i=1,nres-5
            j=i+5
            IF (hb_array(i,j,t).gt.0) THEN
               turn_5(i,t) = 1
 !              turn_5(j,t) = 1
            END IF
         END DO
      END DO
! 5-turns with the opposite direction
      DO t=1,nframes
         DO i=6,nres
            j=i-5
            IF (hb_array(i,j,t).gt.0) THEN
 !               turn_5_II(i,t) = 1
               turn_5_II(j,t) = 1
            END IF
         END DO
      END DO
! 6-turns
      DO t=1,nframes
         DO i=1,nres-6
            j=i+6
            IF (hb_array(i,j,t).gt.0) THEN
               turn_6(i,t) = 1
 !              turn_6(j,t) = 1
            END IF
         END DO
      END DO
! 6-turns with the opposite direction
      DO t=1,nframes
         DO i=7,nres
            j=i-6
            IF (hb_array(i,j,t).gt.0) THEN
 !              turn_6_II(i,t) = 1
               turn_6_II(j,t) = 1
            END IF
         END DO
      END DO
! 7-turns
      DO t=1,nframes
         DO i=1,nres-7
            j=i+7
            IF (hb_array(i,j,t).gt.0) THEN
               turn_7(i,t) = 1
 !              turn_7(j,t) = 1
            END IF
         END DO
      END DO
! 7-turns with the opposite direction
      DO t=1,nframes
         DO i=8,nres
            j=i-7
            IF (hb_array(i,j,t).gt.0) THEN
 !              turn_7_II(i,t) = 1
               turn_7_II(j,t) = 1
            END IF
         END DO
      END DO

      DO t=1,nframes
         DO i=1,nres
               total_turn_2_at_i(i) = total_turn_2_at_i(i) + turn_2(i,t)
               total_turn_2_II_at_i(i) = total_turn_2_II_at_i(i) + turn_2_II(i,t)
               total_turn_3_at_i(i) = total_turn_3_at_i(i) + turn_3(i,t)
               total_turn_3_II_at_i(i) = total_turn_3_II_at_i(i) + turn_3_II(i,t)
               total_turn_4_at_i(i) = total_turn_4_at_i(i) + turn_4(i,t)
               total_turn_4_II_at_i(i) = total_turn_4_II_at_i(i) + turn_4_II(i,t)
               total_turn_5_at_i(i) = total_turn_5_at_i(i) + turn_5(i,t)
               total_turn_5_II_at_i(i) = total_turn_5_II_at_i(i) + turn_5_II(i,t)
               total_turn_6_at_i(i) = total_turn_6_at_i(i) + turn_6(i,t)
               total_turn_6_II_at_i(i) = total_turn_6_II_at_i(i) + turn_6_II(i,t)
               total_turn_7_at_i(i) = total_turn_7_at_i(i) + turn_7(i,t)
               total_turn_7_II_at_i(i) = total_turn_7_II_at_i(i) + turn_7_II(i,t)
         END DO
      END DO

      OPEN(1,file='turns_at_each_residue.prn',form='formatted')

      WRITE(1,1147)'resnum',( i, i=1,nres )
      WRITE(1,1147)'turn_2',( total_turn_2_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_3',( total_turn_3_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_4',( total_turn_4_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_5',( total_turn_5_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_6',( total_turn_6_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_7',( total_turn_7_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_2_II',( total_turn_2_II_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_3_II',( total_turn_3_II_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_4_II',( total_turn_4_II_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_5_II',( total_turn_5_II_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_6_II',( total_turn_6_II_at_i(i), i=1,nres )
      WRITE(1,1147)'turn_7_II',( total_turn_7_II_at_i(i), i=1,nres )

      CLOSE(1)

      ALLOCATE (total_turn_2_at_i_normalized(nres))
      ALLOCATE (total_turn_3_at_i_normalized(nres))
      ALLOCATE (total_turn_4_at_i_normalized(nres))
      ALLOCATE (total_turn_5_at_i_normalized(nres))
      ALLOCATE (total_turn_6_at_i_normalized(nres))
      ALLOCATE (total_turn_7_at_i_normalized(nres))
      ALLOCATE (total_turn_2_II_at_i_normalized(nres))
      ALLOCATE (total_turn_3_II_at_i_normalized(nres))
      ALLOCATE (total_turn_4_II_at_i_normalized(nres))
      ALLOCATE (total_turn_5_II_at_i_normalized(nres))
      ALLOCATE (total_turn_6_II_at_i_normalized(nres))
      ALLOCATE (total_turn_7_II_at_i_normalized(nres))

      DO i=1,nres
         total_turn_2_at_i_normalized(i) = 0.0_dbl
         total_turn_3_at_i_normalized(i) = 0.0_dbl
         total_turn_4_at_i_normalized(i) = 0.0_dbl
         total_turn_5_at_i_normalized(i) = 0.0_dbl
         total_turn_6_at_i_normalized(i) = 0.0_dbl
         total_turn_7_at_i_normalized(i) = 0.0_dbl
         total_turn_2_II_at_i_normalized(i) = 0.0_dbl
         total_turn_3_II_at_i_normalized(i) = 0.0_dbl
         total_turn_4_II_at_i_normalized(i) = 0.0_dbl
         total_turn_5_II_at_i_normalized(i) = 0.0_dbl
         total_turn_6_II_at_i_normalized(i) = 0.0_dbl
         total_turn_7_II_at_i_normalized(i) = 0.0_dbl
     END DO

      DO i=1,nres
         total_turn_2_at_i_normalized(i) = (DBLE(total_turn_2_at_i(i)))/(DBLE(nframes))
         total_turn_3_at_i_normalized(i) = (DBLE(total_turn_3_at_i(i)))/(DBLE(nframes))
         total_turn_4_at_i_normalized(i) = (DBLE(total_turn_4_at_i(i)))/(DBLE(nframes))
         total_turn_5_at_i_normalized(i) = (DBLE(total_turn_5_at_i(i)))/(DBLE(nframes))
         total_turn_6_at_i_normalized(i) = (DBLE(total_turn_6_at_i(i)))/(DBLE(nframes))
         total_turn_7_at_i_normalized(i) = (DBLE(total_turn_7_at_i(i)))/(DBLE(nframes))
         total_turn_2_II_at_i_normalized(i) = (DBLE(total_turn_2_II_at_i(i)))/(DBLE(nframes))
         total_turn_3_II_at_i_normalized(i) = (DBLE(total_turn_3_II_at_i(i)))/(DBLE(nframes))
         total_turn_4_II_at_i_normalized(i) = (DBLE(total_turn_4_II_at_i(i)))/(DBLE(nframes))
         total_turn_5_II_at_i_normalized(i) = (DBLE(total_turn_5_II_at_i(i)))/(DBLE(nframes))
         total_turn_6_II_at_i_normalized(i) = (DBLE(total_turn_6_II_at_i(i)))/(DBLE(nframes))
         total_turn_7_II_at_i_normalized(i) = (DBLE(total_turn_7_II_at_i(i)))/(DBLE(nframes))
      END DO

      OPEN(2,file='turns_at_each_residue_normalized.prn',form='formatted')

      WRITE(2,1148)'resnum',( REAL(i), i=1,nres )
      WRITE(2,1148)'turn_2',( total_turn_2_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_3',( total_turn_3_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_4',( total_turn_4_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_5',( total_turn_5_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_6',( total_turn_6_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_7',( total_turn_7_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_2_II',( total_turn_2_II_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_3_II',( total_turn_3_II_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_4_II',( total_turn_4_II_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_5_II',( total_turn_5_II_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_6_II',( total_turn_6_II_at_i_normalized(i), i=1,nres )
      WRITE(2,1148)'turn_7_II',( total_turn_7_II_at_i_normalized(i), i=1,nres )

      CLOSE(2)

      DO i=1,nres
         total_turn_2_normalized = total_turn_2_normalized + total_turn_2_at_i_normalized(i)
         total_turn_3_normalized = total_turn_3_normalized + total_turn_3_at_i_normalized(i)
         total_turn_4_normalized = total_turn_4_normalized + total_turn_4_at_i_normalized(i)
         total_turn_5_normalized = total_turn_5_normalized + total_turn_5_at_i_normalized(i)
         total_turn_6_normalized = total_turn_6_normalized + total_turn_6_at_i_normalized(i)
         total_turn_7_normalized = total_turn_7_normalized + total_turn_7_at_i_normalized(i)
         total_turn_2_II_normalized = total_turn_2_II_normalized + total_turn_2_II_at_i_normalized(i)
         total_turn_3_II_normalized = total_turn_3_II_normalized + total_turn_3_II_at_i_normalized(i)
         total_turn_4_II_normalized = total_turn_4_II_normalized + total_turn_4_II_at_i_normalized(i)
         total_turn_5_II_normalized = total_turn_5_II_normalized + total_turn_5_II_at_i_normalized(i)
         total_turn_6_II_normalized = total_turn_6_II_normalized + total_turn_6_II_at_i_normalized(i)
         total_turn_7_II_normalized = total_turn_7_II_normalized + total_turn_7_II_at_i_normalized(i)
      END DO

      OPEN(5,file='total_sec_struct_types_normalized.prn',form='formatted')

      WRITE(5,1148)'turn_2',((2.0_dbl*total_turn_2_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_3',((2.0_dbl*total_turn_3_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_4',((2.0_dbl*total_turn_4_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_5',((2.0_dbl*total_turn_5_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_6',((2.0_dbl*total_turn_6_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_7',((2.0_dbl*total_turn_7_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_2_II',((2.0_dbl*total_turn_2_II_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_3_II',((2.0_dbl*total_turn_3_II_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_4_II',((2.0_dbl*total_turn_4_II_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_5_II',((2.0_dbl*total_turn_5_II_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_6_II',((2.0_dbl*total_turn_6_II_normalized)/DBLE(nres))
      WRITE(5,1148)'turn_7_II',((2.0_dbl*total_turn_7_II_normalized)/DBLE(nres))
      WRITE(5,1148)'helix_3',(total_helix_3_normalized/DBLE(nres))
      WRITE(5,1148)'helix_4',(total_helix_4_normalized/DBLE(nres))
      WRITE(5,1148)'helix_5',(total_helix_5_normalized/DBLE(nres))
      WRITE(5,1148)'parallel_bridge',(total_parallel_bridge_normalized/DBLE(nres))
      WRITE(5,1148)'antiparallel_bridge',(total_antiparallel_bridge_normalized/DBLE(nres))
      WRITE(5,1148)'parallel_ladder',(total_parallel_ladder_normalized/DBLE(nres))
      WRITE(5,1148)'antiparallel_ladder',(total_antiparallel_ladder_normalized/DBLE(nres))

      CLOSE(5)

      DEALLOCATE (total_turn_2_at_i_normalized)
      DEALLOCATE (total_turn_3_at_i_normalized)
      DEALLOCATE (total_turn_4_at_i_normalized)
      DEALLOCATE (total_turn_5_at_i_normalized)
      DEALLOCATE (total_turn_6_at_i_normalized)
      DEALLOCATE (total_turn_7_at_i_normalized)
      DEALLOCATE (total_turn_2_II_at_i_normalized)
      DEALLOCATE (total_turn_3_II_at_i_normalized)
      DEALLOCATE (total_turn_4_II_at_i_normalized)
      DEALLOCATE (total_turn_5_II_at_i_normalized)
      DEALLOCATE (total_turn_6_II_at_i_normalized)
      DEALLOCATE (total_turn_7_II_at_i_normalized)

 END IF

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !                                                                                 !
 ! Here we write the Ca(i) - Ca(i+3) distance distributions                        !
 ! VARIABLE : inter_Ca_dist                                                        !
 !                                                                                 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IF (inter_Ca_dist.eq.1) THEN

 ! Allocating and setting the necessary arrays to zero

    ALLOCATE (nint_dCa3(nres-2,nframes))
    ALLOCATE (nint_dCa4(nres-3,nframes))
    ALLOCATE (nint_dCa5(nres-4,nframes))
    ALLOCATE (dCa3_dist(nres-2,60))
    ALLOCATE (dCa4_dist(nres-3,60))
    ALLOCATE (dCa5_dist(nres-4,60))
    ALLOCATE (dCa3_normalized(nres-2,60))
    ALLOCATE (dCa4_normalized(nres-3,60))
    ALLOCATE (dCa5_normalized(nres-4,60)) 

    DO t=1,nframes
       DO i=1,nres-2
          nint_dCa3(i,t) = 0
       END DO
       DO i=1,nres-3
          nint_dCa4(i,t) = 0
       END DO
       DO i=1,nres-4
          nint_dCa5(i,t) = 0
       END DO
    END DO

    DO t=1,60
       DO i=1,nres-2
          dCa3_dist(i,t) = 0
          dCa3_normalized(i,t) = 0.00000
       END DO
       DO i=1,nres-3
          dCa4_dist(i,t) = 0
          dCa4_normalized(i,t) = 0.00000
       END DO
       DO i=1,nres-4
          dCa5_dist(i,t) = 0
          dCa5_normalized(i,t) = 0.00000
       END DO
    END DO

 ! Calculate the dCa distribution for each possible contact

 ! First for dCa3 (these are the contacts between i and i+2)

      DO t=1,nframes
         DO i=1,nres-2
            nint_dCa3(i,t) = NINT(dCa3(i,t)*2.00000)
         END DO
      END DO

      DO i=1,nres-2
         DO t=1,nframes
            dCa3_dist(i,nint_dCa3(i,t)) = dCa3_dist(i,nint_dCa3(i,t)) + 1
         END DO
      END DO

      DO t=1,60
         DO i=1,nres-2
            dCa3_normalized(i,t) = DBLE(dCa3_dist(i,t))/DBLE(nframes)
         END DO
      END DO

     3707 FORMAT(' ',13F10.5)

 ! Write out the dCa3 distribution

    OPEN(1,file='dCa3_distribution.prn',form='formatted')

      WRITE(1,3707)1.0,( REAL(i), i=1,nres-2 )  
      DO t=1,60
         eed_list = REAL(t)/2.0
         WRITE(1,3707)eed_list,( dCa3_normalized(i,t), i=1,nres-2 )
      END DO                           

    CLOSE(1)

 ! For dCa4 (these are the contacts between i and i+3)

      DO t=1,nframes
         DO i=1,nres-3
            nint_dCa4(i,t) = NINT(dCa4(i,t)*2.00000)
         END DO
      END DO

      DO i=1,nres-3
         DO t=1,nframes
            dCa4_dist(i,nint_dCa4(i,t)) = dCa4_dist(i,nint_dCa4(i,t)) + 1
         END DO
      END DO

      DO t=1,60
         DO i=1,nres-3
            dCa4_normalized(i,t) = DBLE(dCa4_dist(i,t))/DBLE(nframes)
         END DO
      END DO

 ! Write out the dCa4 distribution

    OPEN(1,file='dCa4_distribution.prn',form='formatted')

      WRITE(1,3707)1.0,( REAL(i), i=1,nres-3 )
      DO t=1,60
         eed_list = REAL(t)/2.0
         WRITE(1,3707)eed_list,( dCa4_normalized(i,t), i=1,nres-3 )
      END DO                           

    CLOSE(1)

 ! For dCa5 (these are the contacts between i and i+4)
      
      DO t=1,nframes
         DO i=1,nres-4
            nint_dCa5(i,t) = NINT(dCa5(i,t)*2.00000)
         END DO
      END DO

      DO i=1,nres-4
         DO t=1,nframes
            dCa5_dist(i,nint_dCa5(i,t)) = dCa5_dist(i,nint_dCa5(i,t)) + 1
         END DO
      END DO

      DO t=1,60
         DO i=1,nres-4
            dCa5_normalized(i,t) = DBLE(dCa5_dist(i,t))/DBLE(nframes)
         END DO
      END DO

 ! Write out the dCa5 distribution

    OPEN(1,file='dCa5_distribution.prn',form='formatted')

      WRITE(1,3707)1.0,( REAL(i), i=1,nres-4 )
      DO t=1,60
         eed_list = REAL(t)/2.0
         WRITE(1,3707)eed_list,( dCa5_normalized(i,t), i=1,nres-4 )
      END DO

    CLOSE(1)

    DEALLOCATE (nint_dCa3)
    DEALLOCATE (nint_dCa4)
    DEALLOCATE (nint_dCa5)
    DEALLOCATE (dCa3_dist)
    DEALLOCATE (dCa4_dist)
    DEALLOCATE (dCa5_dist)
    DEALLOCATE (dCa3_normalized)
    DEALLOCATE (dCa4_normalized)
    DEALLOCATE (dCa5_normalized)

 END IF

      END PROGRAM dssp

 !--function torsion

        REAL FUNCTION torsion (veca, vecb, vecc, vecd)
        IMPLICIT NONE
        REAL :: veca(3),vecb(3),vecc(3),vecd(3)
        REAL :: vecp(3),vecq(3),vecr(3),vecs(3),vecu(3),vecv(3),vecw(3)              
        REAL :: pi
       ! REAL :: torsion
        INTEGER :: i,j,k
        REAL :: dot,sizeu,sizev,dotrw
        pi=3.1415927

           DO i=1,3 
              vecp(i) = veca(i) - vecb(i)
              vecq(i) = vecc(i) - vecb(i)
              vecr(i) = vecb(i) - vecc(i)
              vecs(i) = vecd(i) - vecc(i)
           END DO
 
           DO i=1,3
              j=i+1
              IF (j.eq.4) THEN
                 j=1
              END IF
              k=6-i-j
              vecu(i)=(vecp(j)*vecq(k))-(vecp(k)*vecq(j))
           END DO

           DO i=1,3
              j=i+1
              IF (j.eq.4) THEN
                 j=1
              END IF
              k=6-i-j
              vecv(i)=(vecr(j)*vecs(k))-(vecr(k)*vecs(j))
           END DO

        dot = (vecu(1)*vecv(1))+(vecu(2)*vecv(2))+(vecu(3)*vecv(3))

        sizeu = sqrt((vecu(1)*vecu(1))+(vecu(2)*vecu(2))+(vecu(3)*vecu(3)))
        sizev = sqrt((vecv(1)*vecv(1))+(vecv(2)*vecv(2))+(vecv(3)*vecv(3)))

        torsion = (acos(dot/(sizeu*sizev)))*(180.0000/pi)

           DO i=1,3
              j=i+1
              IF (j.eq.4) THEN
                 j=1
              END IF
              k=6-i-j
              vecw(i)=(vecu(j)*vecv(k))-(vecu(k)*vecv(j))
           END DO

        dotrw = (vecr(1)*vecw(1)) + (vecr(2)*vecw(2)) + (vecr(3)*vecw(3))

        IF (dotrw.gt.0.0) THEN
           torsion = -torsion
        END IF

        RETURN
        END FUNCTION
