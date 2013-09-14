!********************************************************************** 
!  EXTERNAL RMC MOVE SUBROUTINE                                         
!                                                                       
!  This subroutine is called by DISCUS during RMC refinements if the    
!  mode is set to external using the command 'set mode,external' in     
!  the RMC sublevel.                                                    
!                                                                       
!  All parameters passed to the subroutine are output values, i.e.      
!  must be set by the subroutine. All information about crystal,        
!  RMC settings, etc. are supplied via MODULES through the        
!  USE statements. Do not change any of those values !!             
!                                                                       
!********************************************************************** 
!  1. Output parameters of 'rmc_genmove_ext'                            
!  -----------------------------------------                            
!                                                                       
!  l   laccept         : TRUE if the generated move is valid,           
!                        FALSE otherwise. Just generating ONE move in   
!                        this routine and use 'laccept' to indicate a   
!                        valid configuration. This will keep 'number of 
!                        generated moves' correct.                      
!  i   natoms          : Number of atoms modified                       
!  r   p_new(3,natoms) : New positions x,y,z of selected atoms          
!  i   i_new(natoms)   : New atom type of selected atoms                
!  i   isel(natoms)    : Number of selected atoms within crystal        
!                                                                       
!********************************************************************** 
!  2.   List of variables defining crystal and RMC settings             
!  --------------------------------------------------------             
!                                                                       
!  2.1. crystal_mod.f90                                                 
!  ----------------                                                     
!                                                                       
!  i   cr_natoms          : Number of atoms in crystal                  
!  i   cr_ncatoms         : Number of atoms per unit cell               
!  i   cr_nscat           : Number of different atom types              
!                           Note: VOID is typ 0 !                       
!  i   cr_icc(3)          : Number of unit cells in x,y,z direction     
!  i   cr_iscat(atoms)    : Gives atom typ for atom 'iatom'             
!                                                                       
!  r   cr_a0(3)           : Contains lattice constants a,b,c            
!  r   cr_win(3)          : Contains angles alfa,beta,gamma             
!  r   cr_gten(3,3)       : Contains metric tensor                      
!  r   cr_rten(3,3)       : Contains reciprocal metric tensor           
!  r   cr_pos(3,iatom)    : Position x,y,z of atom 'iatom' in crystal   
!  r   cr_dim(3,2)        : Crystal dimension (x,y,z | min,max)         
!                                                                       
!  c16 cr_spcgr           : CHARACTER containing space group symbol     
!  c4  cr_at_lis(ityp)    : CHARACTER containing atom name for atom     
!                           type 'ityp'.                                
!                                                                       
!  2.2. rmc_mod.f90                                                     
!  ------------                                                         
!                                                                       
!  r   rmc_maxmove(3,ityp)   : Value sx,sy,sz for atom type 'ityp' set  
!                              with 'set move,..' command.              
!  r   rmc_mindist(ityp,jtyp): Minimal allowed distance between atom    
!                              types 'ityp' and 'jtyp'.                 
!  l   rmc_allowed(ityp)     : TRUE is atom type 'ityp' is selected for 
!                              RMC refinement ('sele' command)          
!                                                                       
!  2.3. chem_mod.f90                                                    
!  -------------                                                        
!                                                                       
!  r   chem_ave_pos(3,isite) : Average position x,y,z of site 'isite'   
!  r   chem_ave_sig(3,isite) : Standard deviation sx,sy,sz of 'isite'   
!                                                                       
!********************************************************************** 
!  3. Useful subroutines supplied by DISCUS                             
!  ----------------------------------------                             
!                                                                       
!  real function ran1(idum):                                            
!  -------------------------                                            
!    Gives random number [0,1]. The parameter 'idum' is defined in      
!    MODULE random_mod. The random generator is initialised before RMC level 
!    is entered.                                                        
!                                                                       
!  real function gasdev(sigma):                                         
!  ----------------------------                                         
!    Gives Gaussian distributed random number. FWHM of the distribution 
!    is given by the parameter 'sigma'.                                 
!                                                                       
!  real function do_bang(lspace,u(3),v(3),w(3)):                        
!  ---------------------------------------------                        
!    Gives the angle between vector v-u and v-w in space defined by     
!    lspace (.true. = real space , .false. = reciprocal space)          
!                                                                       
!  real function do_blen(lspace,u(3),v(3)):                             
!  ----------------------------------------                             
!    Calculates the length of the vector v-u in the space defined by    
!    lspace (.true. = real space , .false. = reciprocal space)          
!                                                                       
!  subroutine indextocell (iatom,icell(3),isite)                        
!  ---------------------------------------------                        
!    Gives number of unit cell 'icell(3)' and atom site number 'isite'  
!    for the atom index 'iatom'. This routine works only properly, if   
!    the supplied structure is given in the order specified in the      
!    the DISCUS manual.                                                 
!                                                                       
!  subroutine celltoindex (icell(3),isite,iatom)                        
!  ---------------------------------------------                        
!    Gives the atom index for unit cell 'icell(3)' and site 'isite'.    
!    Again the structure must be given in the 'DISCUS order' to be able 
!    to use this routine.                                               
!                                                                       
!********************************************************************** 
!********************************************************************** 
!                                                                       
      SUBROUTINE rmc_genmove_ext (laccept, natoms, p_new, i_new, isel) 
!                                                                       
!------ We recommend 'implicit none' to avoid the accidental use        
!------ of DISCUS variables defined in the include files                
!                                                                       
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE rmc_mod 
      USE errlist_mod 
      USE random_mod
      IMPLICIT none 
!                                                                       
!------ Include files containing array sizes, crystal definitions,      
!------ error and RMC variables                                         
!                                                                       
       
!                                                                       
!------ Declaration of the variables passed down from DISCUS.           
!------ All information must be passed on via these variables, do       
!------ NOT change any other DISCUS variables !                         
!                                                                       
      REAL p_new (3, rmc_max_atom) 
      INTEGER i_new (rmc_max_atom) 
      INTEGER isel (rmc_max_atom) 
      INTEGER natoms 
      LOGICAL laccept 
!                                                                       
!********************************************************************** 
!  Here is as an example the part of the code for shifting a single     
!  atom. This is exactly what mode SHIFT does. Remove the code, or      
!  modify it to your needs. For more examples, have a look at the       
!  subroutine 'rmc_genmove' in the source file 'rmc.f'.                 
!********************************************************************** 
!                                                                       
!------ declaration of local variables, functions (-> implicit none)    
!                                                                       
      INTEGER i 
      REAL ran1, gasdev 
!                                                                       
!------ we will modify just ONE atom                                    
!                                                                       
      natoms = 1 
!                                                                       
!------ we select one atom at random,                                   
!------ 'idum' is defined in random_mod
!                                                                       
      isel (1) = int (ran1 (idum) * cr_natoms) + 1 
!                                                                       
!------ set new atom type number, here is is unchanged                  
!                                                                       
      i_new (1) = cr_iscat (isel (1) ) 
!                                                                       
!------ Check if the selected atom is on the list of                    
!------ selected atoms set by the command 'sele'                        
!                                                                       
      laccept = rmc_allowed (i_new (1) ) 
!                                                                       
!------ if we have a valid atom, move it. The real function             
!------ 'gasdev' generates Gaussian distributed random number.          
!------ 'rmc_maxmove' contains the FWHM of the distribution             
!------ set by the command 'set move,...'                               
!                                                                       
      IF (laccept) then 
         DO i = 1, 3 
         p_new (i, 1) = cr_pos (i, isel (1) ) + gasdev (rmc_maxmove (i, &
         i_new (1) ) )                                                  
         ENDDO 
      ENDIF 
!                                                                       
!------ We don't have to worry here about moving atoms too close        
!------ together, this is checked in the main RMC routine.              
!                                                                       
      END SUBROUTINE rmc_genmove_ext                
