MODULE domaindis_mod
!
!     Common block for the microdomain distributions
!
!
USE config_mod
!
SAVE
!
INTEGER, PARAMETER  ::  MD_DOMAIN_CUBE      =  -1
INTEGER, PARAMETER  ::  MD_DOMAIN_CYLINDER  =  -2
INTEGER, PARAMETER  ::  MD_DOMAIN_SPHERE    =  -3
INTEGER, PARAMETER  ::  MD_DOMAIN_FUZZY     =  -4
!
INTEGER                              ::  md_ori_n   = 0
REAL                                 ::  md_sep_fuz = 0.5
INTEGER                              ::  mv_orient  = 0
INTEGER                              ::  mc_num     = 0
INTEGER                              ::  mc_type    = 1
!
END MODULE domaindis_mod
