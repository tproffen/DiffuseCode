MODULE spcgr_mod
!+
!     Variables needed for the spacegroups
!-
PUBLIC
PRIVATE i
SAVE
!
INTEGER, PARAMETER  ::  SPCGR_MAX  =  314
!
CHARACTER(LEN=16), DIMENSION(SPCGR_MAX)     ::  spcgr_name      ! Normal space group name
CHARACTER(LEN=16), DIMENSION(16:74    ,1:6) ::  spcgr_name_set  ! alternative setting in orthorhombic only
INTEGER          , DIMENSION(SPCGR_MAX,2)   ::  spcgr_num       ! Space group number
INTEGER          , DIMENSION(SPCGR_MAX)     ::  spcgr_syst      ! Crystal system to which spce group belongs to
INTEGER                                     :: i
!
!                                   abc        bac        cab        cba        bca        acb
DATA (spcgr_name_set(16,i),i=1,6) /'P222   ', 'P222   ', 'P222   ', 'P222   ', 'P222   ', 'P222   '/
DATA (spcgr_name_set(17,i),i=1,6) /'P2221  ', 'P2221  ', 'P2122  ', 'P2122  ', 'P2212  ', 'P2212  '/
DATA (spcgr_name_set(18,i),i=1,6) /'P21212 ', 'P21212 ', 'P22121 ', 'P22121 ', 'P21221 ', 'P21221 '/
DATA (spcgr_name_set(19,i),i=1,6) /'P212121', 'P212121', 'P212121', 'P212121', 'P212121', 'P212121'/
DATA (spcgr_name_set(20,i),i=1,6) /'C2221  ', 'C2221  ', 'A2122  ', 'A2122  ', 'B2212  ', 'B2212  '/
DATA (spcgr_name_set(21,i),i=1,6) /'C222   ', 'C222   ', 'A222   ', 'A222   ', 'B222   ', 'B222   '/
DATA (spcgr_name_set(22,i),i=1,6) /'F222   ', 'F222   ', 'F222   ', 'F222   ', 'F222   ', 'F222   '/
DATA (spcgr_name_set(23,i),i=1,6) /'I222   ', 'I222   ', 'I222   ', 'I222   ', 'I222   ', 'I222   '/
DATA (spcgr_name_set(24,i),i=1,6) /'I212121', 'I212121', 'I212121', 'I212121', 'I212121', 'I212121'/
!
!                                   abc        bac        cab        cba        bca        acb
DATA (spcgr_name_set(25,i),i=1,6) /'Pmm2   ', 'Pmm2   ', 'P2mm   ', 'P2mm   ', 'Pm2m   ', 'Pm2m   '/
DATA (spcgr_name_set(26,i),i=1,6) /'Pmc21  ', 'Pcm21  ', 'P21ma  ', 'P21am  ', 'Pb21m  ', 'Pm21b  '/
DATA (spcgr_name_set(27,i),i=1,6) /'Pcc2   ', 'Pcc2   ', 'P2aa   ', 'P2aa   ', 'Pb2b   ', 'Pb2b   '/
DATA (spcgr_name_set(28,i),i=1,6) /'Pma2   ', 'Pbm2   ', 'P2mb   ', 'P2cm   ', 'Pc2m   ', 'Pm2a   '/
DATA (spcgr_name_set(29,i),i=1,6) /'Pca21  ', 'Pbc21  ', 'P21ab  ', 'P21ca  ', 'Pc21b  ', 'Pb21a  '/
DATA (spcgr_name_set(30,i),i=1,6) /'Pnc2   ', 'Pcn2   ', 'P2na   ', 'P2an   ', 'Pb2n   ', 'Pn2b   '/
DATA (spcgr_name_set(31,i),i=1,6) /'Pmn21  ', 'Pnm21  ', 'P21mn  ', 'P21nm  ', 'Pn21m  ', 'Pm21n  '/
DATA (spcgr_name_set(32,i),i=1,6) /'Pba2   ', 'Pba2   ', 'P2cb   ', 'P2cb   ', 'Pc2a   ', 'Pc2a   '/
DATA (spcgr_name_set(33,i),i=1,6) /'Pna21  ', 'Pbn21  ', 'P21nb  ', 'P21cn  ', 'Pc21n  ', 'Pn21a  '/
DATA (spcgr_name_set(34,i),i=1,6) /'Pnn2   ', 'Pnn2   ', 'P2nn   ', 'P2nn   ', 'Pn2n   ', 'Pn2n   '/
DATA (spcgr_name_set(35,i),i=1,6) /'Cmm2   ', 'Cmm2   ', 'A2mm   ', 'A2mm   ', 'Bm2m   ', 'Bm2m   '/
DATA (spcgr_name_set(36,i),i=1,6) /'Cmc21  ', 'Ccm21  ', 'A21ma  ', 'A21am  ', 'Bb21m  ', 'Bm21b  '/
DATA (spcgr_name_set(37,i),i=1,6) /'Ccc2   ', 'Ccc2   ', 'A2aa   ', 'A2aa   ', 'Bb2b   ', 'Bb2b   '/
DATA (spcgr_name_set(38,i),i=1,6) /'Amm2   ', 'Bmm2   ', 'B2mm   ', 'C2mm   ', 'Cm2m   ', 'Am2m   '/
DATA (spcgr_name_set(39,i),i=1,6) /'Aem2   ', 'Bme2   ', 'B2em   ', 'C2me   ', 'Cm2e   ', 'Ae2m   '/
DATA (spcgr_name_set(40,i),i=1,6) /'Ama2   ', 'Bbm2   ', 'B2mb   ', 'C2cm   ', 'Cc2m   ', 'Am2a   '/
DATA (spcgr_name_set(41,i),i=1,6) /'Aea2   ', 'Bbe2   ', 'B2eb   ', 'C2ce   ', 'Cc2e   ', 'Ae2a   '/
DATA (spcgr_name_set(42,i),i=1,6) /'Fmm2   ', 'Fmm2   ', 'F2mm   ', 'F2mm   ', 'Fm2m   ', 'Fm2m   '/
DATA (spcgr_name_set(43,i),i=1,6) /'Fdd2   ', 'Fdd2   ', 'F2dd   ', 'F2dd   ', 'Fd2d   ', 'Fd2d   '/
DATA (spcgr_name_set(44,i),i=1,6) /'Imm2   ', 'Imm2   ', 'I2mm   ', 'I2mm   ', 'Im2m   ', 'Im2m   '/
DATA (spcgr_name_set(45,i),i=1,6) /'Iba2   ', 'Iba2   ', 'I2cb   ', 'I2cb   ', 'Ic2a   ', 'Ic2a   '/
DATA (spcgr_name_set(46,i),i=1,6) /'Ima2   ', 'Ibm2   ', 'I2mb   ', 'I2cm   ', 'Ic2m   ', 'Im2a   '/
!
!                                   abc        bac        cab        cba        bca        acb
DATA (spcgr_name_set(47,i),i=1,6) /'Pmmm   ', 'Pmmm   ', 'Pmmm   ', 'Pmmm   ', 'Pmmm   ', 'Pmmm   '/
DATA (spcgr_name_set(48,i),i=1,6) /'Pnnn   ', 'Pnnn   ', 'Pnnn   ', 'Pnnn   ', 'Pnnn   ', 'Pnnn   '/
DATA (spcgr_name_set(49,i),i=1,6) /'Pccm   ', 'Pccm   ', 'Pmaa   ', 'Pmaa   ', 'Pbmb   ', 'Pbmb   '/
DATA (spcgr_name_set(50,i),i=1,6) /'Pban   ', 'Pban   ', 'Pncb   ', 'Pncb   ', 'Pcna   ', 'Pcna   '/
DATA (spcgr_name_set(51,i),i=1,6) /'Pmma   ', 'Pmmb   ', 'Pbmm   ', 'Pcmm   ', 'Pmcm   ', 'Pmam   '/
DATA (spcgr_name_set(52,i),i=1,6) /'Pnna   ', 'Pnnb   ', 'Pbnn   ', 'Pcnn   ', 'Pncn   ', 'Pnan   '/
DATA (spcgr_name_set(53,i),i=1,6) /'Pmna   ', 'Pnmb   ', 'Pbmn   ', 'Pcnm   ', 'Pncm   ', 'Pman   '/
DATA (spcgr_name_set(54,i),i=1,6) /'Pcca   ', 'Pccb   ', 'Pbaa   ', 'Pcaa   ', 'Pbcb   ', 'Pbab   '/
DATA (spcgr_name_set(55,i),i=1,6) /'Pbam   ', 'Pbam   ', 'Pmcb   ', 'Pmcb   ', 'Pcma   ', 'Pcma   '/
DATA (spcgr_name_set(56,i),i=1,6) /'Pccn   ', 'Pccn   ', 'Pnaa   ', 'Pnaa   ', 'Pbnb   ', 'Pbnb   '/
DATA (spcgr_name_set(57,i),i=1,6) /'Pbcm   ', 'Pcam   ', 'Pmca   ', 'Pmab   ', 'Pbma   ', 'Pcmb   '/
DATA (spcgr_name_set(58,i),i=1,6) /'Pnnm   ', 'Pnnm   ', 'Pmnn   ', 'Pmnn   ', 'Pnmn   ', 'Pnmn   '/
DATA (spcgr_name_set(59,i),i=1,6) /'Pmmn   ', 'Pmmn   ', 'Pnmm   ', 'Pnmm   ', 'Pmnm   ', 'Pmnm   '/
DATA (spcgr_name_set(60,i),i=1,6) /'Pbcn   ', 'Pcan   ', 'Pnca   ', 'Pnab   ', 'Pbna   ', 'Pcnb   '/
DATA (spcgr_name_set(61,i),i=1,6) /'Pbca   ', 'Pcab   ', 'Pbca   ', 'Pcab   ', 'Pbca   ', 'Pcab   '/
DATA (spcgr_name_set(62,i),i=1,6) /'Pnma   ', 'Pmnb   ', 'Pbnm   ', 'Pcmn   ', 'Pmcn   ', 'Pnam   '/
DATA (spcgr_name_set(63,i),i=1,6) /'Cmcm   ', 'Ccmm   ', 'Amma   ', 'Amam   ', 'Bbmm   ', 'Bmmb   '/
DATA (spcgr_name_set(64,i),i=1,6) /'Cmce   ', 'Ccme   ', 'Aema   ', 'Aeam   ', 'Bbem   ', 'Bmeb   '/
DATA (spcgr_name_set(65,i),i=1,6) /'Cmmm   ', 'Cmmm   ', 'Ammm   ', 'Ammm   ', 'Bmmm   ', 'Bmmm   '/
DATA (spcgr_name_set(66,i),i=1,6) /'Cccm   ', 'Cccm   ', 'Amaa   ', 'Amaa   ', 'Bbmb   ', 'Bbmb   '/
DATA (spcgr_name_set(67,i),i=1,6) /'Cmme   ', 'Cmmb   ', 'Aemm   ', 'Acmm   ', 'Bmem   ', 'Bmam   '/
DATA (spcgr_name_set(68,i),i=1,6) /'Ccce   ', 'Cccb   ', 'Aeaa   ', 'Acaa   ', 'Bbeb   ', 'Bbab   '/
DATA (spcgr_name_set(69,i),i=1,6) /'Fmmm   ', 'Fmmm   ', 'Fmmm   ', 'Fmmm   ', 'Fmmm   ', 'Fmmm   '/
DATA (spcgr_name_set(70,i),i=1,6) /'Fddd   ', 'Fddd   ', 'Fddd   ', 'Fddd   ', 'Fddd   ', 'Fddd   '/
DATA (spcgr_name_set(71,i),i=1,6) /'Immm   ', 'Immm   ', 'Immm   ', 'Immm   ', 'Immm   ', 'Immm   '/
DATA (spcgr_name_set(72,i),i=1,6) /'Ibam   ', 'Ibam   ', 'Imcb   ', 'Imcb   ', 'Icma   ', 'Icma   '/
DATA (spcgr_name_set(73,i),i=1,6) /'Ibca   ', 'Icab   ', 'Ibca   ', 'Icab   ', 'Ibca   ', 'Icab   '/
DATA (spcgr_name_set(74,i),i=1,6) /'Imma   ', 'Immb   ', 'Ibmm   ', 'Icmm   ', 'Imcm   ', 'Imam   '/
!
END MODULE spcgr_mod
