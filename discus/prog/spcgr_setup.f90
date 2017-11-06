MODULE spcgr_setup_mod
!
CONTAINS
!********************************************************************** 
      SUBROUTINE spcgr_setup 
!-                                                                      
!     Sets up the space group symbol tables.                            
!+                                                                      
      USE spcgr_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER TRICLINIC 
      INTEGER MONOCLINIC_B 
      INTEGER MONOCLINIC_C 
      INTEGER ORTHORHOMBIC 
      INTEGER TETRAGONAL 
      INTEGER TRIGONAL 
      INTEGER RHOMBOHEDRAL 
      INTEGER HEXAGONAL 
      INTEGER CUBIC 
!                                                                       
      PARAMETER (TRICLINIC = 1) 
      PARAMETER (MONOCLINIC_B = 2) 
      PARAMETER (MONOCLINIC_C = 3) 
      PARAMETER (ORTHORHOMBIC = 4) 
      PARAMETER (TETRAGONAL = 5) 
      PARAMETER (TRIGONAL = 6) 
      PARAMETER (RHOMBOHEDRAL = 7) 
      PARAMETER (HEXAGONAL = 8) 
      PARAMETER (CUBIC = 9) 
!                                                                       
      INTEGER i 
!                                                                       
      DO i = 1, SPCGR_MAX 
      spcgr_num (i, 1) = i 
      spcgr_num (i, 2) = 0 
      ENDDO 
!                                                                       
!     The triclinic space groups                                        
!                                                                       
      spcgr_name (1)  = 'P1              ' 
      spcgr_syst (1) = TRICLINIC 
                                                                        
      spcgr_name (2)  = 'P-1             ' 
      spcgr_syst (2) = TRICLINIC 
!                                                                       
!     The monoclinic space groups                                       
!                                                                       
      DO i = 3, 15 
      spcgr_syst (i) = MONOCLINIC_B 
      ENDDO 
                                                                        
      spcgr_name (3)  = 'P2              ' 
      spcgr_name (4)  = 'P21             ' 
      spcgr_name (5)  = 'C2              ' 
      spcgr_name (6)  = 'Pm              ' 
      spcgr_name (7)  = 'Pc              ' 
      spcgr_name (8)  = 'Cm              ' 
      spcgr_name (9)  = 'Cc              ' 
      spcgr_name (10)  = 'P2/m            ' 
      spcgr_name (11)  = 'P21/m           ' 
      spcgr_name (12)  = 'C2/m            ' 
      spcgr_name (13)  = 'P2/c            ' 
      spcgr_name (14)  = 'P21/c           ' 
      spcgr_name (15)  = 'C2/c            ' 
                                                                        
!                                                                       
!     The orthorhombic space groups                                     
!                                                                       
      DO i = 16, 74 
      spcgr_syst (i) = ORTHORHOMBIC 
      ENDDO 
                                                                        
      spcgr_name (16)  = 'P222            ' 
      spcgr_name (17)  = 'P2221           ' 
      spcgr_name (18)  = 'P21212          ' 
      spcgr_name (19)  = 'P212121         ' 
      spcgr_name (20)  = 'C2221           ' 
      spcgr_name (21)  = 'C222            ' 
      spcgr_name (22)  = 'F222            ' 
      spcgr_name (23)  = 'I222            ' 
      spcgr_name (24)  = 'I212121         ' 
      spcgr_name (25)  = 'Pmm2            ' 
      spcgr_name (26)  = 'Pmc21           ' 
      spcgr_name (27)  = 'Pcc2            ' 
      spcgr_name (28)  = 'Pma2            ' 
      spcgr_name (29)  = 'Pca21           ' 
      spcgr_name (30)  = 'Pnc2            ' 
      spcgr_name (31)  = 'Pmn21           ' 
      spcgr_name (32)  = 'Pba2            ' 
      spcgr_name (33)  = 'Pna21           ' 
      spcgr_name (34)  = 'Pnn2            ' 
      spcgr_name (35)  = 'Cmm2            ' 
      spcgr_name (36)  = 'Cmc21           ' 
      spcgr_name (37)  = 'Ccc2            ' 
      spcgr_name (38)  = 'Amm2            ' 
      spcgr_name (39)  = 'Abm2            ' 
      spcgr_name (40)  = 'Ama2            ' 
      spcgr_name (41)  = 'Aba2            ' 
      spcgr_name (42)  = 'Fmm2            ' 
      spcgr_name (43)  = 'Fdd2            ' 
      spcgr_name (44)  = 'Imm2            ' 
      spcgr_name (45)  = 'Iba2            ' 
      spcgr_name (46)  = 'Ima2            ' 
      spcgr_name (47)  = 'Pmmm            ' 
                                                                        
      spcgr_name (48)  = 'Pnnn            ' 
      spcgr_num (48, 2) = 276 
      spcgr_name (276)  = 'Pnnn            ' 
      spcgr_syst (276) = ORTHORHOMBIC 
                                                                        
      spcgr_name (49)  = 'Pccm            ' 
                                                                        
      spcgr_name (50)  = 'Pban            ' 
      spcgr_num (50, 2) = 277 
      spcgr_name (277)  = 'Pban            ' 
      spcgr_syst (277) = ORTHORHOMBIC 
                                                                        
      spcgr_name (51)  = 'Pmma            ' 
      spcgr_name (52)  = 'Pnna            ' 
      spcgr_name (53)  = 'Pmna            ' 
      spcgr_name (54)  = 'Pcca            ' 
      spcgr_name (55)  = 'Pbam            ' 
      spcgr_name (56)  = 'Pccn            ' 
      spcgr_name (57)  = 'Pbcm            ' 
      spcgr_name (58)  = 'Pnnm            ' 
                                                                        
      spcgr_name (59)  = 'Pmmn            ' 
      spcgr_num (59, 2) = 278 
      spcgr_name (278)  = 'Pmmn            ' 
      spcgr_syst (278) = ORTHORHOMBIC 
                                                                        
      spcgr_name (60)  = 'Pbcn            ' 
      spcgr_name (61)  = 'Pbca            ' 
      spcgr_name (62)  = 'Pnma            ' 
      spcgr_name (63)  = 'Cmcm            ' 
      spcgr_name (64)  = 'Cmca            ' 
      spcgr_name (65)  = 'Cmmm            ' 
      spcgr_name (66)  = 'Cccm            ' 
      spcgr_name (67)  = 'Cmma            ' 
                                                                        
      spcgr_name (68)  = 'Ccca            ' 
      spcgr_num (68, 2) = 279 
      spcgr_name (279)  = 'Ccca            ' 
      spcgr_syst (279) = ORTHORHOMBIC 
                                                                        
      spcgr_name (69)  = 'Fmmm            ' 
                                                                        
      spcgr_name (70)  = 'Fddd            ' 
      spcgr_num (70, 2) = 280 
      spcgr_name (280)  = 'Fddd            ' 
      spcgr_syst (280) = ORTHORHOMBIC 
                                                                        
      spcgr_name (71)  = 'Immm            ' 
      spcgr_name (72)  = 'Ibam            ' 
      spcgr_name (73)  = 'Ibca            ' 
      spcgr_name (74)  = 'Imma            ' 
                                                                        
!                                                                       
!     The tetragonal space groups                                       
!                                                                       
      DO i = 75, 142 
      spcgr_syst (i) = TETRAGONAL 
      ENDDO 
                                                                        
      spcgr_name (75)  = 'P4              ' 
      spcgr_name (76)  = 'P41             ' 
      spcgr_name (77)  = 'P42             ' 
      spcgr_name (78)  = 'P43             ' 
      spcgr_name (79)  = 'I4              ' 
      spcgr_name (80)  = 'I41             ' 
      spcgr_name (81)  = 'P-4             ' 
      spcgr_name (82)  = 'I-4             ' 
      spcgr_name (83)  = 'P4/m            ' 
      spcgr_name (84)  = 'P42/m           ' 
                                                                        
      spcgr_name (85)  = 'P4/n            ' 
      spcgr_num (85, 2) = 281 
      spcgr_name (281)  = 'P4/n            ' 
      spcgr_syst (281) = TETRAGONAL 
                                                                        
      spcgr_name (86)  = 'P42/n           ' 
      spcgr_num (86, 2) = 282 
      spcgr_name (282)  = 'P42/n           ' 
      spcgr_syst (282) = TETRAGONAL 
                                                                        
      spcgr_name (87)  = 'I4/m            ' 
                                                                        
      spcgr_name (88)  = 'I41/a           ' 
      spcgr_num (88, 2) = 283 
      spcgr_name (283)  = 'I41/a           ' 
      spcgr_syst (283) = TETRAGONAL 
                                                                        
      spcgr_name (89)  = 'P422            ' 
      spcgr_name (90)  = 'P4212           ' 
      spcgr_name (91)  = 'P4122           ' 
      spcgr_name (92)  = 'P41212          ' 
      spcgr_name (93)  = 'P4222           ' 
      spcgr_name (94)  = 'P42212          ' 
      spcgr_name (95)  = 'P4322           ' 
      spcgr_name (96)  = 'P43212          ' 
      spcgr_name (97)  = 'I422            ' 
      spcgr_name (98)  = 'I4122           ' 
      spcgr_name (99)  = 'P4mm            ' 
      spcgr_name (100)  = 'P4bm            ' 
      spcgr_name (101)  = 'P42cm           ' 
      spcgr_name (102)  = 'P42nm           ' 
      spcgr_name (103)  = 'P4cc            ' 
      spcgr_name (104)  = 'P4nc            ' 
      spcgr_name (105)  = 'P42mc           ' 
      spcgr_name (106)  = 'P42bc           ' 
      spcgr_name (107)  = 'I4mm            ' 
      spcgr_name (108)  = 'I4cm            ' 
      spcgr_name (109)  = 'I41md           ' 
      spcgr_name (110)  = 'I41cd           ' 
      spcgr_name (111)  = 'P-42m           ' 
      spcgr_name (112)  = 'P-42c           ' 
      spcgr_name (113)  = 'P-421m          ' 
      spcgr_name (114)  = 'P-421c          ' 
      spcgr_name (115)  = 'P-4m2           ' 
      spcgr_name (116)  = 'P-4c2           ' 
      spcgr_name (117)  = 'P-4b2           ' 
      spcgr_name (118)  = 'P-4n2           ' 
      spcgr_name (119)  = 'I-4m2           ' 
      spcgr_name (120)  = 'I-4c2           ' 
      spcgr_name (121)  = 'I-42m           ' 
      spcgr_name (122)  = 'I-42d           ' 
      spcgr_name (123)  = 'P4/mmm          ' 
      spcgr_name (124)  = 'P4/mcc          ' 
                                                                        
      spcgr_name (125)  = 'P4/nbm          ' 
      spcgr_num (125, 2) = 284 
      spcgr_name (284)  = 'P4/nbm          ' 
      spcgr_syst (284) = TETRAGONAL 
                                                                        
      spcgr_name (126)  = 'P4/nnc          ' 
      spcgr_num (126, 2) = 285 
      spcgr_name (285)  = 'P4/nnc          ' 
      spcgr_syst (285) = TETRAGONAL 
                                                                        
      spcgr_name (127)  = 'P4/mbm          ' 
      spcgr_name (128)  = 'P4/mnc          ' 
                                                                        
      spcgr_name (129)  = 'P4/nmm          ' 
      spcgr_num (129, 2) = 286 
      spcgr_name (286)  = 'P4/nmm          ' 
      spcgr_syst (286) = TETRAGONAL 
                                                                        
      spcgr_name (130)  = 'P4/ncc          ' 
      spcgr_num (130, 2) = 287 
      spcgr_name (287)  = 'P4/ncc          ' 
      spcgr_syst (287) = TETRAGONAL 
                                                                        
      spcgr_name (131)  = 'P42/mmc         ' 
      spcgr_name (132)  = 'P42/mcm         ' 
                                                                        
      spcgr_name (133)  = 'P42/nbc         ' 
      spcgr_num (133, 2) = 288 
      spcgr_name (288)  = 'P42/nbc         ' 
      spcgr_syst (288) = TETRAGONAL 
                                                                        
      spcgr_name (134)  = 'P42/nnm         ' 
      spcgr_num (134, 2) = 289 
      spcgr_name (289)  = 'P42/nnm         ' 
      spcgr_syst (289) = TETRAGONAL 
                                                                        
      spcgr_name (135)  = 'P42/mbc         ' 
      spcgr_name (136)  = 'P42/mnm         ' 
                                                                        
      spcgr_name (137)  = 'P42/nmc         ' 
      spcgr_num (137, 2) = 290 
      spcgr_name (290)  = 'P42/nmc         ' 
      spcgr_syst (290) = TETRAGONAL 
                                                                        
      spcgr_name (138)  = 'P42/ncm         ' 
      spcgr_num (138, 2) = 291 
      spcgr_name (291)  = 'P42/ncm         ' 
      spcgr_syst (291) = TETRAGONAL 
                                                                        
      spcgr_name (139)  = 'I4/mmm          ' 
      spcgr_name (140)  = 'I4/mcm          ' 
                                                                        
      spcgr_name (141)  = 'I41/amd         ' 
      spcgr_num (141, 2) = 292 
      spcgr_name (292)  = 'I41/amd         ' 
      spcgr_syst (292) = TETRAGONAL 
                                                                        
      spcgr_name (142)  = 'I41/acd         ' 
      spcgr_num (142, 2) = 293 
      spcgr_name (293)  = 'I41/acd         ' 
      spcgr_syst (293) = TETRAGONAL 
                                                                        
!                                                                       
!     The trigonal space groups, including the rhombohedral in          
!       hexagonal cells                                                 
!                                                                       
      DO i = 143, 167 
      spcgr_syst (i) = TRIGONAL 
      ENDDO 
                                                                        
      spcgr_name (143)  = 'P3              ' 
      spcgr_name (144)  = 'P31             ' 
      spcgr_name (145)  = 'P32             ' 
                                                                        
      spcgr_name (146)  = 'R3              ' 
      spcgr_num (146, 2) = 294 
      spcgr_name (294)  = 'R3              ' 
      spcgr_syst (294) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (147)  = 'P-3             ' 
                                                                        
      spcgr_name (148)  = 'R-3             ' 
      spcgr_num (148, 2) = 295 
      spcgr_name (295)  = 'R-3             ' 
      spcgr_syst (295) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (149)  = 'P312            ' 
      spcgr_name (150)  = 'P321            ' 
      spcgr_name (151)  = 'P3112           ' 
      spcgr_name (152)  = 'P3121           ' 
      spcgr_name (153)  = 'P3212           ' 
      spcgr_name (154)  = 'P3221           ' 
                                                                        
      spcgr_name (155)  = 'R32             ' 
      spcgr_num (155, 2) = 296 
      spcgr_name (296)  = 'R32             ' 
      spcgr_syst (296) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (156)  = 'P3m1            ' 
      spcgr_name (157)  = 'P31m            ' 
      spcgr_name (158)  = 'P3c1            ' 
      spcgr_name (159)  = 'P31c            ' 
                                                                        
      spcgr_name (160)  = 'R3m             ' 
      spcgr_num (160, 2) = 297 
      spcgr_name (297)  = 'R3m             ' 
      spcgr_syst (297) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (161)  = 'R3c             ' 
      spcgr_num (161, 2) = 298 
      spcgr_name (298)  = 'R3c             ' 
      spcgr_syst (298) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (162)  = 'P-31m           ' 
      spcgr_name (163)  = 'P-31c           ' 
      spcgr_name (164)  = 'P-3m1           ' 
      spcgr_name (165)  = 'P-3c1           ' 
                                                                        
      spcgr_name (166)  = 'R-3m            ' 
      spcgr_num (166, 2) = 299 
      spcgr_name (299)  = 'R-3m            ' 
      spcgr_syst (299) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (167)  = 'R-3c            ' 
      spcgr_num (167, 2) = 300 
      spcgr_name (300)  = 'R-3c            ' 
      spcgr_syst (300) = RHOMBOHEDRAL 
                                                                        
!                                                                       
!     The hexagonal space groups                                        
!                                                                       
      DO i = 168, 194 
      spcgr_syst (i) = HEXAGONAL 
      ENDDO 
                                                                        
      spcgr_name (168)  = 'P6              ' 
      spcgr_name (169)  = 'P61             ' 
      spcgr_name (170)  = 'P65             ' 
      spcgr_name (171)  = 'P62             ' 
      spcgr_name (172)  = 'P64             ' 
      spcgr_name (173)  = 'P63             ' 
      spcgr_name (174)  = 'P-6             ' 
      spcgr_name (175)  = 'P6/m            ' 
      spcgr_name (176)  = 'P63/m           ' 
      spcgr_name (177)  = 'P622            ' 
      spcgr_name (178)  = 'P6122           ' 
      spcgr_name (179)  = 'P6522           ' 
      spcgr_name (180)  = 'P6222           ' 
      spcgr_name (181)  = 'P6422           ' 
      spcgr_name (182)  = 'P6322           ' 
      spcgr_name (183)  = 'P6mm            ' 
      spcgr_name (184)  = 'P6cc            ' 
      spcgr_name (185)  = 'P63cm           ' 
      spcgr_name (186)  = 'P63mc           ' 
      spcgr_name (187)  = 'P-6m2           ' 
      spcgr_name (188)  = 'P-6c2           ' 
      spcgr_name (189)  = 'P-62m           ' 
      spcgr_name (190)  = 'P-62c           ' 
      spcgr_name (191)  = 'P6/mmm          ' 
      spcgr_name (192)  = 'P6/mcc          ' 
      spcgr_name (193)  = 'P63/mcm         ' 
      spcgr_name (194)  = 'P63/mmc         ' 
                                                                        
!                                                                       
!     The cubic space groups                                            
!                                                                       
      DO i = 195, 230 
      spcgr_syst (i) = CUBIC 
      ENDDO 
                                                                        
      spcgr_name (195)  = 'P23             ' 
      spcgr_name (196)  = 'F23             ' 
      spcgr_name (197)  = 'I23             ' 
      spcgr_name (198)  = 'P213            ' 
      spcgr_name (199)  = 'I213            ' 
      spcgr_name (200)  = 'Pm-3            ' 
                                                                        
      spcgr_name (201)  = 'Pn-3            ' 
      spcgr_num (201, 2) = 301 
      spcgr_name (301)  = 'Pn-3            ' 
      spcgr_syst (301) = CUBIC 
                                                                        
      spcgr_name (202)  = 'Fm-3            ' 
                                                                        
      spcgr_name (203)  = 'Fd-3            ' 
      spcgr_num (203, 2) = 302 
      spcgr_name (302)  = 'Fd-3            ' 
      spcgr_syst (302) = CUBIC 
                                                                        
      spcgr_name (204)  = 'Im-3            ' 
      spcgr_name (205)  = 'Pa-3            ' 
      spcgr_name (206)  = 'Ia-3            ' 
      spcgr_name (207)  = 'P432            ' 
      spcgr_name (208)  = 'P4232           ' 
      spcgr_name (209)  = 'F432            ' 
      spcgr_name (210)  = 'F4132           ' 
      spcgr_name (211)  = 'I432            ' 
      spcgr_name (212)  = 'P4332           ' 
      spcgr_name (213)  = 'P4132           ' 
      spcgr_name (214)  = 'I4132           ' 
      spcgr_name (215)  = 'P-43m           ' 
      spcgr_name (216)  = 'F-43m           ' 
      spcgr_name (217)  = 'I-43m           ' 
      spcgr_name (218)  = 'P-43n           ' 
      spcgr_name (219)  = 'F-43c           ' 
      spcgr_name (220)  = 'I-43d           ' 
      spcgr_name (221)  = 'Pm-3m           ' 
                                                                        
      spcgr_name (222)  = 'Pn-3n           ' 
      spcgr_num (222, 2) = 303 
      spcgr_name (303)  = 'Pn-3n           ' 
      spcgr_syst (303) = CUBIC 
                                                                        
      spcgr_name (223)  = 'Pm-3n           ' 
                                                                        
      spcgr_name (224)  = 'Pn-3m           ' 
      spcgr_num (224, 2) = 304 
      spcgr_name (304)  = 'Pn-3m           ' 
      spcgr_syst (304) = CUBIC 
                                                                        
      spcgr_name (225)  = 'Fm-3m           ' 
      spcgr_name (226)  = 'Fm-3c           ' 
                                                                        
      spcgr_name (227)  = 'Fd-3m           ' 
      spcgr_num (227, 2) = 305 
      spcgr_name (305)  = 'Fd-3m           ' 
      spcgr_syst (305) = CUBIC 
                                                                        
      spcgr_name (228)  = 'Fd-3c           ' 
      spcgr_num (228, 2) = 306 
      spcgr_name (306)  = 'Fd-3c           ' 
      spcgr_syst (306) = CUBIC 
                                                                        
      spcgr_name (229)  = 'Im-3m           ' 
      spcgr_name (230)  = 'Ia-3d           ' 
                                                                        
!                                                                       
!     The monoclinic groups, c-axis unique, other cell choices etc...   
!                                                                       
                                                                        
      spcgr_name (231)  = 'P112            ' 
      spcgr_syst (231) = MONOCLINIC_C 
                                                                        
      spcgr_name (232)  = 'P1121           ' 
      spcgr_syst (232) = MONOCLINIC_C 
                                                                        
      spcgr_name (233)  = 'A121            ' 
      spcgr_syst (233) = MONOCLINIC_B 
                                                                        
      spcgr_name (234)  = 'I121            ' 
      spcgr_syst (234) = MONOCLINIC_B 
                                                                        
      spcgr_name (235)  = 'A112            ' 
      spcgr_syst (235) = MONOCLINIC_C 
                                                                        
      spcgr_name (236)  = 'B112            ' 
      spcgr_syst (236) = MONOCLINIC_C 
                                                                        
      spcgr_name (237)  = 'I112            ' 
      spcgr_syst (237) = MONOCLINIC_C 
                                                                        
      spcgr_name (238)  = 'P11m            ' 
      spcgr_syst (238) = MONOCLINIC_C 
                                                                        
      spcgr_name (239)  = 'P1n1            ' 
      spcgr_syst (239) = MONOCLINIC_B 
                                                                        
      spcgr_name (240)  = 'P1a1            ' 
      spcgr_syst (240) = MONOCLINIC_B 
                                                                        
      spcgr_name (241)  = 'P11a            ' 
      spcgr_syst (241) = MONOCLINIC_C 
                                                                        
      spcgr_name (242)  = 'P11n            ' 
      spcgr_syst (242) = MONOCLINIC_C 
                                                                        
      spcgr_name (243)  = 'P11b            ' 
      spcgr_syst (243) = MONOCLINIC_C 
                                                                        
      spcgr_name (244)  = 'A1m1            ' 
      spcgr_syst (244) = MONOCLINIC_B 
                                                                        
      spcgr_name (245)  = 'I1m1            ' 
      spcgr_syst (245) = MONOCLINIC_B 
                                                                        
      spcgr_name (246)  = 'A11m            ' 
      spcgr_syst (246) = MONOCLINIC_C 
                                                                        
      spcgr_name (247)  = 'B11m            ' 
      spcgr_syst (247) = MONOCLINIC_C 
                                                                        
      spcgr_name (248)  = 'I11m            ' 
      spcgr_syst (248) = MONOCLINIC_C 
                                                                        
      spcgr_name (249)  = 'A1n1            ' 
      spcgr_syst (259) = MONOCLINIC_B 
                                                                        
      spcgr_name (250)  = 'I1a1            ' 
      spcgr_syst (250) = MONOCLINIC_B 
                                                                        
      spcgr_name (251)  = 'A11a            ' 
      spcgr_syst (251) = MONOCLINIC_C 
                                                                        
      spcgr_name (252)  = 'B11n            ' 
      spcgr_syst (252) = MONOCLINIC_C 
                                                                        
      spcgr_name (253)  = 'I11b            ' 
      spcgr_syst (253) = MONOCLINIC_C 
                                                                        
      spcgr_name (254)  = 'P112/m          ' 
      spcgr_syst (254) = MONOCLINIC_C 
                                                                        
      spcgr_name (255)  = 'P1121/m         ' 
      spcgr_syst (255) = MONOCLINIC_C 
                                                                        
      spcgr_name (256)  = 'A12/m1          ' 
      spcgr_syst (256) = MONOCLINIC_B 
                                                                        
      spcgr_name (257)  = 'I12/m1          ' 
      spcgr_syst (257) = MONOCLINIC_B 
                                                                        
      spcgr_name (258)  = 'A112/m          ' 
      spcgr_syst (258) = MONOCLINIC_C 
                                                                        
      spcgr_name (259)  = 'B112/m          ' 
      spcgr_syst (259) = MONOCLINIC_C 
                                                                        
      spcgr_name (260)  = 'I112/m          ' 
      spcgr_syst (260) = MONOCLINIC_C 
                                                                        
      spcgr_name (261)  = 'P12/n1          ' 
      spcgr_syst (261) = MONOCLINIC_B 
                                                                        
      spcgr_name (262)  = 'P12/a1          ' 
      spcgr_syst (262) = MONOCLINIC_B 
                                                                        
      spcgr_name (263)  = 'P112/a          ' 
      spcgr_syst (263) = MONOCLINIC_C 
                                                                        
      spcgr_name (264)  = 'P112/n          ' 
      spcgr_syst (264) = MONOCLINIC_C 
                                                                        
      spcgr_name (265)  = 'P112/b          ' 
      spcgr_syst (265) = MONOCLINIC_C 
                                                                        
      spcgr_name (266)  = 'P121/n1         ' 
      spcgr_syst (266) = MONOCLINIC_B 
                                                                        
      spcgr_name (267)  = 'P121/a1         ' 
      spcgr_syst (267) = MONOCLINIC_B 
                                                                        
      spcgr_name (268)  = 'P1121/a         ' 
      spcgr_syst (268) = MONOCLINIC_C 
                                                                        
      spcgr_name (269)  = 'P1121/n         ' 
      spcgr_syst (269) = MONOCLINIC_C 
                                                                        
      spcgr_name (270)  = 'P1121/b         ' 
      spcgr_syst (270) = MONOCLINIC_C 
                                                                        
      spcgr_name (271)  = 'A12/n1          ' 
      spcgr_syst (271) = MONOCLINIC_B 
                                                                        
      spcgr_name (272)  = 'I12/a1          ' 
      spcgr_syst (272) = MONOCLINIC_B 
                                                                        
      spcgr_name (273)  = 'A112/a          ' 
      spcgr_syst (273) = MONOCLINIC_C 
                                                                        
      spcgr_name (274)  = 'B112/n          ' 
      spcgr_syst (274) = MONOCLINIC_C 
                                                                        
      spcgr_name (275)  = 'I112/b          ' 
      spcgr_syst (275) = MONOCLINIC_C 
                                                                        
      spcgr_name (307)  = 'Pbnm            ' 
      spcgr_syst (307) = ORTHORHOMBIC 
                                                                        
      spcgr_name (308)  = 'Pmnn            ' 
      spcgr_syst (308) = ORTHORHOMBIC 
                                                                        
      spcgr_name (309)  = 'Ibmm            ' 
      spcgr_syst (309) = ORTHORHOMBIC 

      spcgr_name (310)  = 'Aem2            ' 
      spcgr_num (310, 1) =  39 
      spcgr_syst (310) = ORTHORHOMBIC 

      spcgr_name (311)  = 'Aea2            ' 
      spcgr_num (311, 1) =  41 
      spcgr_syst (311) = ORTHORHOMBIC 

      spcgr_name (312)  = 'Cmce            ' 
      spcgr_num (312, 1) =  64 
      spcgr_syst (312) = ORTHORHOMBIC 

      spcgr_name (313)  = 'Cmme            ' 
      spcgr_num (313, 1) =  67 
      spcgr_syst (313) = ORTHORHOMBIC 

      spcgr_name (314)  = 'Ccce            ' 
      spcgr_num (314, 1) =  68 
      spcgr_syst (314) = ORTHORHOMBIC 
                                                                        
      END SUBROUTINE spcgr_setup                    
END MODULE spcgr_setup_mod
