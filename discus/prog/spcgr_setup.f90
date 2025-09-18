MODULE spcgr_setup_mod
!
CONTAINS
!********************************************************************** 
      SUBROUTINE spcgr_setup 
!-                                                                      
!     Sets up the space group symbol tables.                            
!+                                                                      
      USE spcgr_mod 
use blanks_mod
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
      INTEGER i, l, k
!logical       :: lall
!integer, save :: j = 0
character(len=16), dimension(1:SPCGR_MAX) :: spcgr_point_t
!                                                                       
DO i = 1, SPCGR_MAX 
   spcgr_num (i, 1) = i 
   spcgr_num (i, 2) = i 
ENDDO 
spcgr_name = ' '
!                                                                       
!     The triclinic space groups                                        
!                                                                       
      spcgr_name (1)  = 'P 1           ' 
      spcgr_syst (1) = TRICLINIC 
                                                                        
      spcgr_name (2)  = 'P -1          ' 
      spcgr_syst (2) = TRICLINIC 
!                                                                       
!     The monoclinic space groups                                       
!                                                                       
      DO i = 3, 15 
      spcgr_syst (i) = MONOCLINIC_B 
      ENDDO 
                                                                        
      spcgr_name (3)  = 'P 2           ' 
      spcgr_name (4)  = 'P 21          ' 
      spcgr_name (5)  = 'C 2           ' 
      spcgr_name (6)  = 'P m           ' 
      spcgr_name (7)  = 'P c           ' 
      spcgr_name (8)  = 'C m           ' 
      spcgr_name (9)  = 'C c           ' 
      spcgr_name (10)  = 'P 2/m         ' 
      spcgr_name (11)  = 'P 21/m        ' 
      spcgr_name (12)  = 'C 2/m         ' 
      spcgr_name (13)  = 'P 2/c         ' 
      spcgr_name (14)  = 'P 21/c        ' 
      spcgr_name (15)  = 'C 2/c         ' 
                                                                        
!                                                                       
!     The orthorhombic space groups                                     
!                                                                       
      DO i = 16, 74 
      spcgr_syst (i) = ORTHORHOMBIC 
      ENDDO 
                                                                        
      spcgr_name (16)  = 'P 2 2 2         ' 
      spcgr_name (17)  = 'P 2 2 21        ' 
      spcgr_name (18)  = 'P 21 21 2       ' 
      spcgr_name (19)  = 'P 21 21 21      ' 
      spcgr_name (20)  = 'C 2 2 21        ' 
      spcgr_name (21)  = 'C 2 2 2         ' 
      spcgr_name (22)  = 'F 2 2 2         ' 
      spcgr_name (23)  = 'I 2 2 2         ' 
      spcgr_name (24)  = 'I 21 21 21      ' 
      spcgr_name (25)  = 'P m m 2         ' 
      spcgr_name (26)  = 'P m c 21        ' 
      spcgr_name (27)  = 'P c c 2         ' 
      spcgr_name (28)  = 'P m a 2         ' 
      spcgr_name (29)  = 'P c a 21        ' 
      spcgr_name (30)  = 'P n c 2         ' 
      spcgr_name (31)  = 'P m n 21        ' 
      spcgr_name (32)  = 'P b a 2         ' 
      spcgr_name (33)  = 'P n a 21        ' 
      spcgr_name (34)  = 'P n n 2         ' 
      spcgr_name (35)  = 'C m m 2         ' 
      spcgr_name (36)  = 'C m c 21        ' 
      spcgr_name (37)  = 'C c c 2         ' 
      spcgr_name (38)  = 'A m m 2         ' 
      spcgr_name (39)  = 'A b m 2         ' 
      spcgr_name (40)  = 'A m a 2         ' 
      spcgr_name (41)  = 'A b a 2         ' 
      spcgr_name (42)  = 'F m m 2         ' 
      spcgr_name (43)  = 'F d d 2         ' 
      spcgr_name (44)  = 'I m m 2         ' 
      spcgr_name (45)  = 'I b a 2         ' 
      spcgr_name (46)  = 'I m a 2         ' 
      spcgr_name (47)  = 'P m m m         ' 
                                                                        
      spcgr_name (48)  = 'P n n n         ' 
      spcgr_num (48, 2) = 276 
      spcgr_name (276)  = 'P n n n         ' 
      spcgr_syst (276) = ORTHORHOMBIC 
                                                                        
      spcgr_name (49)  = 'P c c m         ' 
                                                                        
      spcgr_name (50)  = 'P b a n         ' 
      spcgr_num (50, 2) = 277 
      spcgr_name (277)  = 'P b a n         ' 
      spcgr_syst (277) = ORTHORHOMBIC 
                                                                        
      spcgr_name (51)  = 'P m m a         ' 
      spcgr_name (52)  = 'P n n a         ' 
      spcgr_name (53)  = 'P m n a         ' 
      spcgr_name (54)  = 'P c c a         ' 
      spcgr_name (55)  = 'P b a m         ' 
      spcgr_name (56)  = 'P c c n         ' 
      spcgr_name (57)  = 'P b c m         ' 
      spcgr_name (58)  = 'P n n m         ' 
                                                                        
      spcgr_name (59)  = 'P m m n         ' 
      spcgr_num (59, 2) = 278 
      spcgr_name (278)  = 'P m m n         ' 
      spcgr_syst (278) = ORTHORHOMBIC 
                                                                        
      spcgr_name (60)  = 'P b c n         ' 
      spcgr_name (61)  = 'P b c a         ' 
      spcgr_name (62)  = 'P n m a         ' 
      spcgr_name (63)  = 'C m c m         ' 
      spcgr_name (64)  = 'C m c a         ' 
      spcgr_name (65)  = 'C m m m         ' 
      spcgr_name (66)  = 'C c c m         ' 
      spcgr_name (67)  = 'C m m a         ' 
                                                                        
      spcgr_name (68)  = 'C c c a         ' 
      spcgr_num (68, 2) = 279 
      spcgr_name (279)  = 'C c c a         ' 
      spcgr_syst (279) = ORTHORHOMBIC 
                                                                        
      spcgr_name (69)  = 'F m m m         ' 
                                                                        
      spcgr_name (70)  = 'F d d d         ' 
      spcgr_num (70, 2) = 280 
      spcgr_name (280)  = 'F d d d         ' 
      spcgr_syst (280) = ORTHORHOMBIC 
                                                                        
      spcgr_name (71)  = 'I m m m         ' 
      spcgr_name (72)  = 'I b a m         ' 
      spcgr_name (73)  = 'I b c a         ' 
      spcgr_name (74)  = 'I m m a         ' 
                                                                        
!                                                                       
!     The tetragonal space groups                                       
!                                                                       
      DO i = 75, 142 
      spcgr_syst (i) = TETRAGONAL 
      ENDDO 
                                                                        
      spcgr_name (75)  = 'P 4           ' 
      spcgr_name (76)  = 'P 41          ' 
      spcgr_name (77)  = 'P 42          ' 
      spcgr_name (78)  = 'P 43          ' 
      spcgr_name (79)  = 'I 4           ' 
      spcgr_name (80)  = 'I 41          ' 
      spcgr_name (81)  = 'P -4          ' 
      spcgr_name (82)  = 'I -4          ' 
      spcgr_name (83)  = 'P 4/m         ' 
      spcgr_name (84)  = 'P 42/m        ' 
                                                                        
      spcgr_name (85)  = 'P 4/n         ' 
      spcgr_num (85, 2) = 281 
      spcgr_name (281)  = 'P 4/n         ' 
      spcgr_syst (281) = TETRAGONAL 
                                                                        
      spcgr_name (86)  = 'P 42/n        ' 
      spcgr_num (86, 2) = 282 
      spcgr_name (282)  = 'P 42/n        ' 
      spcgr_syst (282) = TETRAGONAL 
                                                                        
      spcgr_name (87)  = 'I 4/m         ' 
                                                                        
      spcgr_name (88)  = 'I 41/a        ' 
      spcgr_num (88, 2) = 283 
      spcgr_name (283)  = 'I 41/a        ' 
      spcgr_syst (283) = TETRAGONAL 
                                                                        
      spcgr_name (89)  = 'P 4 2 2         ' 
      spcgr_name (90)  = 'P 4 21 2        ' 
      spcgr_name (91)  = 'P 41 2 2        ' 
      spcgr_name (92)  = 'P 41 21 2       ' 
      spcgr_name (93)  = 'P 42 2 2        ' 
      spcgr_name (94)  = 'P 42 21 2       ' 
      spcgr_name (95)  = 'P 43 2 2        ' 
      spcgr_name (96)  = 'P 43 21 2       ' 
      spcgr_name (97)  = 'I 4 2 2         ' 
      spcgr_name (98)  = 'I 41 2 2        ' 
      spcgr_name (99)  = 'P 4 m m         ' 
      spcgr_name (100)  = 'P 4 b m         ' 
      spcgr_name (101)  = 'P 42 c m        ' 
      spcgr_name (102)  = 'P 42 n m        ' 
      spcgr_name (103)  = 'P 4 c c         ' 
      spcgr_name (104)  = 'P 4 n c         ' 
      spcgr_name (105)  = 'P 42 m c        ' 
      spcgr_name (106)  = 'P 42 b c        ' 
      spcgr_name (107)  = 'I 4 m m         ' 
      spcgr_name (108)  = 'I 4 c m         ' 
      spcgr_name (109)  = 'I 41 m d        ' 
      spcgr_name (110)  = 'I 41 c d        ' 
      spcgr_name (111)  = 'P -4 2 m        ' 
      spcgr_name (112)  = 'P -4 2 c        ' 
      spcgr_name (113)  = 'P -4 21 m       ' 
      spcgr_name (114)  = 'P -4 21 c       ' 
      spcgr_name (115)  = 'P -4 m 2        ' 
      spcgr_name (116)  = 'P -4 c 2        ' 
      spcgr_name (117)  = 'P -4 b 2        ' 
      spcgr_name (118)  = 'P -4 n 2        ' 
      spcgr_name (119)  = 'I -4 m 2        ' 
      spcgr_name (120)  = 'I -4 c 2        ' 
      spcgr_name (121)  = 'I -4 2 m        ' 
      spcgr_name (122)  = 'I -4 2 d        ' 
      spcgr_name (123)  = 'P 4/m m m       ' 
      spcgr_name (124)  = 'P 4/m c c       ' 
                                                                        
      spcgr_name (125)  = 'P 4/n b m       ' 
      spcgr_num (125, 2) = 284 
      spcgr_name (284)  = 'P 4/n b m       ' 
      spcgr_syst (284) = TETRAGONAL 
                                                                        
      spcgr_name (126)  = 'P 4/n n c       ' 
      spcgr_num (126, 2) = 285 
      spcgr_name (285)  = 'P 4/n n c       ' 
      spcgr_syst (285) = TETRAGONAL 
                                                                        
      spcgr_name (127)  = 'P 4/m b m       ' 
      spcgr_name (128)  = 'P 4/m n c       ' 
                                                                        
      spcgr_name (129)  = 'P 4/n m m       ' 
      spcgr_num (129, 2) = 286 
      spcgr_name (286)  = 'P 4/n m m       ' 
      spcgr_syst (286) = TETRAGONAL 
                                                                        
      spcgr_name (130)  = 'P 4/n c c       ' 
      spcgr_num (130, 2) = 287 
      spcgr_name (287)  = 'P 4/n c c       ' 
      spcgr_syst (287) = TETRAGONAL 
                                                                        
      spcgr_name (131)  = 'P 42/m m c      ' 
      spcgr_name (132)  = 'P 42/m c m      ' 
                                                                        
      spcgr_name (133)  = 'P 42/n b c      ' 
      spcgr_num (133, 2) = 288 
      spcgr_name (288)  = 'P 42/n b c      ' 
      spcgr_syst (288) = TETRAGONAL 
                                                                        
      spcgr_name (134)  = 'P 42/n n m      ' 
      spcgr_num (134, 2) = 289 
      spcgr_name (289)  = 'P 42/n n m      ' 
      spcgr_syst (289) = TETRAGONAL 
                                                                        
      spcgr_name (135)  = 'P 42/m b c      ' 
      spcgr_name (136)  = 'P 42/m n m      ' 
                                                                        
      spcgr_name (137)  = 'P 42/n m c      ' 
      spcgr_num (137, 2) = 290 
      spcgr_name (290)  = 'P 42/n m c      ' 
      spcgr_syst (290) = TETRAGONAL 
                                                                        
      spcgr_name (138)  = 'P 42/n c m      ' 
      spcgr_num (138, 2) = 291 
      spcgr_name (291)  = 'P 42/n c m      ' 
      spcgr_syst (291) = TETRAGONAL 
                                                                        
      spcgr_name (139)  = 'I 4/m m m       ' 
      spcgr_name (140)  = 'I 4/m c m       ' 
                                                                        
      spcgr_name (141)  = 'I 41/a m d      ' 
      spcgr_num (141, 2) = 292 
      spcgr_name (292)  = 'I 41/a m d      ' 
      spcgr_syst (292) = TETRAGONAL 
                                                                        
      spcgr_name (142)  = 'I 41/a c d      ' 
      spcgr_num (142, 2) = 293 
      spcgr_name (293)  = 'I 41/a c d      ' 
      spcgr_syst (293) = TETRAGONAL 
                                                                        
!                                                                       
!     The trigonal space groups, including the rhombohedral in          
!       hexagonal cells                                                 
!                                                                       
      DO i = 143, 167 
      spcgr_syst (i) = TRIGONAL 
      ENDDO 
                                                                        
      spcgr_name (143)  = 'P 3           ' 
      spcgr_name (144)  = 'P 31          ' 
      spcgr_name (145)  = 'P 32          ' 
                                                                        
      spcgr_name (146)  = 'R 3           ' 
      spcgr_num (146, 2) = 294 
      spcgr_name (294)  = 'R 3           ' 
      spcgr_syst (294) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (147)  = 'P -3          ' 
                                                                        
      spcgr_name (148)  = 'R -3          ' 
      spcgr_num (148, 2) = 295 
      spcgr_name (295)  = 'R -3          ' 
      spcgr_syst (295) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (149)  = 'P 3 1 2         ' 
      spcgr_name (150)  = 'P 3 2 1         ' 
      spcgr_name (151)  = 'P 31 1 2        ' 
      spcgr_name (152)  = 'P 31 2 1        ' 
      spcgr_name (153)  = 'P 32 1 2        ' 
      spcgr_name (154)  = 'P 32 2 1        ' 
                                                                        
      spcgr_name (155)  = 'R 3 2          ' 
      spcgr_num (155, 2) = 296 
      spcgr_name (296)  = 'R 3 2          ' 
      spcgr_syst (296) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (156)  = 'P 3 m 1         ' 
      spcgr_name (157)  = 'P 3 1 m         ' 
      spcgr_name (158)  = 'P 3 c 1         ' 
      spcgr_name (159)  = 'P 3 1 c         ' 
                                                                        
      spcgr_name (160)  = 'R 3 m          ' 
      spcgr_num (160, 2) = 297 
      spcgr_name (297)  = 'R 3 m          ' 
      spcgr_syst (297) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (161)  = 'R 3 c          ' 
      spcgr_num (161, 2) = 298 
      spcgr_name (298)  = 'R 3 c          ' 
      spcgr_syst (298) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (162)  = 'P -3 1 m        ' 
      spcgr_name (163)  = 'P -3 1 c        ' 
      spcgr_name (164)  = 'P -3 m 1        ' 
      spcgr_name (165)  = 'P -3 c 1        ' 
                                                                        
      spcgr_name (166)  = 'R -3 m         ' 
      spcgr_num (166, 2) = 299 
      spcgr_name (299)  = 'R -3 m         ' 
      spcgr_syst (299) = RHOMBOHEDRAL 
                                                                        
      spcgr_name (167)  = 'R -3 c         ' 
      spcgr_num (167, 2) = 300 
      spcgr_name (300)  = 'R -3 c         ' 
      spcgr_syst (300) = RHOMBOHEDRAL 
                                                                        
!                                                                       
!     The hexagonal space groups                                        
!                                                                       
      DO i = 168, 194 
      spcgr_syst (i) = HEXAGONAL 
      ENDDO 
                                                                        
      spcgr_name (168)  = 'P 6           ' 
      spcgr_name (169)  = 'P 61          ' 
      spcgr_name (170)  = 'P 65          ' 
      spcgr_name (171)  = 'P 62          ' 
      spcgr_name (172)  = 'P 64          ' 
      spcgr_name (173)  = 'P 63          ' 
      spcgr_name (174)  = 'P -6          ' 
      spcgr_name (175)  = 'P 6/m         ' 
      spcgr_name (176)  = 'P 63/m        ' 
      spcgr_name (177)  = 'P 6 2 2         ' 
      spcgr_name (178)  = 'P 61 2 2        ' 
      spcgr_name (179)  = 'P 65 2 2        ' 
      spcgr_name (180)  = 'P 62 2 2        ' 
      spcgr_name (181)  = 'P 64 2 2        ' 
      spcgr_name (182)  = 'P 63 2 2        ' 
      spcgr_name (183)  = 'P 6 m m         ' 
      spcgr_name (184)  = 'P 6 c c         ' 
      spcgr_name (185)  = 'P 63 c m        ' 
      spcgr_name (186)  = 'P 63 m c        ' 
      spcgr_name (187)  = 'P -6 m 2        ' 
      spcgr_name (188)  = 'P -6 c 2        ' 
      spcgr_name (189)  = 'P -6 2 m        ' 
      spcgr_name (190)  = 'P -6 2 c        ' 
      spcgr_name (191)  = 'P 6/m m m       ' 
      spcgr_name (192)  = 'P 6/m c c       ' 
      spcgr_name (193)  = 'P 63/m c m      ' 
      spcgr_name (194)  = 'P 63/m m c      ' 
                                                                        
!                                                                       
!     The cubic space groups                                            
!                                                                       
      DO i = 195, 230 
      spcgr_syst (i) = CUBIC 
      ENDDO 
                                                                        
      spcgr_name (195)  = 'P 2 3          ' 
      spcgr_name (196)  = 'F 2 3          ' 
      spcgr_name (197)  = 'I 2 3          ' 
      spcgr_name (198)  = 'P 21 3         ' 
      spcgr_name (199)  = 'I 21 3         ' 
      spcgr_name (200)  = 'P m -3         ' 
                                                                        
      spcgr_name (201)  = 'P n -3         ' 
      spcgr_num (201, 2) = 301 
      spcgr_name (301)  = 'P n -3         ' 
      spcgr_syst (301) = CUBIC 
                                                                        
      spcgr_name (202)  = 'F m -3         ' 
                                                                        
      spcgr_name (203)  = 'F d -3         ' 
      spcgr_num (203, 2) = 302 
      spcgr_name (302)  = 'F d -3         ' 
      spcgr_syst (302) = CUBIC 
                                                                        
      spcgr_name (204)  = 'I m -3         ' 
      spcgr_name (205)  = 'P a -3         ' 
      spcgr_name (206)  = 'I a -3         ' 
      spcgr_name (207)  = 'P 4 3 2         ' 
      spcgr_name (208)  = 'P 42 3 2        ' 
      spcgr_name (209)  = 'F 4 3 2         ' 
      spcgr_name (210)  = 'F 41 3 2        ' 
      spcgr_name (211)  = 'I 4 3 2         ' 
      spcgr_name (212)  = 'P 43 3 2        ' 
      spcgr_name (213)  = 'P 41 3 2        ' 
      spcgr_name (214)  = 'I 41 3 2        ' 
      spcgr_name (215)  = 'P -4 3 m        ' 
      spcgr_name (216)  = 'F -4 3 m        ' 
      spcgr_name (217)  = 'I -4 3 m        ' 
      spcgr_name (218)  = 'P -4 3 n        ' 
      spcgr_name (219)  = 'F -4 3 c        ' 
      spcgr_name (220)  = 'I -4 3 d        ' 
      spcgr_name (221)  = 'P m -3 m        ' 
                                                                        
      spcgr_name (222)  = 'P n -3 n        ' 
      spcgr_num (222, 2) = 303 
      spcgr_name (303)  = 'P n -3 n        ' 
      spcgr_syst (303) = CUBIC 
                                                                        
      spcgr_name (223)  = 'P m -3 n        ' 
                                                                        
      spcgr_name (224)  = 'P n -3 m        ' 
      spcgr_num (224, 2) = 304 
      spcgr_name (304)  = 'P n -3 m        ' 
      spcgr_syst (304) = CUBIC 
                                                                        
      spcgr_name (225)  = 'F m -3 m        ' 
      spcgr_name (226)  = 'F m -3 c        ' 
                                                                        
      spcgr_name (227)  = 'F d -3 m        ' 
      spcgr_num (227, 2) = 305 
      spcgr_name (305)  = 'F d -3 m        ' 
      spcgr_syst (305) = CUBIC 
                                                                        
      spcgr_name (228)  = 'F d -3 c        ' 
      spcgr_num (228, 2) = 306 
      spcgr_name (306)  = 'F d -3 c        ' 
      spcgr_syst (306) = CUBIC 
                                                                        
      spcgr_name (229)  = 'I m -3 m        ' 
      spcgr_name (230)  = 'I a -3 d        ' 
                                                                        
!                                                                       
!     The monoclinic groups, c-axis unique, other cell choices etc...   
!                                                                       
                                                                        
      spcgr_name (231)  = 'P 1 1 2         ' 
      spcgr_syst (231) = MONOCLINIC_C 
                                                                        
      spcgr_name (232)  = 'P 1 1 21        ' 
      spcgr_syst (232) = MONOCLINIC_C 
                                                                        
      spcgr_name (233)  = 'A 1 2 1         ' 
      spcgr_syst (233) = MONOCLINIC_B 
                                                                        
      spcgr_name (234)  = 'I 1 2 1         ' 
      spcgr_syst (234) = MONOCLINIC_B 
                                                                        
      spcgr_name (235)  = 'A 1 1 2         ' 
      spcgr_syst (235) = MONOCLINIC_C 
                                                                        
      spcgr_name (236)  = 'B 1 1 2         ' 
      spcgr_syst (236) = MONOCLINIC_C 
                                                                        
      spcgr_name (237)  = 'I 1 1 2         ' 
      spcgr_syst (237) = MONOCLINIC_C 
                                                                        
      spcgr_name (238)  = 'P 1 1 m         ' 
      spcgr_syst (238) = MONOCLINIC_C 
                                                                        
      spcgr_name (239)  = 'P 1 n 1         ' 
      spcgr_syst (239) = MONOCLINIC_B 
                                                                        
      spcgr_name (240)  = 'P 1 a 1         ' 
      spcgr_syst (240) = MONOCLINIC_B 
                                                                        
      spcgr_name (241)  = 'P 1 1 a         ' 
      spcgr_syst (241) = MONOCLINIC_C 
                                                                        
      spcgr_name (242)  = 'P 1 1 n         ' 
      spcgr_syst (242) = MONOCLINIC_C 
                                                                         
      spcgr_name (243)  = 'P 1 1 b         ' 
      spcgr_syst (243) = MONOCLINIC_C 
                                                                        
      spcgr_name (244)  = 'A 1 m 1         ' 
      spcgr_syst (244) = MONOCLINIC_B 
                                                                        
      spcgr_name (245)  = 'I 1 m 1         ' 
      spcgr_syst (245) = MONOCLINIC_B 
                                                                        
      spcgr_name (246)  = 'A 1 1 m         ' 
      spcgr_syst (246) = MONOCLINIC_C 
                                                                        
      spcgr_name (247)  = 'B 1 1 m         ' 
      spcgr_syst (247) = MONOCLINIC_C 
                                                                        
      spcgr_name (248)  = 'I 1 1 m         ' 
      spcgr_syst (248) = MONOCLINIC_C 
                                                                        
      spcgr_name (249)  = 'A 1 n 1         ' 
      spcgr_syst (249) = MONOCLINIC_B 
                                                                        
      spcgr_name (250)  = 'I 1 a 1         ' 
      spcgr_syst (250) = MONOCLINIC_B 
                                                                        
      spcgr_name (251)  = 'A 1 1 a         ' 
      spcgr_syst (251) = MONOCLINIC_C 
                                                                        
      spcgr_name (252)  = 'B 1 1 n         ' 
      spcgr_syst (252) = MONOCLINIC_C 
                                                                        
      spcgr_name (253)  = 'I 1 1 b         ' 
      spcgr_syst (253) = MONOCLINIC_C 
                                                                        
      spcgr_name (254)  = 'P 1 1 2/m       ' 
      spcgr_syst (254) = MONOCLINIC_C 
                                                                        
      spcgr_name (255)  = 'P 1 1 21/m      ' 
      spcgr_syst (255) = MONOCLINIC_C 
                                                                        
      spcgr_name (256)  = 'A 1 2/m 1       ' 
      spcgr_syst (256) = MONOCLINIC_B 
                                                                        
      spcgr_name (257)  = 'I 1 2/m 1       ' 
      spcgr_syst (257) = MONOCLINIC_B 
                                                                        
      spcgr_name (258)  = 'A 1 1 2/m       ' 
      spcgr_syst (258) = MONOCLINIC_C 
                                                                        
      spcgr_name (259)  = 'B 1 1 2/m       ' 
      spcgr_syst (259) = MONOCLINIC_C 
                                                                        
      spcgr_name (260)  = 'I 1 1 2/m       ' 
      spcgr_syst (260) = MONOCLINIC_C 
                                                                        
      spcgr_name (261)  = 'P 1 2/n 1       ' 
      spcgr_syst (261) = MONOCLINIC_B 
                                                                        
      spcgr_name (262)  = 'P 1 2/a 1       ' 
      spcgr_syst (262) = MONOCLINIC_B 
                                                                        
      spcgr_name (263)  = 'P 1 1 2/a       ' 
      spcgr_syst (263) = MONOCLINIC_C 
                                                                        
      spcgr_name (264)  = 'P 1 1 2/n       ' 
      spcgr_syst (264) = MONOCLINIC_C 
                                                                        
      spcgr_name (265)  = 'P 1 1 2/b       ' 
      spcgr_syst (265) = MONOCLINIC_C 
                                                                        
      spcgr_name (266)  = 'P 1 21/n 1      ' 
      spcgr_syst (266) = MONOCLINIC_B 
                                                                        
      spcgr_name (267)  = 'P 1 21/a 1      ' 
      spcgr_syst (267) = MONOCLINIC_B 
                                                                        
      spcgr_name (268)  = 'P 1 1 21/a      ' 
      spcgr_syst (268) = MONOCLINIC_C 
                                                                        
      spcgr_name (269)  = 'P 1 1 21/n      ' 
      spcgr_syst (269) = MONOCLINIC_C 
                                                                        
      spcgr_name (270)  = 'P 1 1 21/b      ' 
      spcgr_syst (270) = MONOCLINIC_C 
                                                                        
      spcgr_name (271)  = 'A 1 2/n 1       ' 
      spcgr_syst (271) = MONOCLINIC_B 
                                                                        
      spcgr_name (272)  = 'I 1 2/a 1       ' 
      spcgr_syst (272) = MONOCLINIC_B 
                                                                        
      spcgr_name (273)  = 'A 1 1 2/a       ' 
      spcgr_syst (273) = MONOCLINIC_C 
                                                                        
      spcgr_name (274)  = 'B 1 1 2/n       ' 
      spcgr_syst (274) = MONOCLINIC_C 
                                                                        
      spcgr_name (275)  = 'I 1 1 2/b       ' 
      spcgr_syst (275) = MONOCLINIC_C 
                                                                        
      spcgr_name (307)  = 'P b n m         ' 
      spcgr_syst (307) = ORTHORHOMBIC 
                                                                        
      spcgr_name (308)  = 'P m n n         ' 
      spcgr_syst (308) = ORTHORHOMBIC 
                                                                        
      spcgr_name (309)  = 'I b m m         ' 
      spcgr_syst (309) = ORTHORHOMBIC 

      spcgr_name (310)  = 'A e m 2         ' 
      spcgr_num (310, 1) =  39 
      spcgr_syst (310) = ORTHORHOMBIC 

      spcgr_name (311)  = 'A e a 2         ' 
      spcgr_num (311, 1) =  41 
      spcgr_syst (311) = ORTHORHOMBIC 

      spcgr_name (312)  = 'C m c e         ' 
      spcgr_num (312, 1) =  64 
      spcgr_syst (312) = ORTHORHOMBIC 

      spcgr_name (313)  = 'C m m e         ' 
      spcgr_num (313, 1) =  67 
      spcgr_syst (313) = ORTHORHOMBIC 

      spcgr_name (314)  = 'C c c e         ' 
      spcgr_num (314, 1) =  68 
      spcgr_num (314, 2) = 279 
      spcgr_name (315)  = 'C c c a         ' 
      spcgr_num (314, 1) =  68 
      spcgr_num (315, 2) = 279 
      spcgr_syst (314) = ORTHORHOMBIC 
      spcgr_syst (315) = ORTHORHOMBIC 
!
! create Point group name 
!
! TRICLINIC
spcgr_laue(1:2) = '-1'
spcgr_point(1)  =  '1'
spcgr_point(2)  = '-1'
! MONOCLINIC
spcgr_laue(3:15)   = '2/m'
spcgr_point(3:5)   = '2'
spcgr_point(6:9)   = 'm'
spcgr_point(10:15) = '2/m'
! ORTHORHOMBIC
spcgr_laue(16:74)  = 'm m m'
spcgr_point(16:24) = '2 2 2'
spcgr_point(25:46) = 'm m 2'
spcgr_point(47:74) = 'm m m'
! TETRAGONAL
spcgr_laue(75:88)  = '4/m'
spcgr_point(75:80) = '4'
spcgr_point(81:82) = '-4'
spcgr_point(83:88) = '4/m'
spcgr_laue(89:142) = '4/m m m'
spcgr_point(89:98) = '4 2 2'
spcgr_point(99:110) = '4 m m'
spcgr_point(111:124) = '-4 2 m'
spcgr_point(115:120) = '-4 m 2'
spcgr_point(121:122) = '-4 2 m'
spcgr_point(123:142) = '4/m m m'
! Trigonal
spcgr_laue(143:148)  = '-3'
spcgr_point(143:146) = '3'
spcgr_point(147:148) = '-3'
spcgr_laue(149:167)  = '-3m'
spcgr_point(149    ) = '3 1 2'
spcgr_point(150    ) = '3 2'
spcgr_point(151    ) = '3 1 2'
spcgr_point(152    ) = '3 2'
spcgr_point(153    ) = '3 1 2'
spcgr_point(154    ) = '3 2'
spcgr_point(155    ) = '3 2'
!
spcgr_point(156    ) = '3 m'
spcgr_point(157    ) = '3 1 m'
spcgr_point(158    ) = '3 m'
spcgr_point(159    ) = '3 1 m'
spcgr_point(160    ) = '3 m'
spcgr_point(161    ) = '3 m'
spcgr_point(162    ) = '-3 1 m'
spcgr_point(163    ) = '-3 1 m'
spcgr_point(164:167) = '-3 m'
! HEXAGONAL
spcgr_laue(168:194)  = '6/m'
spcgr_point(168:173) = '6'
spcgr_point(174:174) = '-6'
spcgr_point(175:176) = '6/m'
spcgr_laue(177:194)  = '6/m m m'
spcgr_point(177:182) = '6 2 2'
spcgr_point(183:186) = '6 m m'
spcgr_point(187:188) = '-6 m 2'
spcgr_point(189:190) = '-6 2 m'
spcgr_point(191:194) = '6/m m m'
! CUBIC
spcgr_laue(195:206)  = 'm -3'
spcgr_point(195:199)  = '2 3'
spcgr_point(200:206) = 'm -3'
spcgr_laue(207:230)  = 'm -3 m'
spcgr_point(207:214) = '4 3 2'
spcgr_point(215:222) = '-4 3 m'
spcgr_point(221:230) = 'm -3 m'
!
! SPECIAL
spcgr_laue(231:232)  = '1 1 2/m'
spcgr_point(231:232) = '1 1 2'
!
spcgr_laue(233:234)  = '2/m'
spcgr_point(233:234) = '2'
!
spcgr_laue(235:238)  = '1 1 2/m'
spcgr_point(235:237) = '1 1 2'
spcgr_point(238:238) = '1 1 m'
!
spcgr_laue(239:240)  = '2/m'
spcgr_point(239:240) = 'm'
!
spcgr_laue(241:243)  = '1 1 2/m'
spcgr_point(241:243) = '1 1 m'
!
spcgr_laue(244:245)  = '2/m'
spcgr_point(244:245) = 'm'
!
spcgr_laue(246:248)  = '1 1 2/m'
spcgr_point(246:248) = '1 1 m'
!
spcgr_laue(249:250)  = '2/m'
spcgr_point(249:250) = 'm'
!
spcgr_laue(251:255)  = '1 1 2/m'
spcgr_point(251:253) = '1 1 m'
spcgr_point(254:255) = '1 1 2/m'
!
spcgr_laue(256:257)  = '2/m'
spcgr_point(256:257) = '2/m'
!
spcgr_laue(258:260)  = '1 1 2/m'
spcgr_point(258:260) = '1 1 2/m'
!
spcgr_laue(261:262)  = '2/m'
spcgr_point(261:262) = '2/m'
!
spcgr_laue(263:265)  = '1 1 2/m'
spcgr_point(263:265) = '1 1 2/m'
!
spcgr_laue(266:267)  = '2/m'
spcgr_point(266:267) = '2/m'
!
spcgr_laue(268:270)  = '1 1 2/m'
spcgr_point(268:270) = '1 1 2/m'
!
spcgr_laue(271:272)  = '2/m'
spcgr_point(271:272) = '2/m'
!
spcgr_laue(273:275)  = '1 1 2/m'
spcgr_point(273:275) = '1 1 2/m'
!!
spcgr_laue(276:280)  = 'm m m'
spcgr_point(276:280) = 'm m m'
!!
spcgr_laue(281:283)  = '4/m'
spcgr_point(281:283) = '4/m'
!
spcgr_laue(284:293)  = '4/m m m'
spcgr_point(284:293) = '4/m m m'
!!
spcgr_laue(294:295)  = '-3'
spcgr_point(294:294) = '3'
!
spcgr_point(295:295) = '-3'
!
spcgr_laue(296:300)  = '-3 m'
spcgr_point(296:296) = '3 2'
!
spcgr_point(297:298) = '3 m'
!
spcgr_point(299:300) = '-3 m'
!
!!
spcgr_laue(301:302)  = 'm -3'
spcgr_point(301:302) = 'm -3'
!
spcgr_laue(303:306)  = 'm -3 m'
spcgr_point(303:306) = 'm -3 m'
!
!!
spcgr_laue(307:315)  = 'm m m'
spcgr_point(307:309) = 'm m m'
!
spcgr_point(310:311) = 'm m 2'
!
spcgr_point(312:315) = 'm m m'
!!
! Build a control group of all point groups
spcgr_point_t          = ' '
spcgr_point_t(:)(1:15) = spcgr_name(:)(2:16)
do i=1, SPCGR_MAX
   do l=1,16
      if(spcgr_point_t(i)(l:l) == 'a') spcgr_point_t(i)(l:l) = 'm'
      if(spcgr_point_t(i)(l:l) == 'b') spcgr_point_t(i)(l:l) = 'm'
      if(spcgr_point_t(i)(l:l) == 'c') spcgr_point_t(i)(l:l) = 'm'
      if(spcgr_point_t(i)(l:l) == 'e') spcgr_point_t(i)(l:l) = 'm'
      if(spcgr_point_t(i)(l:l) == 'n') spcgr_point_t(i)(l:l) = 'm'
      if(spcgr_point_t(i)(l:l) == 'd') spcgr_point_t(i)(l:l) = 'm'
   enddo
   loop_21: do
      k = index(spcgr_point_t(i),'21')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 2'
      if(k==0) exit loop_21
   enddo loop_21
   loop_31: do
      k = index(spcgr_point_t(i),'31')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 3'
      if(k==0) exit loop_31
   enddo loop_31
   loop_32: do
      k = index(spcgr_point_t(i),'32')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 3'
      if(k==0) exit loop_32
   enddo loop_32
   loop_41: do
      k = index(spcgr_point_t(i),'41')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 4'
      if(k==0) exit loop_41
   enddo loop_41
   loop_42: do
      k = index(spcgr_point_t(i),'42')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 4'
      if(k==0) exit loop_42
   enddo loop_42
   loop_43: do
      k = index(spcgr_point_t(i),'43')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 4'
      if(k==0) exit loop_43
   enddo loop_43
   loop_61: do
      k = index(spcgr_point_t(i),'61')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 6'
      if(k==0) exit loop_61
   enddo loop_61
   loop_62: do
      k = index(spcgr_point_t(i),'62')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 6'
      if(k==0) exit loop_62
   enddo loop_62
   loop_63: do
      k = index(spcgr_point_t(i),'63')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 6'
      if(k==0) exit loop_63
   enddo loop_63
   loop_64: do
      k = index(spcgr_point_t(i),'64')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 6'
      if(k==0) exit loop_64
   enddo loop_64
   loop_65: do
      k = index(spcgr_point_t(i),'65')
      if(k>0) spcgr_point_t(i)(k:k+1) = ' 6'
      if(k==0) exit loop_65
   enddo loop_65
enddo
do i=143, 167
   k = len_trim(spcgr_point_t(i))
   if(spcgr_point_t(i)(k:k) == '1') spcgr_point_t(i)(k:k) = ' '
enddo
do_231: do i=231, 275
   k = len_trim(spcgr_point_t(i))
   if(k>4) then
      if(spcgr_point_t(i)(k-4:k) == '1 2 1') then
         spcgr_point_t(i)(k-4:k) = '2    '
         cycle do_231
      endif
!   k = len_trim(spcgr_point_t(i))
      if(spcgr_point_t(i)(k-4:k) == '1 m 1') then
         spcgr_point_t(i)(k-4:k) = 'm    '
         cycle do_231
      endif
!  k = len_trim(spcgr_point_t(i))
   elseif(k>6) then
      if(spcgr_point_t(i)(k-6:k) == '1 2/m 1') then
         spcgr_point_t(i)(k-6:k) = '2/m    '
         cycle do_231
      endif
!  k = len_trim(spcgr_point_t(i))
   elseif(k>7) then
      if(spcgr_point_t(i)(k-7:k) == '1  2/m 1') then
         spcgr_point_t(i)(k-7:k) = '2/m     '
         cycle do_231
      endif
   endif
enddo do_231
!   
do i=1, SPCGR_MAX
   l = len_trim(spcgr_point_t(i))
   call rem_dbl_bl(spcgr_point_t(i), l)
   call rem_leading_bl(spcgr_point_t(i), l)
   l = len_trim(spcgr_point(i))
   call rem_dbl_bl(spcgr_point(i), l)
   call rem_leading_bl(spcgr_point(i), l)
enddo 
!
!if(j==33) then
!lall = .true.
!do i=1, 2
!  write(*,'(i4,a17, 2x, a6, 2x, a6, l2,1x,a7)') i, spcgr_name(i),spcgr_point(i)(1:len_trim(spcgr_point(i))), spcgr_laue(i)(1:len_trim(spcgr_laue(i))), spcgr_point(i)==spcgr_point_t(i), spcgr_point_t(i)
!  lall = lall .and. spcgr_point(i)==spcgr_point_t(i)
!enddo
!write(*,*) ' ALL CORRECT ? ', lall
!write(*,*)
!read(*,*) i
!do i=3, 15
!  write(*,'(i4,a17, 2x, a6, 2x, a6, l2,1x,a7)') i, spcgr_name(i),spcgr_point(i)(1:len_trim(spcgr_point(i))), spcgr_laue(i)(1:len_trim(spcgr_laue(i))), spcgr_point(i)==spcgr_point_t(i), spcgr_point_t(i)
!  lall = lall .and. spcgr_point(i)==spcgr_point_t(i)
!enddo
!write(*,*) ' ALL CORRECT ? ', lall
!write(*,*)
!read(*,*) i
!do i=16, 74
!  write(*,'(i4,a17, 2x, a6, 2x, a6, l2,1x,a7)') i, spcgr_name(i),spcgr_point(i)(1:len_trim(spcgr_point(i))), spcgr_laue(i)(1:len_trim(spcgr_laue(i))), spcgr_point(i)==spcgr_point_t(i), spcgr_point_t(i)
!  lall = lall .and. spcgr_point(i)==spcgr_point_t(i)
!enddo
!write(*,*) ' ALL CORRECT ? ', lall
!write(*,*)
!read(*,*) i
!do i=75, 142
!  write(*,'(i4,a17, 2x, a6, 2x, a6, l2,1x,a7)') i, spcgr_name(i),spcgr_point(i)(1:len_trim(spcgr_point(i))), spcgr_laue(i)(1:len_trim(spcgr_laue(i))) , spcgr_point(i)==spcgr_point_t(i), spcgr_point_t(i)
!  lall = lall .and. spcgr_point(i)==spcgr_point_t(i)
!enddo
!write(*,*) ' ALL CORRECT ? ', lall
!write(*,*)
!read(*,*) i
!do i=143,167
!  write(*,'(i4,a17, 2x, a6, 2x, a6, l2,1x,a7)') i, spcgr_name(i),spcgr_point(i)(1:len_trim(spcgr_point(i))), spcgr_laue(i)(1:len_trim(spcgr_laue(i))), spcgr_point(i)==spcgr_point_t(i), spcgr_point_t(i)
!  lall = lall .and. spcgr_point(i)==spcgr_point_t(i)
!enddo
!write(*,*) ' ALL CORRECT ? ', lall
!write(*,*)
!read(*,*) i
!do i=168, 194
!  write(*,'(i4,a17, 2x, a6, 2x, a6, l2,1x,a7)') i, spcgr_name(i),spcgr_point(i)(1:len_trim(spcgr_point(i))), spcgr_laue(i)(1:len_trim(spcgr_laue(i))), spcgr_point(i)==spcgr_point_t(i), spcgr_point_t(i)
!  lall = lall .and. spcgr_point(i)==spcgr_point_t(i)
!enddo
!!read(*,*) i
!write(*,*) ' ALL CORRECT ? ', lall
!write(*,*)
!read(*,*) i
!do i=195, 230
!  write(*,'(i4,a17, 2x, a6, 2x, a6, l2,1x,a7)') i, spcgr_name(i),spcgr_point(i)(1:len_trim(spcgr_point(i))), spcgr_laue(i)(1:len_trim(spcgr_laue(i))), spcgr_point(i)==spcgr_point_t(i), spcgr_point_t(i)
!  lall = lall .and. spcgr_point(i)==spcgr_point_t(i)
!enddo
!write(*,*) ' ALL CORRECT ? ', lall
!write(*,*)
!read(*,*) i
!do i=231,315
!  write(*,'(i4,a17, 2x, a6, 2x, a6, l2,1x,a7)') i, spcgr_name(i),spcgr_point(i)(1:len_trim(spcgr_point(i))), spcgr_laue(i)(1:len_trim(spcgr_laue(i))), spcgr_point(i)==spcgr_point_t(i), spcgr_point_t(i)
!  lall = lall .and. spcgr_point(i)==spcgr_point_t(i)
!enddo
!write(*,*) ' ALL CORRECT ? ', lall
!j = 2
!endif
                                                                        
END SUBROUTINE spcgr_setup                    
END MODULE spcgr_setup_mod
