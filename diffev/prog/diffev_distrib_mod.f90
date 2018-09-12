MODULE diffev_distrib_mod
!
IMPLICIT NONE
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE distrib_even(kid, indiv, nodes, pop_c, nindiv, inode, &
           kid_at_indiv, kid_at_node,                            &
           node_has_kids, node_max_kids, node_finished)
!
! Distrib kids evenly onto the available nodes
!
IMPLICIT NONE
!
INTEGER              , INTENT(INOUT) :: kid      ! We will work on this kid
INTEGER              , INTENT(INOUT) :: indiv    ! We will work on this indiv
INTEGER              , INTENT(IN   ) :: nodes    ! Actual number of nodes
INTEGER              , INTENT(IN   ) :: pop_c    ! Actual number of children
INTEGER              , INTENT(IN   ) :: nindiv   ! Actual number of maximum  indivs
INTEGER              , INTENT(IN   ) :: inode    ! Current node
INTEGER, DIMENSION(1:pop_c)        , INTENT(INOUT) :: kid_at_indiv
INTEGER, DIMENSION(1:pop_c)        , INTENT(INOUT) :: kid_at_node 
INTEGER, DIMENSION(1:nodes,0:pop_c), INTENT(INOUT) :: node_has_kids
INTEGER, DIMENSION(1:nodes)        , INTENT(IN)    :: node_max_kids
LOGICAL, DIMENSION(1:nodes)        , INTENT(INOUT) :: node_finished
!
INTEGER :: j,k
INTEGER :: nkid
INTEGER :: min_indiv
!
kid   = MINLOC(kid_at_indiv,1)                             ! This kid needs most work
indiv = kid_at_indiv(kid) + 1                              ! Next individual calculation to be done
IF(node_has_kids(inode,0) < node_max_kids(inode) .AND.  &  ! This node has room for kids
   kid_at_node(kid) == 0                         .AND.  &  ! Not at maximum kids
   indiv <= nindiv                                  ) THEN ! This kid is not yet on a node
!
   node_has_kids(inode,0) = node_has_kids(inode,0) + 1 
   node_has_kids(inode,node_has_kids(inode,0)) = kid       ! Record this kid at the current node
   kid_at_indiv(kid) = indiv
   kid_at_node (kid) = inode
ELSE
   nkid    = 0
   min_indiv = nindiv !+ 1
   find: DO j=1,node_has_kids(inode,0)                     ! Loop over all kids on this node
      k = node_has_kids(inode,j)                           ! Find a kid on this node that needs work
      IF(kid_at_indiv(k) < min_indiv) THEN                 ! This kid needs more work
         nkid = k
         min_indiv = kid_at_indiv(k)
      ENDIF
   ENDDO find
   IF(nkid>0) THEN                                         ! Found a kid that needs most work
      kid = nkid
      indiv = kid_at_indiv(kid) + 1                        ! Next individual calculation to be done
      kid_at_indiv(kid) = indiv
      kid_at_node (kid) = inode
   ELSE
      node_finished(inode) = .TRUE.
      kid   = 0
      indiv = 0
   ENDIF
ENDIF
!
END SUBROUTINE distrib_even
!
!*******************************************************************************
!
SUBROUTINE distrib_preve(kid, indiv, nodes, pop_c, nindiv, inode, &
           kid_at_indiv, node_has_kids, node_finished)
!
! Distrib kids evenly onto the available nodes
! Part for the second loop, kids are placed onto the previous node
!
IMPLICIT NONE
!
INTEGER              , INTENT(INOUT) :: kid      ! We will work on this kid
INTEGER              , INTENT(INOUT) :: indiv    ! We will work on this indiv
INTEGER              , INTENT(IN   ) :: nodes    ! Actual number of nodes
INTEGER              , INTENT(IN   ) :: pop_c    ! Actual number of children
INTEGER              , INTENT(IN   ) :: nindiv   ! Actual number of maximum  indivs
INTEGER              , INTENT(IN   ) :: inode    ! Current node
INTEGER, DIMENSION(1:pop_c)        , INTENT(INOUT) :: kid_at_indiv
INTEGER, DIMENSION(1:nodes,0:pop_c*nindiv), INTENT(INOUT) :: node_has_kids
LOGICAL, DIMENSION(1:nodes)        , INTENT(INOUT) :: node_finished
!
INTEGER :: j,k
INTEGER :: nkid
INTEGER :: min_indiv
!
kid   = 0
indiv = 0
nkid    = 0
min_indiv = nindiv !+ 1
find: DO j=1,node_has_kids(inode,0)                     ! Loop over all kids on this node
   k = node_has_kids(inode,j)                           ! Find a kid on this node that needs work
   IF(kid_at_indiv(k) < min_indiv) THEN                 ! This kid needs more work
      nkid = k
      min_indiv = kid_at_indiv(k)
   ENDIF
ENDDO find
!
IF(nkid>0) THEN                                         ! Found a kid that needs most work
   kid = nkid
   indiv = kid_at_indiv(kid) + 1                        ! Next individual calculation to be done
   kid_at_indiv(kid) = indiv
ELSE
   node_finished(inode) = .TRUE.
ENDIF
!
END SUBROUTINE distrib_preve
!
!*******************************************************************************
!
SUBROUTINE distrib_sequential(kid, indiv, nodes, pop_c, nindiv,  &
           numcalcs, inode, numsent,              &
           kid_at_indiv, kid_at_node, node_has_kids)
!
! Distributes all kids and indives sequentially onto the nodes and cores
!
IMPLICIT NONE
!
INTEGER              , INTENT(INOUT) :: kid      ! We will work on this kid
INTEGER              , INTENT(INOUT) :: indiv    ! We will work on this indiv
INTEGER              , INTENT(IN   ) :: nodes    ! Actual number of nodes
INTEGER              , INTENT(IN   ) :: pop_c    ! Actual number of children
INTEGER              , INTENT(IN   ) :: nindiv   ! Actual number of maximum  indivs
INTEGER              , INTENT(IN   ) :: numcalcs ! Total number of jobs needed
INTEGER              , INTENT(IN   ) :: inode    ! Current node
INTEGER              , INTENT(IN   ) :: numsent  ! Jobs sent out
INTEGER, DIMENSION(1:pop_c)        , INTENT(INOUT) :: kid_at_indiv
INTEGER, DIMENSION(1:pop_c)        , INTENT(INOUT) :: kid_at_node 
INTEGER, DIMENSION(1:nodes,0:pop_c*nindiv), INTENT(INOUT) :: node_has_kids
!
!
kid   = 0
indiv = 0
IF(numsent<numcalcs) THEN
   kid    = MOD( numsent,  pop_c) + 1
   indiv  =      numsent / pop_c  + 1
   kid_at_indiv(kid) = indiv
   kid_at_node (kid) = inode
   node_has_kids(inode,kid) = kid
   node_has_kids(inode,0  ) = kid
ENDIF
!
END SUBROUTINE distrib_sequential
!
END MODULE diffev_distrib_mod
