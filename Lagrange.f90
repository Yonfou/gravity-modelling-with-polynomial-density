
SUBROUTINE Lagrange(N,N_MAX,Xnd,CSum)
!-------------------------------------------------------------------------------------------  
  IMPLICIT NONE
  
INTEGER(KIND=4) :: N,N_MAX
REAL   (KIND=8) :: Xnd(N_MAX),CSum

INTEGER(KIND=4) :: II,JJ,rank,Nleft
REAL   (KIND=8) :: Fac,AMUL
INTEGER(KIND=4),ALLOCATABLE :: T(:)
!------------------------------------------------------------------------------------------- 


       !---------------------------------------------
       IF(N<N_MAX) THEN  
           !-----------------------------------------
           Nleft=N_MAX-N
           !-----------------------------------------
           ALLOCATE(T(Nleft))
           !-----------------------------------------
           Fac=1.D0
           !-----------------------------------------
           DO II=N+1,N_MAX
              !-------------------------------------- 
              Fac=Fac*II
              !-------------------------------------- 
           END DO
           !-----------------------------------------
           DO II=1,Nleft
              !--------------------------------------  
              Fac=Fac/II 
              !--------------------------------------
           END DO
           !-----------------------------------------
           rank=Fac
           !-----------------------------------------
           CSum=0.D0
           !-----------------------------------------
           DO II=0,rank-1
              !-------------------------------------- 
              T(:)=0; CALL ksubset_colex_unrank(II,Nleft,N_MAX,T)
              !-------------------------------------- 
              AMUL=1.D0
              !-------------------------------------- 
              DO JJ=1,Nleft
                 !----------------------------------- 
                 AMUL=AMUL*Xnd(T(JJ))
                 !----------------------------------- 
              END DO
              !-------------------------------------- 
              CSum=CSum+AMUL
              !-------------------------------------- 
           END DO
           !-----------------------------------------
           CSum=CSum*(-1.)**Nleft
           !-----------------------------------------
           DEALLOCATE(T)
           !-----------------------------------------
       ELSE
           !-----------------------------------------
           CSum=1.D0
           !-----------------------------------------
       END IF
       !---------------------------------------------


!------------------------------------------------------------------------------------------- 
                                                                                     RETURN
    END SUBROUTINE Lagrange
!------------------------------------------------------------------------------------------- 

    
    
    
    
    
    
subroutine ksubset_colex_unrank ( rank, k, n, t )

!*****************************************************************************
!
!! KSUBSET_COLEX_UNRANK computes the K subset of given colex rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer RANK, the rank of the K subset.
!
!    Input, integer K, the number of elements each K subset must
!    have.  0 <= K <= N.
!
!    Input, integer N, the number of elements in the master set.
!    N must be positive.
!
!    Output, integer T(K), describes the K subset of the given
!    rank.  T(I) is the I-th element.  The elements must be listed in
!    DESCENDING order.
!
  implicit none

  integer k

  integer i
  integer i4_choose
  integer n
  integer nksub
  integer rank
  integer rank_copy
  integer t(k)
  integer x
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop 1
  end if

  if ( k == 0 ) then
    return
  endif

  if ( k < 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input K is illegal.'
    stop 1
  end if

  call ksubset_enum ( k, n, nksub )

  if ( rank < 0 .or. nksub < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUBSET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop 1
  end if
!
  rank_copy = rank

  x = n

  do i = 1, k

    do while ( rank_copy < i4_choose ( x, k + 1 - i ) )
      x = x - 1
    end do

    t(i) = x + 1
    rank_copy = rank_copy - i4_choose ( x, k + 1 - i )

  end do

  return
    end
subroutine ksubset_enum ( k, n, nksub )

!*****************************************************************************80
!
!! KSUBSET_ENUM enumerates the K element subsets of an N set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer K, the number of elements each K subset must
!    have. 0 <= K <= N.
!
!    Input, integer N, the number of elements in the master set.
!    0 <= N.
!
!    Output, integer NKSUB, the number of distinct elements.
!
  implicit none

  integer i4_choose
  integer k
  integer n
  integer nksub

  nksub = i4_choose ( n, k )

  return
    end
function i4_choose ( n, k )

!*****************************************************************************80
!
!! I4_CHOOSE computes the binomial coefficient C(N,K) as an I4.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer N, K, are the values of N and K.
!
!    Output, integer I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer i
  integer i4_choose
  integer k
  integer mn
  integer mx
  integer n
  integer value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end