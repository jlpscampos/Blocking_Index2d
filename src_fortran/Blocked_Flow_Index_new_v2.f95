MODULE BLOCKED_FLOW_Index_new
  !_________________________________________________________________________________________________________!
  !Purpose : In this module contains a set of subroutines and functions to compute blocking index using the !
  !          methodology developed by Tibaldi & Molteni 1990 and by Pelly & Hoskins 2003                    !
  !          * Here the blocking index are computed in 2 dimensions!                                        !
  !                                                                                                         !
  !Record of Revisions:                                                                                     !
  !     DATE                      NAME                            DESCRIPTION                               !
  !     27/08/2013                J.Leandro P.S.Campos            Original code                             !
  !     06/08/2015                J.Leandro P.S.Campos            Revision of some subroutines              !
  !     07/08/2015                J.Leandro P.S.Campos            Subroutine TM_Index added                 !
  !     10/08/2017                J.Leandro P.S.Campos            Main program changed now read netcdf files!
  !     10/08/2017                J.Leandro P.S.Campos            Added lines in blocking Subroutine        !
  !     26/06/2018                J.Leandro P.S.Campos            Added a new subroutine                    !
  !_________________________________________________________________________________________________________!
  
  IMPLICIT NONE

  REAL,PARAMETER,PRIVATE :: pi = 3.1416                 !number pi
  REAL,PARAMETER,PRIVATE :: a = 6.37E6                  !Earth's radius

  REAL,PARAMETER,PRIVATE :: dphi = 15.0                 !latitudinal blocking extension
  REAL,PARAMETER,PRIVATE :: delta= 4.0                  !Latitudinal allowed variation for blockings  
  REAL,PARAMETER,PRIVATE :: l = 7.5                     !maximun gap between two blocked longitudes
  REAL,PARAMETER,PRIVATE :: w = 15.0                    !Minimum blocking width
  REAL,PARAMETER,PRIVATE :: deltalambda = 5.0           !Smooth factor
  REAL,PARAMETER,PUBLIC  :: p = 10.0                    !Overlap parameter
  
  INTEGER,PARAMETER,PUBLIC :: SINGLE = 4
  
  !_________________________________________________________________________________________________________!
  ! Blocking Variables (Unused in the main routine - 20/07/2018) Some fortran compilers does nor support this
  TYPE :: blocking_PH                                   !blocking variable for PV-Theta method
    REAL :: index                          !-> binary index 
    REAL :: theta                          !-> input variable
    REAL :: dif                            !-> diference
  END TYPE

  TYPE :: blocking_TM                                   !blocking variable for Tibaldi & Molteni's method
    REAL :: index                          !-> binary index
    REAL :: hgt                            !-> input variable
    REAL :: ghgn                           !-> gradient northward central latitude
    REAL :: ghgs                           !-> gradient southward central latitude
  END TYPE
  !________________________________________________________________________________________________________!
  
  CONTAINS
  
    SUBROUTINE PH_Index ( ind, B, tht, xdim, ydim, res )
      !Computing the blocking index as in Pelly and Hoslkins (2003)
      IMPLICIT NONE

      INTEGER,INTENT(in) :: xdim, ydim
      REAL,DIMENSION(xdim,ydim),INTENT(in)  :: tht
      REAL,DIMENSION(xdim,ydim),INTENT(out) :: B, ind
      REAL,INTENT(in) :: res

      INTEGER :: dpi, del
      INTEGER :: i, j, k

      ! Computing the difference between the mean southward and equatorward the central latitude j
      dpi = INT(dphi/res)
      DO i = 1, xdim
        DO j = dpi+1, ydim-(dpi+1)
          B(i,j) = (1./REAL(dpi))*( SUM( tht(i,j-dpi:j) ) - SUM( tht(i,j:j+dpi) ) )
        END DO
      END DO
        
      ! Computing the blocking index if B > 0 -> blocking in the SH
      !                                 B < 0 -> blocking in the NH 
      del = INT(delta/res) ; ind = 0.
      DO i = 1, xdim
        DO j = dpi+del, ydim-(dpi+del)
          
          DO k = j-del, j+del
            IF( j < INT(ydim/2) .and. B(i,k) > 0. ) ind(i,j) = 1
            IF( j > INT(ydim/2) .and. B(i,k) < 0. ) ind(i,j) = 1
          END DO
          
        END DO
      END DO

    END SUBROUTINE PH_Index

  !________________________________________________________________________________________________________!

    SUBROUTINE TM_INDEX( ind, ghgn, ghgs, hgt, xdim, ydim, res )
      ! Blocking index based on 500hPa Geopotential Height 
      IMPLICIT NONE

      INTEGER,INTENT(in) :: xdim, ydim
      REAL,DIMENSION(xdim,ydim),INTENT(inout) :: ghgn, ghgs, ind
      REAL,DIMENSION(xdim,ydim),INTENT(in) :: hgt
      REAL,INTENT(in) :: res

      INTEGER :: dpi, del
      INTEGER :: i, j, k

      ! Computing the GHGN and GHGS parameters
      dpi = INT(dphi/res)
      DO i = 1, xdim
        DO j = dpi+1, ydim-(dpi+1)

          IF( j < INT(ydim/2) )THEN      ! Southern Hemisphere
            GHGN(i,j) = -1*( hgt(i,j)-hgt(i,j+dpi) )/( ANG(j) - ANG(j+dpi) )
            GHGS(i,j) = -1*( hgt(i,j-dpi)-hgt(i,j) )/( ANG(j-dpi) - ANG(j) )
          ELSE IF( j> INT(ydim/2) )THEN  ! Northern Hemisphere
            GHGN(i,j) = ( hgt(i,j+dpi)-hgt(i,j) )/( ANG(j+dpi) - ANG(j) )
            GHGS(i,j) = ( hgt(i,j)-hgt(i,j-dpi) )/( ANG(j) - ANG(j-dpi) )
          END IF
          
        END DO
      END DO

      ! Computing the blocking index
      del = INT(delta/res) ; ind = 0.
      DO i = 1, xdim
        DO j = dpi+del, ydim-(dpi+del)

          DO k = j-del, j+del
            IF( j<INT(ydim/2) .and. (GHGS(i,j)<-10. .and. GHGN(i,j)>0.) ) ind(i,j) = 1.
            IF( j>INT(ydim/2) .and. (GHGS(i,j)>0. .and. GHGN(i,j)<-10.) ) ind(i,j) = 1.
          END DO
          
        END DO
      END DO

      END SUBROUTINE TM_INDEX


      REAL FUNCTION ANG( val )
        IMPLICIT NONE

        INTEGER,INTENT(iN) :: val
        INTEGER :: res = 1.5

        ANG = -90 + res*(val-1)
        RETURN
      END FUNCTION ANG

  !________________________________________________________________________________________________________!

    SUBROUTINE GLBI ( ind, xdim, ydim, res )
      IMPLICIT NONE

      INTEGER,INTENT(in) :: xdim, ydim
      REAL,DIMENSION(xdim,ydim),INTENT(inout) :: ind
      REAL,INTENT(in) :: res

      INTEGER :: i, j, k, par, ww

      ! Fill all the gaps of l degrees between contiguous blocked longitudes
      par = INT(l/res)
      DO j = 1, ydim
        DO i = 1, xdim-par
          
          IF( INT(ind(i,j)) == 1 )THEN
            DO k = i+par, i, -1
              IF( INT(ind(k,j)) == 1 )THEN
                ind(i:k,j ) = 1. ; EXIT  
              END IF
            END DO
          END IF
          
        END DO
      END DO

      !Fill the gaps between +1 or -1 latitude grids (added in 17/08/2016)
      DO i = 1, xdim
        DO j = 2, ydim-1

          IF( ind(i,j-1) == 1 .and. ind(i,j+1) == 1 )THEN
            ind(i,j) = 1.
          END IF

        END DO
      END DO

      ww = INT(w/res)
      ! Group all blocked latitudes with minimum width of W degrees
      DO j = 1, ydim
        DO i = 2, xdim
          IF( ind(i,j) > 0. ) ind(i,j) = ind(i,j) + ind(i-1,j)
        END DO
        DO i = xdim-1, 1, -1
          IF( ind(i,j) > 0 .and. ind(i+1,j) > 0. ) ind(i,j) = ind(i+1,j)
        END DO
      END DO
      WHERE( ind < REAL(ww) ) ind = 0.
      WHERE( ind >=REAL(ww) ) ind = 1.
        
    END SUBROUTINE GLBI
    !________________________________________________________________________________________________________!


    !________________________________________________________________________________________________________!

    SUBROUTINE BLOCKING (blk, ind, xdim, ydim, tdim, res, persistence )
      IMPLICIT NONE

      INTEGER,INTENT(in) :: xdim, ydim, tdim, persistence
      REAL,INTENT(in) :: res
      REAL,DIMENSION(xdim,ydim,tdim),INTENT(in) :: ind
      REAL,DIMENSION(xdim,ydim,tdim),INTENT(inout) :: blk
      REAL,DIMENSION(xdim,ydim,tdim) :: xind
      
      INTEGER :: i, j, t, pp
      REAL :: val

      ! If day t-1 and t+1 are blocked and day t is not blocked, day t is also blocked. (Added 10/08/2017)
      xind = ind 
      DO t = 2, tdim-1
        DO i = 1,xdim
          DO j = 1,ydim
            IF( xind(i,j,t)<1. .and. ( xind(i,j,t-1)>0. .and. xind(i,j,t+1)>0. ) ) xind(i,j,t) = xind(i,j,t-1);
          END DO  
		END DO
      END DO

      ! Find all group of blocked longitudes one day prior and after the blocked flow that overlap at least p degrees.
      !xind = ind ;
      pp = INT(p/res)
      DO t = 2, tdim
        DO j = 1, ydim
          DO i = 1, xdim
            IF( xind(i,j,t)>0. .and. xind(i,j,t-1)>0. ) xind(i,j,t)=xind(i,j,t)+xind(i,j,t-1)
          END DO
        END DO
      END DO

      DO t = tdim-1, 1, -1
        DO j = 1, ydim
          DO i = 1, xdim
            IF( xind(i,j,t) > 0. .and. xind(i,j,t+1) > 0. ) xind(i,j,t) = xind(i,j,t+1)
          END DO
        END DO
      END DO

      DO t = 1, tdim-1
        
        DO j = 1, ydim-pp
          DO i = 1, xdim-pp
            val = MAXVAL( xind(i:i+pp, j:j+pp, t) ) 
            WHERE( xind(i:i+pp, j:j+pp, t)>0. ) xind(i:i+pp, j:j+pp, t) = val
          END DO
        END DO

        DO j = ydim, pp+1, -1
          DO i = xdim, pp+1, -1
            val = MAXVAL( xind(i-pp:i, j-pp:j, t) ) 
            WHERE( xind(i-pp:i, j-pp:j, t)>0. ) xind(i-pp:i, j-pp:j, t) = val
          END DO
        END DO

      END DO

      ! All blocked flow are blocking if last at least persistence days
      blk = 0.
      WHERE( xind >= REAL(persistence) ) blk = 1.
    
    END SUBROUTINE BLOCKING

    
    !________________________________________________________________________________________________________________!


    SUBROUTINE SMOOTH(input,output, xdim, ydim, tdim, res)
      ! This subroutine smooth the data in deltalambda degrees 
      ! Use this before perform the blocked flow index for the
      ! Potential Temperature field
      
      INTEGER,INTENT(in) :: xdim, ydim,tdim
      REAL,INTENT(in) :: res
      REAL,DIMENSION(xdim,ydim,tdim) :: input
      REAL,DIMENSION(xdim,ydim,tdim) :: output
      INTEGER :: i,j,t,smth
 
      smth = INT( (deltalambda/res)/2. )

      DO t = 1,tdim
        DO j = 1,ydim
          DO i = 1,xdim

            IF( i > smth .or. i < (xdim-smth) )THEN
              output(i,j,t) = SUM(input(i-smth:i+smth,j,t))/REAL(2*smth)
            ELSE IF( i < smth )THEN
              output(i,j,t) = (SUM(input(smth-i:i+smth,j,t))+SUM(input(xdim-smth:xdim,j,t)))/REAL(2*smth)
            ELSE IF( i > xdim-smth )THEN
              output(i,j,t) = (SUM(input(i-smth:xdim,j,t))+SUM(input(1:xdim-i,j,t)))/REAL(2*smth)
            ELSE IF( i == smth )THEN
              output(i,j,t) = ( SUM(input(1:smth,j,t))+input(xdim,j,t))/REAL(2*smth)
            ELSE IF( i == xdim-smth )THEN
              output(i,j,t) = ( SUM(input(xdim-smth:smth,j,t))+input(1,j,t))/REAL(2*smth)
            END IF
            
          END DO
        END DO
      END DO
    
    END SUBROUTINE SMOOTH

    !______________________________________________________________________________________!
    SUBROUTINE DAILY_MEANS(input, output, xdim, ydim, tdim, time1day)

       INTEGER,INTENT(in) :: xdim, ydim, tdim
       INTEGER,INTENT(in) :: time1day
       REAL,DIMENSION(xdim,ydim,tdim),INTENT(in) :: input
       REAL,DIMENSION(xdim,ydim,INT(tdim/4)),INTENT(out):: output
       REAL,DIMENSION(xdim,ydim) :: med
       INTEGER :: t, k, i

       i = 1;
       DO t = 1,tdim,time1day
         med = 0
         DO k = t,t+time1day-1
           med = input(:,:,k)+med
         END DO
         output(:,:,i) = med/REAL(time1day)
         i = i+1
       END DO

    END SUBROUTINE DAILY_MEANS
    

END MODULE BLOCKED_FLOW_Index_new    