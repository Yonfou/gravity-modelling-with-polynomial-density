
SUBROUTINE Coordinate(X0,X1,NX,DX,X)
!-------------------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(KIND=4) :: NX,I
  REAL   (KIND=8) :: X0,X1,DX,X(NX)
  
!-------------------------------------------------------------------------------------------  
  
                              
                !-------------------------------------------------------------
                X(:)=0.D0
                !-------------------------------------------------------------
                
                                                        
                !-------------------------------------------------------------
                DX=(X1-X0)/NX
                !------------------------------------------------------------- 
                DO I=1,NX
                   !---------------------------------------------------------- 
                   X(I)=X0+(I-1.D0/2.D0)*DX
                   !----------------------------------------------------------
                END DO
                !-------------------------------------------------------------
  
  
  
!-------------------------------------------------------------------------------------------
                                                                                    RETURN
END SUBROUTINE
!===========================================================================================
    
   
