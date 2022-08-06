
SUBROUTINE FFT_1D(N,F,SIGN)
 
                                                                     IMPLICIT NONE     
                                                                     
  INCLUDE 'FFTW3.F'  
  
!--------------------------------------------------------------------------------------------------

  
   INTEGER(KIND=4)::SIGN,N
   
   COMPLEX(KIND=8)::F(N)
   
   
   COMPLEX(KIND=8),ALLOCATABLE::FO(:)
   
   INTEGER(KIND=8)::FWD,BWD,STATUS
                                                                           
!--------------------------------------------------------------------------------------------------
                                                                                                                                                                        
                                                                                      
 ALLOCATE(FO(N)); FO(:)=0.   


 IF(SIGN==1) THEN
    
      FWD=0; STATUS=0

      CALL DFFTW_PLAN_DFT_1D(FWD,N,F,FO,FFTW_FORWARD,FFTW_ESTIMATE)  
       
      CALL DFFTW_EXECUTE(FWD,F,FO)
        
      CALL DFFTW_DESTROY_PLAN(FWD) 
      
      F(:)=FO(:)
      
 END IF

 IF(SIGN==-1) THEN
        
        BWD=0; STATUS=0
        
        CALL DFFTW_PLAN_DFT_1D(BWD,N,F,FO,FFTW_BACKWARD,FFTW_ESTIMATE)
  
        CALL DFFTW_EXECUTE_DFT(BWD,F,FO) 
        
        CALL DFFTW_DESTROY_PLAN(BWD)
        
        F(:)=FO(:)/N
      
 END IF

 
 DEALLOCATE(FO)       
 
!--------------------------------------------------------------------------------------------------                                                                                      
                                                                                            RETURN
 END SUBROUTINE                                                                                      
!--------------------------------------------------------------------------------------------------

