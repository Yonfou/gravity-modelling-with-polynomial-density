
SUBROUTINE ReadIn()
!-------------------------------------------------------------------------------------------
                                                                             USE Source
                                                                             USE ObserPoint
                                                                             USE GravityField
                                                                             USE Const
                                                                             USE Gauss
  IMPLICIT NONE
  
  INTEGER(KIND=4) :: II,JJ,KK,NN,NT
  REAL   (KIND=8) :: bb,ro_0,ro_1,ro_2,CA,CB,CC,SpC(3),SpR(2)
  CHARACTER(LEN=100):: Modelname
  REAL   (KIND=8),ALLOCATABLE :: T(:)
!-------------------------------------------------------------------------------------------

 
  
                !-------------------------------------------------------------                
                !                         CONSTANTS                           
                !-------------------------------------------------------------
                G0=6.674D-11                                      
                !-------------------------------------------------------------
                PI=3.141592653589793238462643383279D0             
                !-------------------------------------------------------------
                
                
                
                
                !-------------------------------------------------------------                
                !               INPUT THE BASIC PARAMETERS                             
                !-------------------------------------------------------------
                OPEN(1,FILE='Para.txt',STATUS='UNKNOWN')
                    !---------------------------------------------------------
                    READ(1,*);        READ(1,*);        READ(1,*)
                    !---------------------------------------------------------
                    READ(1,*)  Modelname
                    READ(1,*)  NG1,NG2
                    READ(1,*)  Norm
                    READ(1,*)  alpha
                    READ(1,*)  NSL
                    READ(1,*)  NSF
                    READ(1,*)  NPL
                    READ(1,*)  NPF
                    READ(1,*)  NPR
                    READ(1,*)  REarth;      REarth=REarth*1.D3
                    !---------------------------------------------------------
                    READ(1,*)  file_sl
                    READ(1,*)  file_sf
                    READ(1,*)  file_srn
                    READ(1,*)  file_sr
                    READ(1,*)  file_dgr
                    READ(1,*)  file_obs
                    READ(1,*)  file_rho
                    READ(1,*)  file_rst
                    !---------------------------------------------------------
                CLOSE(1)
                !-------------------------------------------------------------
                
                
                WRITE(*,*) 
                WRITE(*,*) 'Current model:  ',TRIM(ADJUSTL(Modelname))
                WRITE(*,*)
                
                

                !-------------------------------------------------------------                
                !           INPUT THE MESH GRID FOR THE SOURCE REGION                             
                !-------------------------------------------------------------
                OPEN(1,FILE=TRIM(ADJUSTL(file_sl)),STATUS='UNKNOWN')
                    !---------------------------------------------------------
                    ALLOCATE(SL(NSL));    SL(:)=0.D0
                    !---------------------------------------------------------
                    DO II=1,NSL
                       !------------------------------------------------------ 
                       READ(1,*) SL(II)
                       !------------------------------------------------------ 
                    END DO
                    !---------------------------------------------------------
                    dSL=SL(2)-SL(1)
                    !---------------------------------------------------------
                CLOSE(1)
                !-------------------------------------------------------------
                OPEN(1,FILE=TRIM(ADJUSTL(file_sf)),STATUS='UNKNOWN')
                    !---------------------------------------------------------
                    ALLOCATE(SF(NSF),dSF(NSF));    SF(:)=0.D0;    dSF(:)=0.D0
                    !---------------------------------------------------------
                    DO II=1,NSF
                       !------------------------------------------------------ 
                       READ(1,*) SF(II),dSF(II)
                       !------------------------------------------------------ 
                    END DO
                    !---------------------------------------------------------
                CLOSE(1)
                !-------------------------------------------------------------
                OPEN(1,FILE=TRIM(ADJUSTL(file_srn)),STATUS='UNKNOWN')
                    !---------------------------------------------------------
                    ALLOCATE(NR(NSF));    NR(:)=0
                    !---------------------------------------------------------
                    DO II=1,NSF
                       !------------------------------------------------------ 
                       READ(1,*) NR(II)
                       !------------------------------------------------------ 
                    END DO
                    !---------------------------------------------------------
                CLOSE(1)
                !-------------------------------------------------------------
                OPEN(1,FILE=TRIM(ADJUSTL(file_sr)),STATUS='UNKNOWN')
                    !---------------------------------------------------------
                    NSR=MAXVAL(NR)
                    !---------------------------------------------------------
                    ALLOCATE(SR(NSR,NSF));   SR(:,:)=0.D0
                    !---------------------------------------------------------
                    ALLOCATE(dSR(NSR,NSF));  dSR(:,:)=0.D0
                    !---------------------------------------------------------
                    DO II=1,NSR
                    DO JJ=1,NSF
                       !------------------------------------------------------ 
                       READ(1,*) SR(II,JJ),dSR(II,JJ)
                       !------------------------------------------------------ 
                    END DO
                    ENDDO
                    !---------------------------------------------------------
                CLOSE(1)
                !-------------------------------------------------------------
                OPEN(1,FILE=TRIM(ADJUSTL(file_dgr)),STATUS='UNKNOWN')
                    !---------------------------------------------------------
                    ALLOCATE(NP(NSR,NSF));    NP(:,:)=0
                    !---------------------------------------------------------
                    DO II=1,NSR
                       !------------------------------------------------------ 
                       READ(1,*) NP(II,:)
                       !------------------------------------------------------ 
                    END DO
                    !---------------------------------------------------------
                CLOSE(1)
                !-------------------------------------------------------------
                
                
                
               
                
                !-------------------------------------------------------------                
                !              INPUT THE OBSERVATION POINTS                             
                !-------------------------------------------------------------
                ALLOCATE(PL(NPL),PF(NPF),PR(NPR))
                !-------------------------------------------------------------
                OPEN(1,FILE=TRIM(ADJUSTL(file_obs)),STATUS='UNKNOWN')
                    !---------------------------------------------------------
                    PL(:)=0.D0;            PF(:)=0.D0;            PR(:)=0.D0  
                    !---------------------------------------------------------
                    DO II=1,NPR
                    DO JJ=1,NPF    
                    DO KK=1,NPL
                       !------------------------------------------------------ 
                       READ(1,*) PR(II),PF(JJ),PL(KK)
                       !------------------------------------------------------
                    END DO
                    END DO
                    END DO
                    !---------------------------------------------------------
                CLOSE(1)
                !-------------------------------------------------------------
                
                
                
        
                !-------------------------------------------------------------                
                !                     ALLOCATE ARRAYS                             
                !-------------------------------------------------------------                
                CALL DeOrAllocate(1)
                !------------------------------------------------------------- 
                
                
                
                
                !-------------------------------------------------------------                
                !                     DENSITY MODEL                             
                !-------------------------------------------------------------                
                ALLOCATE(Ro(MAXVAL(NP)+1,NSR,NSF,NSL));  Ro(:,:,:,:)=0.D0
                !-------------------------------------------------------------
                OPEN(1,FILE=TRIM(ADJUSTL(file_rho)),STATUS='UNKNOWN')
                    !---------------------------------------------------------
                    DO NN=1,SIZE(Ro,1)
                    DO II=1,NSR
                    DO JJ=1,NSF    
                    DO KK=1,NSL
                       !------------------------------------------------------ 
                       READ(1,*) Ro(NN,II,JJ,KK)
                       !------------------------------------------------------
                    END DO
                    END DO
                    END DO
                    END DO
                    !---------------------------------------------------------
                CLOSE(1)
                !-------------------------------------------------------------
                
                 
                
                
                !-------------------------------------------------------------                
                !                      GAUSS NODES                             
                !------------------------------------------------------------- 
                CALL GaussNodes(NG1,AP1,TP1,NG2,AP2,TP2,APX); APX(:)=APX(:)*G0/4.D0     
                !------------------------------------------------------------- 
                
                
                
                
                !-------------------------------------------------------------
                !             ANALYTICAL SOLUTION FOR THE SHELL                         
                !-------------------------------------------------------------   
                NT=0;             ALLOCATE(T(0:NT));             T(0)=Ro(1,1,1,1)
                !-------------------------------------------------------------
                SpC(:)=0.D0;  SpR(1)=SR(1,1)-0.5*dSR(1,1);  SpR(2)=SR(NSR,1)+0.5*dSR(NSR,1)
                !-------------------------------------------------------------
                IF(TRIM(ADJUSTL(Modelname))=='Shell') CALL SperShell(1,NT,T,CA,CB,CC,SpC,SpR,PR,PF,PL,NPR,NPF,NPL,GV0,GA0,GG0)
                !-------------------------------------------------------------
                DEALLOCATE(T)
                !-------------------------------------------------------------
                
                
                
                !-------------------------------------------------------------
                !             ANALYTICAL SOLUTION FOR THE PREM MODEL                         
                !-------------------------------------------------------------
                IF(TRIM(ADJUSTL(Modelname))=='PREM') CALL PREM(PR,PF,PL,NPR,NPF,NPL,GV0,GA0,GG0)
                !-------------------------------------------------------------
                
                
                
                
!-------------------------------------------------------------------------------------------
                                                                                    RETURN
END SUBROUTINE
!===========================================================================================