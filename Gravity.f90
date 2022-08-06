
SUBROUTINE Gravity()
!-------------------------------------------------------------------------------------------
                                                                             USE Source
                                                                             USE ObserPoint
                                                                             USE GravityField
                                                                             USE Const,ONLY:err
  IMPLICIT NONE
  
  INTEGER(KIND=4) :: I,J,II,JJ,KK,T0,T1,SGN0,SGN1,NP0,ISR(2)
  REAL   (KIND=8) :: SF0,SR0,dSF0,dSR0
  REAL   (KIND=8) :: RoC(NSL),F(NPL),FF(10,NPL)
  REAL   (KIND=8),ALLOCATABLE :: MV(:,:),MA(:,:,:),MG(:,:,:),WM(:,:)
!-------------------------------------------------------------------------------------------  

                
                
                     CALL SYSTEM_CLOCK(T0)
          
                                                                   
                    !-------------------------------------------------------------
                    GV(:,:,:)=0.D0;       GA(:,:,:,:)=0.D0;       GG(:,:,:,:)=0.D0
                    !-------------------------------------------------------------  
                    
                    
                    !-------------------------------------------------------------  
                    IF(PL(1)==SL(1)) THEN
                        
                        WRITE(*,*)
                        WRITE(*,*) 'Fast algorithm based on symmetric Toeplitz matrices...'
                        WRITE(*,*)
                       
                       
                        !---------------------------------------------------------   
                        ISR=MAXLOC(SR)
                        !---------------------------------------------------------   
                        DO I=1,NPR
                           !------------------------------------------------------  
                           WRITE(*,*) 'Observation surface: ', I,NPR 
                           !------------------------------------------------------  
                           IF(PR(I)>=MAXVAL(SR)+0.5*dSR(ISR(1),ISR(2))+1.D3) err=1.D-12
                           IF(PR(I)<MAXVAL(SR)+0.5*dSR(ISR(1),ISR(2))+1.D3)  err=1.D-16
                           !------------------------------------------------------ 
                        DO J=1,NPF
                           !------------------------------------------------------ 
                           FF(:,:)=0.D0
                           !------------------------------------------------------
                           DO JJ=1,NSF 
                               !---------------------------------------------------
                               SF0=SF(JJ);       dSF0=dSF(JJ)
                               !--------------------------------------------------- 
                           DO II=1,NR(JJ)
                               !--------------------------------------------------- 
                               NP0=NP(II,JJ);    SR0=SR(II,JJ);    dSR0=dSR(II,JJ)
                               !---------------------------------------------------
                               ALLOCATE(WM(0:NP0,0:NP0))
                               !---------------------------------------------------
                               ALLOCATE(MV(NSL,0:NP0),MA(NSL,0:NP0,3),MG(NSL,0:NP0,6))
                               !---------------------------------------------------
                               CALL Lagrange_Coeff(NP0,SR0,dSR0,WM)
                               !---------------------------------------------------  
                               CALL WeightCoef0(NP0,Norm,NPL,NSL,PR(I),PF(J),PL,SR0,SF0,SL,dSR0,dSF0,dSL,MV,MA,MG,WM)
                               !--------------------------------------------------- 
                               DO KK=0,NP0                  
                                  !------------------------------------------------ 
                                  RoC(:)=Ro(KK+1,II,JJ,:)
                                  !------------------------------------------------ 
                                  CALL  BCE0(NSL,MV(:,KK),RoC,F);                FF(1,:)=FF(1,:)+F(:)
                                  !------------------------------------------------ 
                                  CALL  BCE0(NSL,MA(:,KK,1),RoC,F);              FF(2,:)=FF(2,:)+F(:)
                                  CALL  BCE1(NSL,MA(:,KK,2),-MA(:,KK,2),RoC,F);  FF(3,:)=FF(3,:)+F(:)
                                  CALL  BCE0(NSL,MA(:,KK,3),RoC,F);              FF(4,:)=FF(4,:)+F(:)
                                  !------------------------------------------------ 
                                  CALL  BCE0(NSL,MG(:,KK,1),RoC,F);              FF(5,:)=FF(5,:)+F(:)
                                  CALL  BCE0(NSL,MG(:,KK,2),RoC,F);              FF(6,:)=FF(6,:)+F(:)
                                  CALL  BCE0(NSL,MG(:,KK,3),RoC,F);              FF(7,:)=FF(7,:)+F(:)
                                  CALL  BCE1(NSL,MG(:,KK,4),-MG(:,KK,4),RoC,F);  FF(8,:)=FF(8,:)+F(:)
                                  CALL  BCE0(NSL,MG(:,KK,5),RoC,F);              FF(9,:)=FF(9,:)+F(:)
                                  CALL  BCE1(NSL,MG(:,KK,6),-MG(:,KK,6),RoC,F);  FF(10,:)=FF(10,:)+F(:)
                                  !------------------------------------------------ 
                               END DO
                               !--------------------------------------------------- 
                               DEALLOCATE(WM,MV,MA,MG)
                               !--------------------------------------------------- 
                           END DO
                           !------------------------------------------------------- 
                           END DO
                           !------------------------------------------------------- 
                           GV(I,J,:)=FF(1,:)
                           GA(:,I,J,:)=FF(2:4,:)
                           GG(:,I,J,:)=FF(5:10,:)
                           !------------------------------------------------------- 
                        END DO
                        END DO
                        !---------------------------------------------------------   
                        
                    ELSE 
                        
                        WRITE(*,*)
                        WRITE(*,*) 'Fast algorithm based on non-symmetric Toeplitz matrices...'
                        WRITE(*,*)
                        
                        !---------------------------------------------------------   
                        ISR=MAXLOC(SR)
                        !---------------------------------------------------------   
                        DO I=1,NPR
                           !------------------------------------------------------  
                           WRITE(*,*) 'Observation surface: ', I,NPR 
                           !------------------------------------------------------  
                           IF(PR(I)>=MAXVAL(SR)+0.5*dSR(ISR(1),ISR(2))+1.D3) err=1.D-12
                           IF(PR(I)<MAXVAL(SR)+0.5*dSR(ISR(1),ISR(2))+1.D3)  err=1.D-16
                           !------------------------------------------------------ 
                        DO J=1,NPF
                           !------------------------------------------------------
                           FF(:,:)=0.D0
                           !------------------------------------------------------
                           DO JJ=1,NSF 
                               !---------------------------------------------------
                               SF0=SF(JJ);       dSF0=dSF(JJ)
                               !--------------------------------------------------- 
                           DO II=1,NR(JJ)
                               !--------------------------------------------------- 
                               NP0=NP(II,JJ);    SR0=SR(II,JJ);    dSR0=dSR(II,JJ)
                               !---------------------------------------------------
                               ALLOCATE(WM(0:NP0,0:NP0))
                               ALLOCATE(MV(NPL+NSL,0:NP0))
                               ALLOCATE(MA(NPL+NSL,0:NP0,3))
                               ALLOCATE(MG(NPL+NSL,0:NP0,6))
                               !---------------------------------------------------
                               CALL Lagrange_Coeff(NP0,SR0,dSR0,WM)
                               !---------------------------------------------------  
                               CALL WeightCoef1(NP0,Norm,NPL,NSL,PR(I),PF(J),PL,SR0,SF0,SL,dSR0,dSF0,dSL,MV,MA,MG,WM)
                               !--------------------------------------------------- 
                               DO KK=0,NP0                  
                                  !------------------------------------------------ 
                                  RoC(:)=Ro(KK+1,II,JJ,:)
                                  !------------------------------------------------ 
                                  CALL  BCE3(NSL,NPL,MV(:,KK),RoC,F);    FF(1,:)=FF(1,:)+F(:)
                                  !------------------------------------------------ 
                                  CALL  BCE3(NSL,NPL,MA(:,KK,1),RoC,F);  FF(2,:)=FF(2,:)+F(:)
                                  CALL  BCE3(NSL,NPL,MA(:,KK,2),RoC,F);  FF(3,:)=FF(3,:)+F(:)
                                  CALL  BCE3(NSL,NPL,MA(:,KK,3),RoC,F);  FF(4,:)=FF(4,:)+F(:)
                                  !------------------------------------------------ 
                                  CALL  BCE3(NSL,NPL,MG(:,KK,1),RoC,F);  FF(5,:)=FF(5,:)+F(:)
                                  CALL  BCE3(NSL,NPL,MG(:,KK,2),RoC,F);  FF(6,:)=FF(6,:)+F(:)
                                  CALL  BCE3(NSL,NPL,MG(:,KK,3),RoC,F);  FF(7,:)=FF(7,:)+F(:)
                                  CALL  BCE3(NSL,NPL,MG(:,KK,4),RoC,F);  FF(8,:)=FF(8,:)+F(:)
                                  CALL  BCE3(NSL,NPL,MG(:,KK,5),RoC,F);  FF(9,:)=FF(9,:)+F(:)
                                  CALL  BCE3(NSL,NPL,MG(:,KK,6),RoC,F);  FF(10,:)=FF(10,:)+F(:)
                                  !------------------------------------------------ 
                               END DO
                               !--------------------------------------------------- 
                               DEALLOCATE(WM,MV,MA,MG)
                               !--------------------------------------------------- 
                           END DO
                           !------------------------------------------------------- 
                           END DO
                           !------------------------------------------------------- 
                           GV(I,J,:)=FF(1,:)
                           GA(:,I,J,:)=FF(2:4,:)
                           GG(:,I,J,:)=FF(5:10,:)
                           !------------------------------------------------------- 
                        END DO
                        END DO
                        !---------------------------------------------------------   
                        
                    END IF
                    !--------------------------------------------------------------
                
                    
                 
                    CALL SYSTEM_CLOCK(T1)
  
                    Comp_time=(T1-T0)/10000.
                    
                    WRITE(*,*)
                    WRITE(*,*) 'Time for modelling��s����', Comp_time
                    WRITE(*,*)


                
!-------------------------------------------------------------------------------------------
                                                                                    RETURN
END SUBROUTINE
!===========================================================================================