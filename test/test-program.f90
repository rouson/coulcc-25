PROGRAM CCTEST
!-----------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: spi=>int32, dpf=>real64 &
                              &,stdin=>input_unit,stdout=>output_unit
!-----------------------------------------------------------------------
  USE COULCC_M                                               
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
      LOGICAL :: WHIT                                                      
      INTEGER(spi) :: I,L,NL,MODE,IFAIL,MD,IH,KFN,KFIN
      COMPLEX(dpf) :: X,ETA,ZLMIN,ZL,WS
      COMPLEX(dpf),DIMENSION(201) :: FC,GC,FCP,GCP,SIG               
!-----------------------------------------------------------------------
      REAL(dpf) :: RERR
      INTEGER(spi),DIMENSION(2) :: NFP,KASE
      INTEGER(spi),DIMENSION(3) :: NPQ
      COMMON /STEED/ RERR,NFP,NPQ,KASE                                  
!-----------------------------------------------------------------------
      CHARACTER(len=20) :: NOTE                                         
      CHARACTER(len=4)  :: IRREG,REG                                  
!-----------------------------------------------------------------------
      INTEGER(spi),PARAMETER,DIMENSION(8) :: LDISP=[1,2,4,11,31,101,301,1001]
      REAL(dpf),   PARAMETER :: ZERO=0._dpf
      REAL(dpf),   PARAMETER :: HALF=0.5_dpf
      REAL(dpf),   PARAMETER :: ONE=1._dpf
      REAL(dpf),   PARAMETER :: FOUR=4._dpf
      COMPLEX(dpf),PARAMETER :: CI=CMPLX(0._dpf,1._dpf,KIND=dpf)
      REAL(dpf),   PARAMETER :: PI=FOUR*ATAN(ONE)
!-----------------------------------------------------------------------
      CHARACTER(len=4),DIMENSION(2,4,3) :: WHO                          
      DATA WHO / 'F','G','j','y','J','Y','I','K' ,         &           
       &         'F','H+','j','h(1)','J','H(1)','?','?' ,  &           
       &         'F','H-','j','h(2)','J','H(2)','?','?' /               
!-----------------------------------------------------------------------
 1000 FORMAT('1TEST OF THE CONTINUED-FRACTION COULOMB & BESSEL ROUTINES'/)
 1010 FORMAT(/'0X =',F12.4,F9.4,', ETA =',2F8.3,', ZLMIN =',2F8.3,'  NL=',I4,  '  MODE = ',I3,'  KFN =',I2,8X,A20)
 1020 FORMAT(' COULCC  :: IFAIL =',I4,'  RERR =',1P,E10.1,'  ITS =',I6,I4,2I6,I4,2I2/)
 1030 FORMAT(' ZL =',2F8.3,' :: FC = ',A4,' =',1P,2D20.12,',  GC = ',A4,' =',2D20.12)
 1040 FORMAT(24X,' FC''=       ',1P,2D20.12,',  GC'' =',6X ,2D20.12)    
 1050 FORMAT(24X,' SIG=   ',   2F20.12)                                 
!-----------------------------------------------------------------------
      WRITE(STDOUT,1000)                                                     
!-----------------------------------------------------------------------
  10  READ(STDIN,*,END=40) X,ETA,ZLMIN,NL,MODE,KFN,WHIT,NOTE                
      IF (NL.LE.0) GO TO 40                                              
      IFAIL = 1                                                         
      MD = MOD(ABS(MODE),10)                                            
      IH =     ABS(MODE)/10                                             
      KFIN = KFN                                                        
      IF (WHIT) KFN=MAX(KFN,0)                                           
      WRITE(STDOUT,1010) X,ETA,ZLMIN,NL,MODE,KFN,NOTE                        
!-----------------------------------------------------------------------
      CALL COULCC(X,ETA,ZLMIN,NL,FC,GC,FCP,GCP,SIG,MODE,KFN,IFAIL)      
!-----------------------------------------------------------------------
      WRITE(STDOUT,1020) IFAIL,RERR,NFP,NPQ,KASE                             
      IF (IFAIL.LT.0) GO TO 30                                           
      DO 20 I=1,8                                                       
        L = LDISP(I)                                                      
        IF (L.GT.NL-IFAIL) GO TO 20                                        
        ZL = ZLMIN + L - 1                                             
        IF (KFN.NE.0) SIG(L) = ZERO                                     
        IRREG = WHO(2,MAX(KFN+1,1),IH+1)                               
        REG = WHO(1,MAX(KFN+1,1),1)                                  
        IF (WHIT) THEN                                                  
          IRREG = 'WHIT'                                              
          WS = EXP(-HALF*PI*(ETA - CI*ZL)  - CI*SIG(L))               
          GC(L)  = WS * GC(L)                                         
          IF (MD.EQ.1) GCP(L) = CI*WS * GCP(L)                         
          FC(L)  = CI/WS * FC(L)                                      
          IF (MOD(MD,2).EQ.1) FCP(L)  = FCP(L) / WS                    
          REG = 'WH-F'                                               
        END IF                                                        
        WRITE(STDOUT,1030) ZL,REG,FC(L),IRREG,GC(L)                            
        IF (MD.EQ.1)WRITE(STDOUT,1040) FCP(L),GCP(L)                 
        IF (SIG(L).NE.ZERO.AND.KFIN.EQ.0) WRITE(STDOUT,1050) SIG(L)      
  20  CONTINUE                                                          
  30  CONTINUE                                                          
      GO TO 10                                                          
  40  STOP                                                              
!-----------------------------------------------------------------------
END PROGRAM CCTEST
