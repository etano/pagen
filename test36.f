CCCTEST OF COULCC36-f90
      use CSTEED
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 X,ETA,ZLMIN,FC(501),GC(501),FCP(501),GCP(501),
     X           SIG(501),ZL,WS,CI
      LOGICAL WHIT
!      INTEGER NFP,N11,NPQ(2),N20,KAS(2)
      CHARACTER*20 NOTE
      CHARACTER*4 WHO(2,4,3),IRREG,REG
      DATA ZERO,HALF,ONE,FOUR,CI / 0D+0,0.5D+0,1D+0,4D+0,(0D+0,1D+0) /
      DATA WHO / 'F','G','j','y','J','Y','I','K' ,
     X           'F','H+','j','h(1)','J','H(1)','?','?' ,
     X           'F','H-','j','h(2)','J','H(2)','?','?' /
C
      PI = FOUR * ATAN(ONE)
      WRITE(6,1000)
C
C   10 READ(5,*,END=40) X,ETA,ZLMIN,NL,MODE,KFN,WHIT,NOTE
   10 X=0.4472135955
      ETA=0.111803398875
      ZLMIN=0.
      NL=1
      MODE=3
      KFN=0
      WHIT=.FALSE.
      NOTE='hi'
      IF(NL.LE.0) GO TO 40
      IFAIL = 1
      MD = MOD(ABS(MODE),10)
      IH =     ABS(MODE)/10
      KFIN = KFN
      IF(WHIT) KFN=MAX(KFN,0)
      WRITE(6,1010) X,ETA,ZLMIN,NL,MODE,KFN,NOTE
C coulcc.coulcc(r2k,fo1,-.5,lmax,1,0,ifail)
      CALL COULCC(X,ETA,ZLMIN,NL,FC,GC,FCP,GCP,SIG,MODE,KFN,IFAIL)
C
      WRITE(6,1020) IFAIL,RERR,NFP,N11,NPQ,N20,KAS
      IF(IFAIL.LT.0) GO TO 30
      DO 20 I=1,NL
      L = I
      WRITE(6,*) L
      IF(L.GT.NL-IFAIL) GO TO 20
         ZL = ZLMIN + L - 1
         IF(KFN.NE.0) SIG(L) = ZERO
         IRREG = WHO(2,MAX(KFN+1,1),IH+1)
           REG = WHO(1,MAX(KFN+1,1),1)
         IF(WHIT) THEN
            IRREG = 'WHIT'
            WS = EXP(-HALF*PI*(ETA - CI*ZL)  - CI*SIG(L))
            GC(L)  = WS * GC(L)
            IF(MD.EQ.1) GCP(L) = CI*WS * GCP(L)
            FC(L)  = CI/WS * FC(L)
            IF(MOD(MD,2).EQ.1) FCP(L)  = FCP(L) / WS
             REG = 'WH-F'
           ENDIF
      WRITE(6,1030) ZL,REG,FC(L),IRREG,GC(L)
      WRITE(6,1040)            FCP(L),GCP(L)
      WRITE(6,1050) SIG(L)
   20 CONTINUE
   30 CONTINUE
   40 STOP
 1000 FORMAT('1TEST OF THE CONTINUED-FRACTION COULOMB & BESSEL ROUTINES'
     X /)
 1010 FORMAT(/'0X =',F12.4,F9.4,', ETA =',2F8.3,', ZLMIN =',2F8.3,'  NL
     X=',I4,  '  MODE = ',I3,'  KFN =',I2,8X,A20)
 1020 FORMAT(' COULCC  :: IFAIL =',I4,
     X   '  RERR =',1P,E10.1,'  ITS =',I6,I4,2I6,I4,2I2/)
 1030 FORMAT(' ZL =',2F8.3,' :: FC = ',A4,' =',1P,2D20.12,',  GC = ',
     X   A4,' =',2D20.12)
 1040 FORMAT(24X,' FC''=       ',1P,2D20.12,',  GC'' =',6X ,2D20.12)
 1050 FORMAT(24X,' SIG=   ',   2G20.12)
      END

