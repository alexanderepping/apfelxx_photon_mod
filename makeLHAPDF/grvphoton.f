* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*      G R V - P H O T O N - P A R A M E T R I Z A T I O N S      *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE :                *
*              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/31             *
*             PUBLISHED IN PHYS. REV. D46 (1992) 1973             *
*                                                                 *
*    THE OUTPUT IS ALWAYS   1./ ALPHA(EM) * X * PARTON DENSITY    *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE PARTON DISTRIBUTIONS   *
*   FOR Q ** 2 BETWEEN MU ** 2 (=  0.25 / 0.30  GEV ** 2  IN LO   *
*   / HO) AND  1.E6 GEV ** 2  AND FOR X BETWEEN  1.E-5  AND  1.   *
*                                                                 *
*              HEAVY QUARK THRESHOLDS  Q(H) = M(H) :              *
*         M(C)  =  1.5,  M(B)  =  4.5,  M(T)  =  100  GEV         *
*                                                                 *
*      CORRESPONDING LAMBDA(F) VALUES FOR F ACTIVE FLAVOURS :     *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,   LAMBDA(6)  =  0.082  GEV     *
*      HO :   LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.131,   LAMBDA(6)  =  0.053  GEV     *
*                                                                 *
*      HO DISTRIBUTIONS REFER TO THE DIS(GAMMA) SCHEME, SEE :     *
*              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/26             *
*              PUBLISHED IN PHYS. REV. D45 (1992) 3986            *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
       SUBROUTINE GRVGLO (X, Q2, UL, DL, SL, CL, BL, GL)
       IMPLICIT real (A - Z)
       MU2  = 0.25
       LAM2 = 0.232 * 0.232
       S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2))
       SS = SQRT (S)
       S2 = S * S
C...X * U = X * UBAR :
       AL =  1.717
       BE =  0.641
       AK =  0.500 - 0.176 * S
       BK = 15.00  - 5.687 * SS - 0.552 * S2
       AG =  0.235 + 0.046 * SS
       BG =  0.082 - 0.051 * S  + 0.168 * S2
       C  =   0.0  + 0.459 * S
       D  =  0.354 - 0.061 * S
       E  =  4.899 + 1.678 * S
       ES =  2.046 + 1.389 * S
       UL = FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * D = X * DBAR :
       AL =  1.549
       BE =  0.782
       AK =  0.496 + 0.026 * S
       BK =  0.685 - 0.580 * SS + 0.608 * S2
       AG =  0.233 + 0.302 * S
       BG =   0.0  - 0.818 * S  + 0.198 * S2
       C  =  0.114 + 0.154 * S
       D  =  0.405 - 0.195 * S  + 0.046 * S2
       E  =  4.807 + 1.226 * S
       ES =  2.166 + 0.664 * S
       DL  =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * G :
       AL =  0.676
       BE =  1.089
       AK =  0.462 - 0.524 * SS
       BK =  5.451              - 0.804 * S2
       AG =  0.535 - 0.504 * SS + 0.288 * S2
       BG =  0.364 - 0.520 * S
       C  = -0.323              + 0.115 * S2
       D  =  0.233 + 0.790 * S  - 0.139 * S2
       E  =  0.893 + 1.968 * S
       ES =  3.432 + 0.392 * S
       GL =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * S = X * SBAR :
       SF =   0.0
       AL =  1.609
       BE =  0.962
       AK =  0.470              - 0.099 * S2
       BK =  3.246
       AG =  0.121 - 0.068 * SS
       BG = -0.090 + 0.074 * S
       C  =  0.062 + 0.034 * S
       D  =   0.0  + 0.226 * S  - 0.060 * S2
       E  =  4.288 + 1.707 * S
       ES =  2.122 + 0.656 * S
       SL =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * C = X * CBAR :
       SF =  0.888
       AL =  0.970
       BE =  0.545
       AK =  1.254 - 0.251 * S
       BK =  3.932              - 0.327 * S2
       AG =  0.658 + 0.202 * S
       BG = -0.699
       C  =  0.965
       D  =   0.0  + 0.141 * S  - 0.027 * S2
       E  =  4.911 + 0.969 * S
       ES =  2.796 + 0.952 * S
       CL =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * B = X * BBAR :
       SF =  1.351
       AL =  1.016
       BE =  0.338
       AK =  1.961 - 0.370 * S
       BK =  0.923 + 0.119 * S
       AG =  0.815 + 0.207 * S
       BG = -2.275
       C  =  1.480
       D  = -0.223 + 0.173 * S
       E  =  5.426 + 0.623 * S
       ES =  3.819 + 0.901 * S
       BL =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       RETURN
       END
C
C
       SUBROUTINE GRVGHO (X, Q2, UH, DH, SH, CH, BH, GH)
       IMPLICIT real (A - Z)
       MU2  = 0.3
       LAM2 = 0.248 * 0.248
       S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2))
       SS = SQRT (S)
       S2 = S * S
C...X * U = X * UBAR :
       AL =  0.583
       BE =  0.688
       AK =  0.449 - 0.025 * S  - 0.071 * S2
       BK =  5.060 - 1.116 * SS
       AG =  0.103
       BG =  0.319 + 0.422 * S
       C  =  1.508 + 4.792 * S  - 1.963 * S2
       D  =  1.075 + 0.222 * SS - 0.193 * S2
       E  =  4.147 + 1.131 * S
       ES =  1.661 + 0.874 * S
       UH =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * D = X * DBAR :
       AL =  0.591
       BE =  0.698
       AK =  0.442 - 0.132 * S  - 0.058 * S2
       BK =  5.437 - 1.916 * SS
       AG =  0.099
       BG =  0.311 - 0.059 * S
       C  =  0.800 + 0.078 * S  - 0.100 * S2
       D  =  0.862 + 0.294 * SS - 0.184 * S2
       E  =  4.202 + 1.352 * S
       ES =  1.841 + 0.990 * S
       DH  =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * G :
       AL =  1.161
       BE =  1.591
       AK =  0.530 - 0.742 * SS + 0.025 * S2
       BK =  5.662
       AG =  0.533 - 0.281 * SS + 0.218 * S2
       BG =  0.025 - 0.518 * S  + 0.156 * S2
       C  = -0.282              + 0.209 * S2
       D  =  0.107 + 1.058 * S  - 0.218 * S2
       E  =   0.0  + 2.704 * S
       ES =  3.071 - 0.378 * S
       GH =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * S = X * SBAR :
       SF =   0.0
       AL =  0.635
       BE =  0.456
       AK =  1.770 - 0.735 * SS - 0.079 * S2
       BK =  3.832
       AG =  0.084 - 0.023 * S
       BG =  0.136
       C  =  2.119 - 0.942 * S  + 0.063 * S2
       D  =  1.271 + 0.076 * S  - 0.190 * S2
       E  =  4.604 + 0.737 * S
       ES =  1.641 + 0.976 * S
       SH =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * C = X * CBAR :
       SF =  0.820
       AL =  0.926
       BE =  0.152
       AK =  1.142 - 0.175 * S
       BK =  3.276
       AG =  0.504 + 0.317 * S
       BG = -0.433
       C  =  3.334
       D  =  0.398 + 0.326 * S  - 0.107 * S2
       E  =  5.493 + 0.408 * S
       ES =  2.426 + 1.277 * S
       CH =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * B = X * BBAR :
       SF =  1.297
       AL =  0.969
       BE =  0.266
       AK =  1.953 - 0.391 * S
       BK =  1.657 - 0.161 * S
       AG =  1.076 + 0.034 * S
       BG = -2.015
       C  =  1.662
       D  =  0.353 + 0.016 * S
       E  =  5.713 + 0.249 * S
       ES =  3.456 + 0.673 * S
       BH =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       RETURN
       END
C
C
       SUBROUTINE GRVGH0 (X, Q2, U0, D0, S0, C0, B0, G0)
       IMPLICIT real (A - Z)
       MU2  = 0.3
       LAM2 = 0.248 * 0.248
       S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2))
       SS = SQRT (S)
       S2 = S * S
C...X * U = X * UBAR :
       AL =  1.447
       BE =  0.848
       AK =  0.527 + 0.200 * S  - 0.107 * S2
       BK =  7.106 - 0.310 * SS - 0.786 * S2
       AG =  0.197 + 0.533 * S
       BG =  0.062 - 0.398 * S  + 0.109 * S2
       C  =          0.755 * S  - 0.112 * S2
       D  =  0.318 - 0.059 * S
       E  =  4.225 + 1.708 * S
       ES =  1.752 + 0.866 * S
       U0 =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * D = X * DBAR :
       AL =  1.424
       BE =  0.770
       AK =  0.500 + 0.067 * SS - 0.055 * S2
       BK =  0.376 - 0.453 * SS + 0.405 * S2
       AG =  0.156 + 0.184 * S
       BG =   0.0  - 0.528 * S  + 0.146 * S2
       C  =  0.121 + 0.092 * S
       D  =  0.379 - 0.301 * S  + 0.081 * S2
       E  =  4.346 + 1.638 * S
       ES =  1.645 + 1.016 * S
       D0  =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * G :
       AL =  0.661
       BE =  0.793
       AK =  0.537 - 0.600 * SS
       BK =  6.389              - 0.953 * S2
       AG =  0.558 - 0.383 * SS + 0.261 * S2
       BG =   0.0  - 0.305 * S
       C  = -0.222              + 0.078 * S2
       D  =  0.153 + 0.978 * S  - 0.209 * S2
       E  =  1.429 + 1.772 * S
       ES =  3.331 + 0.806 * S
       G0 =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * S = X * SBAR :
       SF =   0.0
       AL =  1.578
       BE =  0.863
       AK =  0.622 + 0.332 * S  - 0.300 * S2
       BK =  2.469
       AG =  0.211 - 0.064 * SS - 0.018 * S2
       BG = -0.215 + 0.122 * S
       C  =  0.153
       D  =   0.0  + 0.253 * S  - 0.081 * S2
       E  =  3.990 + 2.014 * S
       ES =  1.720 + 0.986 * S
       S0 =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * C = X * CBAR :
       SF =  0.820
       AL =  0.929
       BE =  0.381
       AK =  1.228 - 0.231 * S
       BK =  3.806             - 0.337 * S2
       AG =  0.932 + 0.150 * S
       BG = -0.906
       C  =  1.133
       D  =   0.0  + 0.138 * S  - 0.028 * S2
       E  =  5.588 + 0.628 * S
       ES =  2.665 + 1.054 * S
       C0 =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C...X * B = X * BBAR :
       SF =  1.297
       AL =  0.970
       BE =  0.207
       AK =  1.719 - 0.292 * S
       BK =  0.928 + 0.096 * S
       AG =  0.845 + 0.178 * S
       BG = -2.310
       C  =  1.558
       D  = -0.191 + 0.151 * S
       E  =  6.089 + 0.282 * S
       ES =  3.379 + 1.062 * S
       B0 =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       RETURN
       END
C
C
C       FUNCTION F(X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
C       IMPLICIT REAL (A - Z)
C       SX = SQRT (X)
C       LX = ALOG (1./X)
C       F  = (X**AK * (AG + BG * SX + C * X**BK)  +  S**AL
C     1       * EXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D
C       RETURN
C       END
C
C      renamed function F to FONE
       FUNCTION FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       IMPLICIT real (A - Z)
       SX = SQRT (X)
       LX = ALOG (1./X)
       FONE  = (X**AK * (AG + BG * SX + C * X**BK)  +  S**AL
     1    * EXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D
       RETURN
       END
C
       FUNCTION FS(X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES)
       IMPLICIT real (A - Z)
       IF (S .LE. SF) THEN
          FS = 0.0
       ELSE
          SX = SQRT (X)
          LX = ALOG (1./X)
          DS = S - SF
          FS = (DS * X**AK * (AG + BG * SX + C * X**BK) + DS**AL
     1         * EXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D
       END IF
       RETURN
       END
