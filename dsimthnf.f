      SUBROUTINE DSIMTHNF(N, WTH, HA, S0, PR, PSI1, PSI2, PSI3, PSI4,
     &                    P0, HG, EC, BETA, ALPHAE, RD, LAMBDA, MI,
     &                    RHO1, RHO2, RHO3, RHO4, RHOF,
     &                    K1, K2, K3, K4, KF,
     &                    CP1, CP2, CP3, CP4, CPF,
     &                    SIG1, SIG2, SIG3, SIG4, SIGF,
     &                    TOL, MAXIT, W,
     &                    F, P, HT, G, S, Q1, Q2)

      INTEGER N, I, MAXIT, ITER, NN
      DOUBLE PRECISION WTH, HA, S0, PR, PSI1, PSI2, PSI3, PSI4
      DOUBLE PRECISION P0, HG, EC, BETA, ALPHAE, RD, LAMBDA, MI
      DOUBLE PRECISION RHO1, RHO2, RHO3, RHO4, RHOF
      DOUBLE PRECISION K1, K2, K3, K4, KF
      DOUBLE PRECISION CP1, CP2, CP3, CP4, CPF
      DOUBLE PRECISION SIG1, SIG2, SIG3, SIG4, SIGF
      DOUBLE PRECISION TOL, Q1, Q2
      DOUBLE PRECISION X(1000), F(1000), P(1000), HT(1000)
      DOUBLE PRECISION G(1000), S(1000)
      DOUBLE PRECISION FOLD(1000), POLD(1000), HTOLD(1000)
      DOUBLE PRECISION SOLD(1000), GOLD(1000)
      DOUBLE PRECISION T1, T2, T3, T4, T5
      DOUBLE PRECISION H, D0, D00, D000, D0000
      DOUBLE PRECISION D01, D02, D03, D04
      DOUBLE PRECISION KNF, KK, DKHNF, DKHNF4
      DOUBLE PRECISION SNF, SGMA, DSGMA, DSGMA4
      DOUBLE PRECISION W, CHECK1, CHECK2, CHECK3, CHECK4, CHECK5, CHEQ
      DOUBLE PRECISION A, B, C, D, A0, A1, A2, DD, D1, D2, D3
      DOUBLE PRECISION TEMP1, TEMP2

C --- Geometry and Initialization
      H = WTH / (N - 1)
      NN = 3

      DO I = 1, N
         X(I) = (I - 1) * H
      END DO

C --- Thermophysical property ratios
      D0 = RHO1 / RHOF
      D00 = RHO2 / RHOF
      D000 = RHO3 / RHOF
      D0000 = RHO4 / RHOF

      D01 = (RHO1 * CP1) / (RHOF * CPF)
      D02 = (RHO2 * CP2) / (RHOF * CPF)
      D03 = (RHO3 * CP3) / (RHOF * CPF)
      D04 = (RHO4 * CP4) / (RHOF * CPF)

      T1 = ((1.0D0 - PSI1)**2.5D0) * ((1.0D0 - PSI2)**2.5D0) *
     &     ((1.0D0 - PSI3)**2.5D0) / ((1.0D0 - PSI4)**2.5D0)

      T2 = (1.0D0 - PSI4) * ((1.0D0 - PSI3) *
     &     ((1.0D0 - PSI2) * (1.0D0 - PSI1 + PSI1*D0) + PSI2*D00)
     &     + PSI3*D000) + PSI4*D0000

      T3 = (1.0D0 - PSI4) * ((1.0D0 - PSI3) *
     &     ((1.0D0 - PSI2) * (1.0D0 - PSI1 + PSI1*D01) + PSI2*D02)
     &     + PSI3*D03) + PSI4*D04

C --- Effective thermal conductivity
      KNF = (K1 + (NN - 1)*KF - (NN - 1)*PSI1*(KF - K1)) /
     &      (K1 + (NN - 1)*KF + PSI1*(KF - K1))

      KK = KNF * (K2 + (NN - 1)*KNF - (NN - 1)*PSI2*(KNF - K2)) /
     &     (K2 + (NN - 1)*KNF + PSI2*(KNF - K2))

      DKHNF = KK * (K3 + (NN - 1)*KK - (NN - 1)*PSI3*(KK - K3)) /
     &       (K3 + (NN - 1)*KK + PSI3*(KK - K3))

      TEMP1 = K4 + (NN - 1)*DKHNF - (NN - 1)*PSI4*(DKHNF - K4)
      TEMP2 = K4 + (NN - 1)*DKHNF + PSI4*(DKHNF - K4)
      DKHNF4 = DKHNF * TEMP1 / TEMP2

      T4 = DKHNF4

C --- Effective electrical conductivity
      SNF = (SIG1 + 2*SIGF - 2*PSI1*(SIGF - SIG1)) /
     &      (SIG1 + 2*SIGF + PSI1*(SIGF - SIG1))

      SGMA = SNF * (SIG2 + 2*SNF - 2*PSI2*(SNF - SIG2)) /
     &       (SIG2 + 2*SNF + PSI2*(SNF - SIG2))

      DSGMA = SGMA * (SIG3 + 2*SGMA - 2*PSI3*(SGMA - SIG3)) /
     &        (SIG3 + 2*SGMA + PSI3*(SGMA - SIG3))

      DSGMA4 = DSGMA * (SIG4 + 2*DSGMA - 2*PSI4*(DSGMA - SIG4)) /
     &         (SIG4 + 2*DSGMA + PSI4*(DSGMA - SIG4))

      T5 = DSGMA4

C --- Initial conditions
      F(1) = S0
      F(N) = 0.0D0
      P(1) = 1.0D0
      P(N) = 0.0D0
      HT(1) = 1.0D0
      HT(N) = 0.0D0
      G(1) = 0.0D0
      G(N) = 0.0D0
      S(1) = 0.0D0
      S(N) = 1.0D0

      DO I = 2, N - 1
         F(I) = 0.0D0
         P(I) = 0.0D0
         HT(I) = 0.0D0
         G(I) = 0.0D0
         S(I) = 0.0D0
      END DO

      DO I = 1, N
         FOLD(I) = F(I)
         POLD(I) = P(I)
         HTOLD(I) = HT(I)
         GOLD(I) = G(I)
         SOLD(I) = S(I)
      END DO

C --- SOR iteration loop
      DO ITER = 1, MAXIT

C --- Momentum equation
         DO I = 2, N - 1
            D1 = 1.0D0 + 1.0D0 / BETA
            D2 = 1.0D0 + 2.0D0 * ALPHAE * X(I)
            A = 4*D1*D2*T1/T2 + 2*H*H*(P(I) + (T5/T2)*HA**2
     &          - (T1/T2)*P0)
            B = 2*D1*D2*T1/T2 + H*(F(I) + 2*ALPHAE*D1*T1/T2)
            C = 2*D1*D2*T1/T2 - H*(F(I) + 2*ALPHAE*D1*T1/T2)
            DD = MI*(1.0D0/T3) *
     &           (0.5D0*(G(I+1)-G(I-1))**2
     &           - 2*G(I)*(G(I+1)-2*G(I)+G(I-1))
     &           - 2*H*H)

            P(I) = W/A*(B*P(I+1) + C*P(I-1) + DD) + (1.0D0 - W)*P(I)
         END DO

C --- Energy equation
         DO I = 2, N - 1
            D2 = 1.0D0 + 2.0D0 * ALPHAE * X(I)
            D3 = 1.0D0 + 4.0D0 * RD / 3.0D0
            A = 4*D2*D3*T4/T3 + 2*H*H*PR*(P(I) + HG/T3)
            B = 2*D2*D3*T4/T3 - H*(PR*F(I) + 2*ALPHAE*T4/T3)
            C = 2*D2*D3*T4/T3 + H*(PR*F(I) + 2*ALPHAE*T4/T3)
            D = D1*D2*PR*EC*((P(I+1)-P(I-1))/(2*H))**2 +
     &          HA*HA*PR*EC*P(I)**2
            HT(I) = W/A*(C*HT(I+1) + B*HT(I-1) + 2*H*H*D)
     &               + (1.0D0 - W)*HT(I)
         END DO

C --- Velocity (Simpson)
         F(2) = F(1) + (H/24.0D0)*(9.0D0*P(1) + 19.0D0*P(2)
     &          - 5.0D0*P(3) + P(4))
         DO I = 2, N - 1
            F(I+1) = F(I-1) + (H/3.0D0)*(P(I-1) + 4.0D0*P(I) + P(I+1))
         END DO

C --- Induced magnetic field (S)
         S(1) = (1.0D0/3.0D0)*(4.0D0*S(2) - S(3))
         DO I = 2, N - 1
            A0 = 4.0D0 * LAMBDA / T3
            A1 = 2.0D0 * LAMBDA / T3 + H*F(I)
            A2 = 2.0D0 * LAMBDA / T3 - H*F(I)
            DD = -2.0D0 * G(I) * (F(I+1) - 2.0D0*F(I) + F(I-1))
            S(I) = (W/A0)*(A1*S(I+1) + A2*S(I-1) + DD)
     &             + (1.0D0 - W)*S(I)
         END DO

C --- G from S (Simpson)
         G(2) = G(1) + (H/24.0D0)*(9.0D0*S(1) + 19.0D0*S(2)
     &          - 5.0D0*S(3) + S(4))
         DO I = 2, N - 1
            G(I+1) = G(I-1) + (H/3.0D0)*(S(I-1) + 4.0D0*S(I) + S(I+1))
         END DO

C --- Convergence check
         CHECK1 = 0.0D0
         CHECK2 = 0.0D0
         CHECK3 = 0.0D0
         CHECK4 = 0.0D0
         CHECK5 = 0.0D0

         DO I = 1, N
            CHECK1 = MAX(CHECK1, ABS(F(I) - FOLD(I)))
            CHECK2 = MAX(CHECK2, ABS(P(I) - POLD(I)))
            CHECK3 = MAX(CHECK3, ABS(HT(I) - HTOLD(I)))
            CHECK4 = MAX(CHECK4, ABS(G(I) - GOLD(I)))
            CHECK5 = MAX(CHECK5, ABS(S(I) - SOLD(I)))
         END DO

         TEMP1 = MAX(MAX(CHECK1, CHECK2), CHECK3)
         CHEQ  = MAX(MAX(CHECK4, CHECK5), TEMP1)

         IF (CHEQ .LT. TOL) GOTO 100

         DO I = 1, N
            FOLD(I) = F(I)
            POLD(I) = P(I)
            HTOLD(I) = HT(I)
            GOLD(I) = G(I)
            SOLD(I) = S(I)
         END DO

      END DO

100   CONTINUE

C --- Final quantities
      Q1 = ((-P(3) + 4*P(2) - 3*P(1)) / (2*H)) *
     &     (1.0D0 + 1.0D0/BETA) *
     &     ( (1.0D0 - PSI1)**(-2.5D0) ) *
     &     ( (1.0D0 - PSI2)**(-2.5D0) ) *
     &     ( (1.0D0 - PSI3)**(-2.5D0) ) *
     &     ( (1.0D0 - PSI4)**(-2.5D0) )

      Q2 = -((-HT(3) + 4*HT(2) - 3*HT(1)) / (2*H)) *
     &      (1.0D0 + 4.0D0*RD/3.0D0) * T4

      RETURN
      END
