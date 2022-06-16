using DifferentialEquations
using ModelingToolkit
using Plots
Plots.gr(lw=2)

# fact is a normalized-Hill activation function
function act(X, W, n, EC_50)
    beta = ((EC_50^n) - 1)/(2*(EC_50^n) - 1)
    K = (beta-1)^(1/n)
    
    return W*((beta*(X^n))/((K^n)+(X^n)))
end

inhib(act_n_hill_fn, W) = W - act_n_hill_fn

AND(a, b) = a*b

OR(a, b) = a + b - a*b

# create ODE System using ModelingToolkit

# --- species ---
@parameters(
    max_AC, tau_AC,             # 1. AC
    max_B1AR, tau_B1AR,         # 2 .B1AR
    max_B1ARPA, tau_B1ARPA,     # 3. B1ARPA
    max_B1ARPG, tau_B1ARPG,     # 4. B1ARPG
    max_cAMP, tau_cAMP,         # 5. cAMP
    max_Fsk, tau_Fsk,           # 6. Fsk
    max_Gbg, tau_Gbg,           # 7. Gbg (in Excel) / Gsbg (in PDF)
    max_GRK, tau_GRK,           # 8. GRK
    max_Gsa, tau_Gsa,           # 9. Gsa
    max_GsaBg, tau_GsaBg,       # 10. GsaBg
    max_IBMX, tau_IBMX,         # 11. IBMX
    max_ICa, tau_ICa,           # 12. ICa
    max_IKs, tau_IKs,           # 13. IKs
    max_Inhib1, tau_Inhib1,     # 14. Inhib1
    max_NE, tau_NE,             # 15. NE
    max_PDE, tau_PDE,           # 16. PDE
    max_PKA, tau_PKA,           # 17. PKA
    max_PKAC, tau_PKAC,         # 18. PKAC
    max_PKAR, tau_PKAR,         # 19. PKAR
    max_PKI, tau_PKI,           # 20. PKI
    max_PLB, tau_PLB,           # 21. PLB
    max_PP1, tau_PP1,           # 22. PP1
    max_PP2A, tau_PP2A,         # 23. PP2A
    max_RyR, tau_RyR,           # 24. RyR
    max_TnI, tau_TnI )          # 25. TnI

# --- reactions ---
@parameters(
    r1_W, r1_n, r1_EC_50,      #  r1 = [ => NE ]
    r2_W, r2_n, r2_EC_50,      #  r2 = [ => Fsk ]
    r3_W, r3_n, r3_EC_50,      #  r3 = [ => IBMX ]
    r4_W, r4_n, r4_EC_50,      #  r4 = [ => PKI ]
    r5_W, r5_n, r5_EC_50,      #  r5 = [ => PP2A ]
    r6_W, r6_n, r6_EC_50,      #  r6 = [ NE & !B1ARPG & !B1ARPA => B1AR ]
    r7_W, r7_n, r7_EC_50,      #  r7 = [ B1AR => GRK ]
    r8_W, r8_n, r8_EC_50,      #  r8 = [ GRK => B1ARPG ]
    r9_W, r9_n, r9_EC_50,      #  r9 = [ B1AR => GsaBg ]
    r10_W, r10_n, r10_EC_50,   # r10 = [ GsaBg => Gsa ]
    r11_W, r11_n, r11_EC_50,   # r11 = [ GsaBg => Gbg ]
    r12_W, r12_n, r12_EC_50,   # r12 = [ Gsa => AC ]
    r13_W, r13_n, r13_EC_50,   # r13 = [ Fsk => AC ]
    r14_W, r14_n, r14_EC_50,   # r14 = [ AC => cAMP ]
    r15_W, r15_n, r15_EC_50,   # r15 = [ !IBMX => PDE ]
    r16_W, r16_n, r16_EC_50,   # r16 = [ !PDE => cAMP ]
    r17_W, r17_n, r17_EC_50,   # r17 = [ cAMP => PKA ]
    r18_W, r18_n, r18_EC_50,   # r18 = [ PKA => PKAR ]
    r19_W, r19_n, r19_EC_50,   # r19 = [ PKA => PKAC ]
    r20_W, r20_n, r20_EC_50,   # r20 = [ PKAC => IKs ]
    r21_W, r21_n, r21_EC_50,   # r21 = [ !PKI => PKAC ]
    r22_W, r22_n, r22_EC_50,   # r22 = [ PKAC => B1ARPA ]
    r23_W, r23_n, r23_EC_50,   # r23 = [ PKAC => TnI ]
    r24_W, r24_n, r24_EC_50,   # r24 = [ PKAC => PLB ]
    r25_W, r25_n, r25_EC_50,   # r25 = [ PKAC => Inhib1 ]
    r26_W, r26_n, r26_EC_50,   # r26 = [ PKAC => RyR ]
    r27_W, r27_n, r27_EC_50,   # r27 = [ PKAC => ICa ]
    r28_W, r28_n, r28_EC_50,   # r28 = [ !Inhib1 => PP1 ]
    r29_W, r29_n, r29_EC_50,   # r29 = [ !PP1 => IKs ]
    r30_W, r30_n, r30_EC_50,   # r30 = [ !PP1 => RyR ]
    r31_W, r31_n, r31_EC_50,   # r31 = [ !PP1 => PLB ]
    r32_W, r32_n, r32_EC_50,   # r32 = [ !PP1 => ICa ]
    r33_W, r33_n, r33_EC_50,   # r33 = [ !PP2A => ICa ]
    r34_W, r34_n, r34_EC_50,   # r34 = [ !PP2A => TnI ]
    r35_W, r35_n, r35_EC_50,   # r35 = [ !PP2A => Inhib1 ]
    r36_W, r36_n, r36_EC_50 )  # r36 = [ !PP2A => RyR ]

@variables(t, 
    AC(t),    B1AR(t),    B1ARPA(t),   B1ARPG(t),
    cAMP(t),  Fsk(t),     Gbg(t),      GRK(t),
    Gsa(t),   GsaBg(t),   IBMX(t),     ICa(t),
    IKs(t),   Inhib1(t),  NE(t),       PDE(t),
    PKA(t),   PKAC(t),    PKAR(t),     PKI(t),
    PLB(t),   PP1(t),     PP2A(t),     RyR(t), TnI(t) )
    

@variables(
    r6_fact_NE(t), r6_fact_B1ARPG(t), r6_fact_B1ARPA(t), r6_finhib_B1ARPG(t), r6_finhib_B1ARPA(t),
    r7_fact_B1AR(t),
    r8_fact_GRK(t),
    r9_fact_B1AR(t),
    r10_fact_GsaBg(t),
    r11_fact_GsaBg(t),
    r12_fact_Gsa(t),
    r13_fact_Fsk(t),
    r14_fact_AC(t),
    r15_fact_IBMX(t), r15_finhib_IBMX(t),
    r16_fact_PDE(t), r16_finhib_PDE(t),
    r17_fact_cAMP(t),
    r18_fact_PKA(t),
    r19_fact_PKA(t),
    r20_fact_PKAC(t),
    r21_fact_PKI(t), r21_finhib_PKI(t),
    r22_fact_PKAC(t),
    r23_fact_PKAC(t),
    r24_fact_PKAC(t),
    r25_fact_PKAC(t),
    r26_fact_PKAC(t),
    r27_fact_PKAC(t),
    r28_fact_Inhib1(t), r28_finhib_Inhib1(t),
    r29_fact_PP1(t), r29_finhib_PP1(t),
    r30_fact_PP1(t), r30_finhib_PP1(t),
    r31_fact_PP1(t), r31_finhib_PP1(t),
    r32_fact_PP1(t), r32_finhib_PP1(t),
    r33_fact_PP2A(t), r33_finhib_PP2A(t),
    r34_fact_PP2A(t), r34_finhib_PP2A(t),
    r35_fact_PP2A(t), r35_finhib_PP2A(t),
    r36_fact_PP2A(t), r36_finhib_PP2A(t),
)

Dt = Differential(t)

eqsFull = [
    # reaction 6 = [ NE & !B1ARPG & !B1ARPA => B1AR ]
    r6_fact_NE ~ act(NE, r6_W, r6_n, r6_EC_50), 
    r6_fact_B1ARPG ~ act(B1ARPG, r6_W, r6_n, r6_EC_50),
    r6_fact_B1ARPA ~ act(B1ARPA, r6_W, r6_n, r6_EC_50),
    r6_finhib_B1ARPG ~ inhib(r6_fact_B1ARPG, r6_W),
    r6_finhib_B1ARPA ~ inhib(r6_fact_B1ARPA, r6_W),
    # reaction 7 = [ B1AR => GRK ]
    r7_fact_B1AR ~ act(B1AR, r7_W, r7_n, r7_EC_50),
    # reaction 8 = [ GRK => B1ARPG ]
    r8_fact_GRK ~ act(GRK, r8_W, r8_n, r8_EC_50),
    # reaction 9 = [ B1AR => GsaBg ]
    r9_fact_B1AR ~ act(B1AR, r9_W, r9_n, r9_EC_50),
    # reaction 10 = [ GsaBg => Gsa ]
    r10_fact_GsaBg ~ act(GsaBg, r10_W, r10_n, r10_EC_50),
    # reaction 11 = [ GsaBg => Gbg ]
    r11_fact_GsaBg ~ act(GsaBg, r11_W, r11_n, r11_EC_50),
    # reaction 12 = [ Gsa => AC ]
    r12_fact_Gsa ~ act(Gsa, r12_W, r12_n, r12_EC_50),
    # reaction 13 = [ Fsk => AC ]
    r13_fact_Fsk ~ act(Fsk, r13_W, r13_n, r13_EC_50),
    # reaction 14 = [ AC => cAMP ]
    r14_fact_AC ~ act(AC, r14_W, r14_n, r14_EC_50),
    # reaction 15 = [ !IBMX => PDE ]
    r15_fact_IBMX ~ act(IBMX, r15_W, r15_n, r15_EC_50),
    r15_finhib_IBMX ~ inhib(r15_fact_IBMX, r15_W),
    # reaction 16 = [ !PDE => cAMP ]
    r16_fact_PDE ~ act(PDE, r16_W, r16_n, r16_EC_50),
    r16_finhib_PDE ~ inhib(r16_fact_PDE, r16_W),
    # reaction 17 = [ cAMP => PKA ]
    r17_fact_cAMP ~ act(cAMP, r17_W, r17_n, r17_EC_50),
    # reaction 18 = [ PKA => PKAR ]
    r18_fact_PKA ~ act(PKA, r18_W, r18_n, r18_EC_50),
    # reaction 19 = [ PKA => PKAC ]
    r19_fact_PKA ~ act(PKA, r19_W, r19_n, r19_EC_50),
    # reaction 20 = [ PKAC => IKs ]
    r20_fact_PKAC ~ act(PKAC, r20_W, r20_n, r20_EC_50),
    # reaction 21 = [ !PKI => PKAC ]
    r21_fact_PKI ~ act(PKI, r21_W, r21_n, r21_EC_50),
    r21_finhib_PKI ~ inhib(r21_fact_PKI, r21_W),
    # reaction 22 = [ PKAC => B1ARPA ]
    r22_fact_PKAC ~ act(PKAC, r22_W, r22_n, r22_EC_50),
    # reaction 23 = [ PKAC => TnI ]
    r23_fact_PKAC ~ act(PKAC, r23_W, r23_n, r23_EC_50),
    # reaction 24 = [ PKAC => PLB ]
    r24_fact_PKAC ~ act(PKAC, r24_W, r24_n, r24_EC_50),
    # reaction 25 = [ PKAC => Inhib1 ]
    r25_fact_PKAC ~ act(PKAC, r25_W, r25_n, r25_EC_50),
    # reaction 26 = [ PKAC => RyR ]
    r26_fact_PKAC ~ act(PKAC, r26_W, r26_n, r26_EC_50),
    # reaction 27 = [ PKAC => ICa ]
    r27_fact_PKAC ~ act(PKAC, r27_W, r27_n, r27_EC_50),
    # reaction 28 = [ !Inhib1 => PP1 ]
    r28_fact_Inhib1 ~ act(Inhib1, r28_W, r28_n, r28_EC_50),
    r28_finhib_Inhib1 ~ inhib(r28_fact_Inhib1, r28_W),
    # reaction 29 = [ !PP1 => IKs ]
    r29_fact_PP1 ~ act(PP1, r29_W, r29_n, r29_EC_50),
    r29_finhib_PP1 ~ inhib(r29_fact_PP1, r29_W),
    # reaction 30 = [ !PP1 => RyR ]
    r30_fact_PP1 ~ act(PP1, r30_W, r30_n, r30_EC_50),
    r30_finhib_PP1 ~ inhib(r30_fact_PP1, r30_W),
    # reaction 31 = [ !PP1 => PLB ]
    r31_fact_PP1 ~ act(PP1, r31_W, r31_n, r31_EC_50),
    r31_finhib_PP1 ~ inhib(r31_fact_PP1, r31_W),
    # reaction 32 = [ !PP1 => ICa ]
    r32_fact_PP1 ~ act(PP1, r32_W, r32_n, r32_EC_50),
    r32_finhib_PP1 ~ inhib(r32_fact_PP1, r32_W),
    # reaction 33 = [ !PP2A => ICa ]
    r33_fact_PP2A ~ act(PP2A, r33_W, r33_n, r33_EC_50),
    r33_finhib_PP2A ~ inhib(r33_fact_PP2A, r33_W),
    # reaction 34 = [ !PP2A => TnI ]
    r34_fact_PP2A ~ act(PP2A, r34_W, r34_n, r34_EC_50),
    r34_finhib_PP2A ~ inhib(r34_fact_PP2A, r34_W),
    # reaction 35 = [ !PP2A => Inhib1 ]
    r35_fact_PP2A ~ act(PP2A, r35_W, r35_n, r35_EC_50),
    r35_finhib_PP2A ~ inhib(r35_fact_PP2A, r35_W),
    # reaction 36 = [ !PP2A => RyR ]
    r36_fact_PP2A ~ act(PP2A, r36_W, r36_n, r36_EC_50),
    r36_finhib_PP2A ~ inhib(r36_fact_PP2A, r36_W),

    Dt(NE)   ~  (r1_W*max_NE - NE)      / tau_NE,
    Dt(Fsk)  ~  (r2_W*max_Fsk - Fsk)    / tau_Fsk,
    Dt(IBMX) ~  (r3_W*max_IBMX - IBMX)  / tau_IBMX,
    Dt(PKI)  ~  (r4_W*max_PKI - PKI)    / tau_PKI,
    Dt(PP2A) ~  (r5_W*max_PP2A - PP2A)  / tau_PP2A,

    Dt(B1AR) ~ (r6_W*AND(r6_fact_NE, AND(r6_finhib_B1ARPG, r6_finhib_B1ARPA))*max_B1AR - B1AR) / tau_B1AR,
    Dt(GRK) ~ (r7_W*r7_fact_B1AR*max_GRK - GRK) / tau_GRK,
    Dt(B1ARPG) ~ (r8_W*r8_fact_GRK*max_B1ARPG - B1ARPG) / tau_B1ARPG,
    Dt(GsaBg) ~ (r9_W*r9_fact_B1AR*max_GsaBg - GsaBg) / tau_GsaBg,
    Dt(Gsa) ~ (r10_W*r10_fact_GsaBg*max_Gsa - Gsa) / tau_Gsa,
    Dt(Gbg) ~ (r11_W*r11_fact_GsaBg*max_Gbg - Gbg) / tau_Gbg,
    Dt(AC) ~ (OR(r12_W*r12_fact_Gsa, r13_W*r13_fact_Fsk)*max_AC - AC) / tau_AC,
    Dt(cAMP) ~ (OR(r14_W*r14_fact_AC, r16_W*r16_finhib_PDE)*max_cAMP - cAMP) / tau_cAMP,
    Dt(PDE) ~ (r15_W*r15_finhib_IBMX*max_PDE - PDE) / tau_PDE,
    Dt(PKA) ~ (r17_W*r17_fact_cAMP*max_PKA - PKA) / tau_PKA,
    Dt(PKAR) ~ (r18_W*r18_fact_PKA*max_PKAR - PKAR) / tau_PKAR,
    Dt(PKAC) ~ (OR(r19_W*r19_fact_PKA, r21_W*r21_finhib_PKI)*max_PKAC - PKAC) / tau_PKAC,
    Dt(IKs) ~ (OR(r20_W*r20_fact_PKAC, r29_W*r29_finhib_PP1)*max_IKs - IKs) / tau_IKs,
    Dt(B1ARPA) ~ (r22_W*r22_fact_PKAC*max_B1ARPA - B1ARPA) / tau_B1ARPA,
    Dt(TnI) ~ (OR(r23_W*r23_fact_PKAC, r34_W*r34_finhib_PP2A)*max_TnI - TnI) / tau_TnI,
    Dt(PLB) ~ (OR(r24_W*r24_fact_PKAC, r31_W*r31_finhib_PP1)*max_PLB - PLB) / tau_PLB,
    Dt(Inhib1) ~ (OR(r25_W*r25_fact_PKAC, r35_W*r35_finhib_PP2A)*max_Inhib1 - Inhib1) / tau_Inhib1,
    Dt(RyR) ~ (OR(r26_W*r26_fact_PKAC, OR(r30_W*r30_finhib_PP1, r36_W*r36_finhib_PP2A))*max_RyR - RyR) / tau_RyR,
    Dt(ICa) ~ (OR(r27_W*r27_fact_PKAC, OR(r32_W*r32_finhib_PP1, r33_W*r33_finhib_PP2A))*max_ICa - ICa) / tau_ICa,
    Dt(PP1) ~ (r28_W*r28_finhib_Inhib1*max_PP1 - PP1) / tau_PP1,
]


# assign full ODESystem to fix parameter index/position
@named fullsys = ODESystem( eqsFull, t, 
    [ AC,   B1AR, B1ARPA, B1ARPG, cAMP,  
      Fsk,  Gbg,  GRK,    Gsa,    GsaBg,
      IBMX, ICa,  IKs,    Inhib1, NE,
      PDE,  PKA,  PKAC,   PKAR,   PKI,
      PLB,  PP1,  PP2A,   RyR,    TnI ],
    [ max_AC, tau_AC,             # 1. AC
      max_B1AR, tau_B1AR,         # 2 .B1AR
      max_B1ARPA, tau_B1ARPA,     # 3. B1ARPA
      max_B1ARPG, tau_B1ARPG,     # 4. B1ARPG
      max_cAMP, tau_cAMP,         # 5. cAMP
      max_Fsk, tau_Fsk,           # 6. Fsk
      max_Gbg, tau_Gbg,           # 7. Gbg (in Excel) / Gsbg (in PDF)
      max_GRK, tau_GRK,           # 8. GRK
      max_Gsa, tau_Gsa,           # 9. Gsa
      max_GsaBg, tau_GsaBg,       # 10. GsaBg
      max_IBMX, tau_IBMX,         # 11. IBMX
      max_ICa, tau_ICa,           # 12. ICa
      max_IKs, tau_IKs,           # 13. IKs
      max_Inhib1, tau_Inhib1,     # 14. Inhib1
      max_NE, tau_NE,             # 15. NE
      max_PDE, tau_PDE,           # 16. PDE
      max_PKA, tau_PKA,           # 17. PKA
      max_PKAC, tau_PKAC,         # 18. PKAC
      max_PKAR, tau_PKAR,         # 19. PKAR
      max_PKI, tau_PKI,           # 20. PKI
      max_PLB, tau_PLB,           # 21. PLB
      max_PP1, tau_PP1,           # 22. PP1
      max_PP2A, tau_PP2A,         # 23. PP2A
      max_RyR, tau_RyR,           # 24. RyR
      max_TnI, tau_TnI,           # 25. TnI
      r1_W, r1_n, r1_EC_50,      #  r1 = [ => NE ]
      r2_W, r2_n, r2_EC_50,      #  r2 = [ => Fsk ]
      r3_W, r3_n, r3_EC_50,      #  r3 = [ => IBMX ]
      r4_W, r4_n, r4_EC_50,      #  r4 = [ => PKI ]
      r5_W, r5_n, r5_EC_50,      #  r5 = [ => PP2A ]
      r6_W, r6_n, r6_EC_50,      #  r6 = [ NE & !B1ARPG & !B1ARPA => B1AR ]
      r7_W, r7_n, r7_EC_50,      #  r7 = [ B1AR => GRK ]
      r8_W, r8_n, r8_EC_50,      #  r8 = [ GRK => B1ARPG ]
      r9_W, r9_n, r9_EC_50,      #  r9 = [ B1AR => GsaBg ]
      r10_W, r10_n, r10_EC_50,   # r10 = [ GsaBg => Gsa ]
      r11_W, r11_n, r11_EC_50,   # r11 = [ GsaBg => Gbg ]
      r12_W, r12_n, r12_EC_50,   # r12 = [ Gsa => AC ]
      r13_W, r13_n, r13_EC_50,   # r13 = [ Fsk => AC ]
      r14_W, r14_n, r14_EC_50,   # r14 = [ AC => cAMP ]
      r15_W, r15_n, r15_EC_50,   # r15 = [ !IBMX => PDE ]
      r16_W, r16_n, r16_EC_50,   # r16 = [ !PDE => cAMP ]
      r17_W, r17_n, r17_EC_50,   # r17 = [ cAMP => PKA ]
      r18_W, r18_n, r18_EC_50,   # r18 = [ PKA => PKAR ]
      r19_W, r19_n, r19_EC_50,   # r19 = [ PKA => PKAC ]
      r20_W, r20_n, r20_EC_50,   # r20 = [ PKAC => IKs ]
      r21_W, r21_n, r21_EC_50,   # r21 = [ !PKI => PKAC ]
      r22_W, r22_n, r22_EC_50,   # r22 = [ PKAC => B1ARPA ]
      r23_W, r23_n, r23_EC_50,   # r23 = [ PKAC => TnI ]
      r24_W, r24_n, r24_EC_50,   # r24 = [ PKAC => PLB ]
      r25_W, r25_n, r25_EC_50,   # r25 = [ PKAC => Inhib1 ]
      r26_W, r26_n, r26_EC_50,   # r26 = [ PKAC => RyR ]
      r27_W, r27_n, r27_EC_50,   # r27 = [ PKAC => ICa ]
      r28_W, r28_n, r28_EC_50,   # r28 = [ !Inhib1 => PP1 ]
      r29_W, r29_n, r29_EC_50,   # r29 = [ !PP1 => IKs ]
      r30_W, r30_n, r30_EC_50,   # r30 = [ !PP1 => RyR ]
      r31_W, r31_n, r31_EC_50,   # r31 = [ !PP1 => PLB ]
      r32_W, r32_n, r32_EC_50,   # r32 = [ !PP1 => ICa ]
      r33_W, r33_n, r33_EC_50,   # r33 = [ !PP2A => ICa ]
      r34_W, r34_n, r34_EC_50,   # r34 = [ !PP2A => TnI ]
      r35_W, r35_n, r35_EC_50,   # r35 = [ !PP2A => Inhib1 ]
      r36_W, r36_n, r36_EC_50 ]  # r36 = [ !PP2A => RyR ]
)

fullSys = structural_simplify(fullsys)

#= 
# NOTE: mentioned in article, help for tuning
#
# From these three constraints, 
#
#   fact(0) = 0, fact(EC50) = 0.5 and fact(1) = 1.
#
# derived --> beta = (EC_50^n - 1)/(2*EC_50^n -1), K = (beta-1)^(1/n)
#
# further constrained fact(X) = 1 for X ≥ 1 to ensure that species activities are limited to YMAX.
#     As default parameters, we used W = 1, EC50 = 0.5, n = 1.4, τ = 1, and YMAX = 1.
#
# τ is the time constant for a given species
# W is the reaction weight (constrained to 0 ≤ W ≤ 1)
# YMAX is the maximal fractional activation allowing simulations of knock-down (YMAX < 1) or overexpression (YMAX > 1). (protein expression)
# Typical default reaction and species parameter values are W = 1, EC50 = 0.5, n = 1.4, τ = 1, and YMAX = 1
=#

# initial condition
u0 = [ AC => 0.0,   B1AR => 0.0, B1ARPA => 0.0, B1ARPG => 0.0, cAMP => 0.0,  
       Fsk => 0.0,  Gbg => 0.0,  GRK => 0.0,    Gsa => 0.0,    GsaBg => 0.0,
       IBMX => 0.0, ICa => 0.0,  IKs => 0.0,    Inhib1 => 0.0, NE => 0.0,
       PDE => 0.0,  PKA => 0.0,  PKAC => 0.0,   PKAR => 0.0,   PKI => 0.0,
       PLB => 0.0,  PP1 => 0.0,  PP2A => 0.0,   RyR => 0.0,    TnI => 0.0 ]

# parameters
# 
# NOTE: if weight is to small, may cause solver difficult to solve answer
# by testing, choose weight = 0.99(ON_state), 0.02(OFF_state)

params = [
    # --- species ---
    max_AC => 1.0, tau_AC => 1.0,             # 1. AC
    max_B1AR => 1.0, tau_B1AR => 1.0,         # 2 .B1AR
    max_B1ARPA => 1.0, tau_B1ARPA => 1.0,     # 3. B1ARPA
    max_B1ARPG => 1.0, tau_B1ARPG => 1.0,     # 4. B1ARPG
    max_cAMP => 1.0, tau_cAMP => 1.0,         # 5. cAMP
    max_Fsk => 1.0, tau_Fsk => 1.0,           # 6. Fsk
    max_Gbg => 1.0, tau_Gbg => 1.0,           # 7. Gbg (in Excel) / Gsbg (in PDF)
    max_GRK => 1.0, tau_GRK => 1.0,           # 8. GRK
    max_Gsa => 1.0, tau_Gsa => 1.0,           # 9. Gsa
    max_GsaBg => 1.0, tau_GsaBg => 1.0,       # 10. GsaBg
    max_IBMX => 1.0, tau_IBMX => 1.0,         # 11. IBMX
    max_ICa => 1.0, tau_ICa => 1.0,           # 12. ICa
    max_IKs => 1.0, tau_IKs => 1.0,           # 13. IKs
    max_Inhib1 => 1.0, tau_Inhib1 => 1.0,     # 14. Inhib1
    max_NE => 1.0, tau_NE => 1.0,             # 15. NE
    max_PDE => 1.0, tau_PDE => 1.0,           # 16. PDE
    max_PKA => 1.0, tau_PKA => 1.0,           # 17. PKA
    max_PKAC => 1.0, tau_PKAC => 1.0,         # 18. PKAC
    max_PKAR => 1.0, tau_PKAR => 1.0,         # 19. PKAR
    max_PKI => 1.0, tau_PKI => 1.0,           # 20. PKI
    max_PLB => 1.0, tau_PLB => 1.0,           # 21. PLB
    max_PP1 => 1.0, tau_PP1 => 1.0,           # 22. PP1
    max_PP2A => 1.0, tau_PP2A => 1.0,         # 23. PP2A
    max_RyR => 1.0, tau_RyR => 1.0,           # 24. RyR
    max_TnI => 1.0, tau_TnI => 1.0,           # 25. TnI
    
    # --- reactions ---
    # %% input %%
    r1_W => 0.02, r1_n => 1.4, r1_EC_50 => 0.5,      #  r1 = [ => NE ]
    r2_W => 0.02, r2_n => 1.4, r2_EC_50 => 0.5,      #  r2 = [ => Fsk ]
    r3_W => 0.02, r3_n => 1.4, r3_EC_50 => 0.5,      #  r3 = [ => IBMX ]
    r4_W => 0.02, r4_n => 1.4, r4_EC_50 => 0.5,      #  r4 = [ => PKI ]
    r5_W => 0.02, r5_n => 1.4, r5_EC_50 => 0.5,      #  r5 = [ => PP2A ]
    # %% middle %%
    r6_W => 0.98, r6_n => 1.4, r6_EC_50 => 0.5,      #  r6 = [ NE & !B1ARPG & !B1ARPA => B1AR ]
    r7_W => 0.98, r7_n => 1.4, r7_EC_50 => 0.5,      #  r7 = [ B1AR => GRK ]
    r8_W => 0.98, r8_n => 1.4, r8_EC_50 => 0.5,      #  r8 = [ GRK => B1ARPG ]
    r9_W => 0.98, r9_n => 1.4, r9_EC_50 => 0.5,      #  r9 = [ B1AR => GsaBg ]
    r10_W => 0.98, r10_n => 1.4, r10_EC_50 => 0.5,   # r10 = [ GsaBg => Gsa ]
    r11_W => 0.98, r11_n => 1.4, r11_EC_50 => 0.5,   # r11 = [ GsaBg => Gbg ]
    r12_W => 0.98, r12_n => 1.4, r12_EC_50 => 0.5,   # r12 = [ Gsa => AC ]
    r13_W => 0.98, r13_n => 1.4, r13_EC_50 => 0.5,   # r13 = [ Fsk => AC ]
    r14_W => 0.98, r14_n => 1.4, r14_EC_50 => 0.5,   # r14 = [ AC => cAMP ]
    r15_W => 0.98, r15_n => 1.4, r15_EC_50 => 0.5,   # r15 = [ !IBMX => PDE ]
    r16_W => 0.98, r16_n => 1.4, r16_EC_50 => 0.5,   # r16 = [ !PDE => cAMP ]
    r17_W => 0.98, r17_n => 1.4, r17_EC_50 => 0.5,   # r17 = [ cAMP => PKA ]
    r18_W => 0.98, r18_n => 1.4, r18_EC_50 => 0.5,   # r18 = [ PKA => PKAR ]
    r19_W => 0.98, r19_n => 1.4, r19_EC_50 => 0.5,   # r19 = [ PKA => PKAC ]
    r20_W => 0.98, r20_n => 1.4, r20_EC_50 => 0.5,   # r20 = [ PKAC => IKs ]
    r21_W => 0.98, r21_n => 1.4, r21_EC_50 => 0.5,   # r21 = [ !PKI => PKAC ]
    r22_W => 0.98, r22_n => 1.4, r22_EC_50 => 0.5,   # r22 = [ PKAC => B1ARPA ]
    r23_W => 0.98, r23_n => 1.4, r23_EC_50 => 0.5,   # r23 = [ PKAC => TnI ]
    r24_W => 0.98, r24_n => 1.4, r24_EC_50 => 0.5,   # r24 = [ PKAC => PLB ]
    r25_W => 0.98, r25_n => 1.4, r25_EC_50 => 0.5,   # r25 = [ PKAC => Inhib1 ]
    r26_W => 0.98, r26_n => 1.4, r26_EC_50 => 0.5,   # r26 = [ PKAC => RyR ]
    r27_W => 0.98, r27_n => 1.4, r27_EC_50 => 0.5,   # r27 = [ PKAC => ICa ]
    r28_W => 0.98, r28_n => 1.4, r28_EC_50 => 0.5,   # r28 = [ !Inhib1 => PP1 ]
    r29_W => 0.98, r29_n => 1.4, r29_EC_50 => 0.5,   # r29 = [ !PP1 => IKs ]
    r30_W => 0.98, r30_n => 1.4, r30_EC_50 => 0.5,   # r30 = [ !PP1 => RyR ]
    r31_W => 0.98, r31_n => 1.4, r31_EC_50 => 0.5,   # r31 = [ !PP1 => PLB ]
    r32_W => 0.98, r32_n => 1.4, r32_EC_50 => 0.5,   # r32 = [ !PP1 => ICa ]
    r33_W => 0.98, r33_n => 1.4, r33_EC_50 => 0.5,   # r33 = [ !PP2A => ICa ]
    r34_W => 0.98, r34_n => 1.4, r34_EC_50 => 0.5,   # r34 = [ !PP2A => TnI ]
    r35_W => 0.98, r35_n => 1.4, r35_EC_50 => 0.5,   # r35 = [ !PP2A => Inhib1 ]
    r36_W => 0.98, r36_n => 1.4, r36_EC_50 => 0.5 ]  # r36 = [ !PP2A => RyR ]

# time range
tend = 60.0
tspan = (0.0, tend)


# initial ODE Problem
prob = ODEProblem(fullSys, u0, tspan, params)

#= 
INDEX_TABLE for sol1

    %%% input %%%
    1. NE
    2. Fsk
    3. IBMX
    4. PKI
    5. PP2A

    %%% middle %%%
    6. B1AR
    7. GRK
    8. B1ARPG
    9. GsaBg
    10. Gsa
    11. Gbg
    12. AC
    13. cAMP
    14. PDE
    15. PKA
    16. PKAR
    17. PKAC
    18. IKs
    19. B1ARPA
    20. TnI
    21. PLB
    22. Inhib1
    23. RyR
    24. ICa
    25. PP1

=#

# solve problem
@time sol1 = solve(prob, QNDF())

# plot all together
plot(sol1, vars=([1, 10, 13, 21]),
     xticks = 0:10:tend, yticks = 0:0.1:1, x_tickfontsize = 6,
     xlims=(0, tend), ylims=(0.0, 1.1),
     title="cardiac β-adrenergic signaling (static)", legend=:outerbottomright, 
     xlabel="Time (sec)", ylabel="Fractioanl activation"
)