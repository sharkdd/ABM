#ifndef Model_H
#define Model_H

#include <string>

// +=+ Named constants +=======================================================+

// The model is stratified to exclude the population who debut before
// (ANTE) PrEP rollout ends from people who debut afterwards (POST).
// This lets us track stratified outcomes.
enum Track {ANTE=0, POST=1, TRACKS=2};

// Some variables, like the force of infection, are indexed by sex and
// circumcision status. The index MU indicates uncircumcised men, MC
// indicates circumcised men. SEXES is the number of distinct sexes,
// SEXCIRC is the number of distinct sex and circumcision states.
enum Sex {F=0, M=1, SEXES=2, MU=1, MC=2, SEXCIRC=3};

// Sexual activity levels
enum Level {LEAST=0, LOW=1, MID=2, HIGH=3, LEVELS=4};

// Population is divided into 5-year age bands, except that the
// population aged >=55 is stored in a single age band (a=40, b=8)
// Age 15 corresponds to a=0, age 54 corresponds to a=39
enum Age {BANDS=9, AGES=41};

// b = band[a] sets b to the age band corresponding to discrete age a.
const int band[] = {
  0,0,0,0,0,
  1,1,1,1,1,
  2,2,2,2,2,
  3,3,3,3,3,
  4,4,4,4,4,
  5,5,5,5,5,
  6,6,6,6,6,
  7,7,7,7,7,
  8
};

// P: Primary (acute) infection
// L: Latent (chronic) infection
// A: AIDS
enum Stage {P1=0, P2=1, L1=2, L2=3, L3=4, A1=5, STAGES=6};

// NONE: not on ART or PrEP
// PREP_HIGH: on PrEP with high drug concentrations
// PREP_POOR: on PrEP with low drug concentrations
// ART1_NEW: on first-line ART, first year
// ART1_OLD: on first-line ART, subsequent years
// ART1_NA: non-adherent to first-line ART
// ART2: on second-line ART
// ART2_NA: non-adherent to second-line ART
enum Drug {NONE=0, PREP_HIGH=1, PREP_POOR=2, PREP_STATES=3, ART1_NEW=3, ART1_OLD=4, ART1_NA=5, ART2=6, ART2_NA=7, ARV_STATES=8};

// By design, variants 0...4 may be present at transmission. The
// implementation may take advantage of this ordering.
// R: DR from ART exposure, PrEP-sensitive 
// C: DR from ART exposure, PrEP cross-resistant
// S2 is acquired resistance to second-line ART.
// TRANSMIT: number of transmissible variants (WT, R1, C1, Q1, S1)
enum Virus {WT=0,
	    R1=1,   C1=2,   Q1=3,   S1=4,
	    R2=5,   C2=6,   Q2=7,   S2=8,
	    WR1= 9, WC1=10, WQ1=11, WS1=12,
	    WR2=13, WC2=14, WQ2=15, WS2=16,
	    VARIANTS=17, TRANSMIT=5};

// HIV disease state consists of disease stage (X susceptible,
// Y1,...,Y6 infected), ARV status (N,H,P,E,L,F,T) and virus variant
// (WT, {R,C,Q,S,r,c,q}x{1,2}+s1). Acquired resistance to first-line
// ART is further stratified by second-line susceptibility. Since
// people of any sex, age or activity level can be in all of these
// states, there's no need to substratify further
enum State {
  // susceptible, N: off PrEP, H: high PrEP drug levels, P: poor PrEP drug levels
  XN=0, XH, XP,

  // infected, N: not on ARVs (S2/s2 represented seperately below)
  Y1NWT, Y2NWT, Y3NWT, Y4NWT, Y5NWT, Y6NWT, // WT
  Y1NR1, Y2NR1, Y3NR1, Y4NR1, Y5NR1, Y6NR1, // R1
  Y1NC1, Y2NC1, Y3NC1, Y4NC1, Y5NC1, Y6NC1, // C1
  Y1NQ1, Y2NQ1, Y3NQ1, Y4NQ1, Y5NQ1, Y6NQ1, // Q1
  Y1NS1, Y2NS1, Y3NS1, Y4NS1, Y5NS1, Y6NS1, // S1
  Y4NR2S, Y5NR2S, Y6NR2S, // R2, sensitive to second-line ART
  Y4NC2S, Y5NC2S, Y6NC2S, // C2, sensitive to second-line ART
  Y4NR2R, Y5NR2R, Y6NR2R, // R2, resistant to second-line ART
  Y4NC2R, Y5NC2R, Y6NC2R, // C2, resistant to second-line ART
  Y1NQ2, Y2NQ2, Y3NQ2, Y4NQ2, Y5NQ2, Y6NQ2, // Q2
  Y1Nr1, Y2Nr1, Y3Nr1, Y4Nr1, Y5Nr1, Y6Nr1, // r1
  Y1Nc1, Y2Nc1, Y3Nc1, Y4Nc1, Y5Nc1, Y6Nc1, // c1
  Y1Nq1, Y2Nq1, Y3Nq1, Y4Nq1, Y5Nq1, Y6Nq1, // q1
  Y1Ns1, Y2Ns1, Y3Ns1, Y4Ns1, Y5Ns1, Y6Ns1, // s1
  Y4Nr2S, Y5Nr2S, Y6Nr2S, // r2, sensitive to second-line ART
  Y4Nc2S, Y5Nc2S, Y6Nc2S, // c2, sensitive to second-line ART
  Y4Nr2R, Y5Nr2R, Y6Nr2R, // r2, resistant to second-line ART
  Y4Nc2R, Y5Nc2R, Y6Nc2R, // c2, resistant to second-line ART
  Y1Nq2, Y2Nq2, Y3Nq2, Y4Nq2, Y5Nq2, Y6Nq2, // q2

  // infected H: high PrEP levels [R2, C2, r2, c2, s2 ineligible]
  Y1HWT, Y2HWT, Y3HWT, Y4HWT, Y5HWT, Y6HWT, // WT
  Y1HR1, Y2HR1, Y3HR1, Y4HR1, Y5HR1, Y6HR1, // R1
  Y1HC1, Y2HC1, Y3HC1, Y4HC1, Y5HC1, Y6HC1, // C1
  Y1HQ1, Y2HQ1, Y3HQ1, Y4HQ1, Y5HQ1, Y6HQ1, // Q1
  Y1HS1, Y2HS1, Y3HS1, Y4HS1, Y5HS1, Y6HS1, // S1
  Y1HQ2, Y2HQ2, Y3HQ2, Y4HQ2, Y5HQ2, Y6HQ2, // Q2
  Y1Hr1, Y2Hr1, Y3Hr1, Y4Hr1, Y5Hr1, Y6Hr1, // r1
  Y1Hc1, Y2Hc1, Y3Hc1, Y4Hc1, Y5Hc1, Y6Hc1, // c1
  Y1Hq1, Y2Hq1, Y3Hq1, Y4Hq1, Y5Hq1, Y6Hq1, // q1
  Y1Hs1, Y2Hs1, Y3Hs1, Y4Hs1, Y5Hs1, Y6Hs1, // s1

  // infected P: poor PrEP levels [R2, C2, r2, c2, s2 ineligible]
  Y1PWT, Y2PWT, Y3PWT, Y4PWT, Y5PWT, Y6PWT, // WT
  Y1PR1, Y2PR1, Y3PR1, Y4PR1, Y5PR1, Y6PR1, // R1
  Y1PC1, Y2PC1, Y3PC1, Y4PC1, Y5PC1, Y6PC1, // C1
  Y1PQ1, Y2PQ1, Y3PQ1, Y4PQ1, Y5PQ1, Y6PQ1, // Q1
  Y1PS1, Y2PS1, Y3PS1, Y4PS1, Y5PS1, Y6PS1, // S1
  Y1PQ2, Y2PQ2, Y3PQ2, Y4PQ2, Y5PQ2, Y6PQ2, // Q2
  Y1Pr1, Y2Pr1, Y3Pr1, Y4Pr1, Y5Pr1, Y6Pr1, // r1
  Y1Pc1, Y2Pc1, Y3Pc1, Y4Pc1, Y5Pc1, Y6Pc1, // c1
  Y1Pq1, Y2Pq1, Y3Pq1, Y4Pq1, Y5Pq1, Y6Pq1, // q1
  Y1Ps1, Y2Ps1, Y3Ps1, Y4Ps1, Y5Ps1, Y6Ps1, // s1

  // infected E: early (first-year) ART
  Y4EWT, Y5EWT, Y6EWT, // WT
  Y4ER1, Y5ER1, Y6ER1, Y4Er1, Y5Er1, Y6Er1, // R1 and r1
  Y4EC1, Y5EC1, Y6EC1, Y4Ec1, Y5Ec1, Y6Ec1, // C1 and c1
  Y4EQ1, Y5EQ1, Y6EQ1, Y4Eq1, Y5Eq1, Y6Eq1, // Q1 and q1
  Y4ES1, Y5ES1, Y6ES1, Y4Es1, Y5Es1, Y6Es1, // S1 and s1
  Y4EQ2, Y5EQ2, Y6EQ2, Y4Eq2, Y5Eq2, Y6Eq2, // Q2 and q2

  // infected L: late (subsequent years) ART
  Y4LWT, Y5LWT, Y6LWT, // WT
  Y4LR1, Y5LR1, Y6LR1, Y4Lr1, Y5Lr1, Y6Lr1, // R1 and r1
  Y4LC1, Y5LC1, Y6LC1, Y4Lc1, Y5Lc1, Y6Lc1, // C1 and c1
  Y4LQ1, Y5LQ1, Y6LQ1, Y4Lq1, Y5Lq1, Y6Lq1, // Q1 and q1
  Y4LS1, Y5LS1, Y6LS1, Y4Ls1, Y5Ls1, Y6Ls1, // S1 and s1
  Y4LR2S, Y5LR2S, Y6LR2S, // R2, sensitive to second-line ART
  Y4LC2S, Y5LC2S, Y6LC2S, // C2, sensitive to second-line ART
  Y4LR2R, Y5LR2R, Y6LR2R, // R2, resistant to second-line ART
  Y4LC2R, Y5LC2R, Y6LC2R, // C2, resistant to second-line ART
  Y4LQ2, Y5LQ2, Y6LQ2, Y4Lq2, Y5Lq2, Y6Lq2, // Q2 and q2

  // infected F: Failed ART due to non-adherence
  Y4FWT, Y5FWT, Y6FWT, // WT
  Y4FR1, Y5FR1, Y6FR1, Y4Fr1, Y5Fr1, Y6Fr1, // R1 and r1
  Y4FC1, Y5FC1, Y6FC1, Y4Fc1, Y5Fc1, Y6Fc1, // C1 and c1
  Y4FQ1, Y5FQ1, Y6FQ1, Y4Fq1, Y5Fq1, Y6Fq1, // Q1 and q1
  Y4FS1, Y5FS1, Y6FS1, Y4Fs1, Y5Fs1, Y6Fs1, // S1 and s1
  Y4FQ2, Y5FQ2, Y6FQ2, Y4Fq2, Y5Fq2, Y6Fq2, // Q2 and q2

  // infected T: second-line ART
  Y4TWT, Y5TWT, Y6TWT, // WT
  Y4TR1, Y5TR1, Y6TR1, Y4Tr1, Y5Tr1, Y6Tr1, // R1 and r1
  Y4TC1, Y5TC1, Y6TC1, Y4Tc1, Y5Tc1, Y6Tc1, // C1 and c1
  Y4TQ1, Y5TQ1, Y6TQ1, Y4Tq1, Y5Tq1, Y6Tq1, // Q1 and q1
  Y4TS1, Y5TS1, Y6TS1, Y4Ts1, Y5Ts1, Y6Ts1, // S1 and s1
  Y4TQ2, Y5TQ2, Y6TQ2, Y4Tq2, Y5Tq2, Y6Tq2, // Q2 and q2
  Y4TR2S, Y5TR2S, Y6TR2S, // R2, sensitive to second-line ART
  Y4TC2S, Y5TC2S, Y6TC2S, // C2, sensitive to second-line ART
  Y4TR2R, Y5TR2R, Y6TR2R, // R2, TDR to second-line ART
  Y4TC2R, Y5TC2R, Y6TC2R, // C2, TDR to second-line ART
  Y4Tr2S, Y5Tr2S, Y6Tr2S, // r2, sensitive to second-line ART
  Y4Tc2S, Y5Tc2S, Y6Tc2S, // c2, sensitive to second-line ART
  Y4Tr2R, Y5Tr2R, Y6Tr2R, // r2, TDR to second-line ART
  Y4Tc2R, Y5Tc2R, Y6Tc2R, // c2, TDR to second-line ART
  Y4TS2, Y5TS2, Y6TS2,    // S2 (acquired resistance to second-line ART)

  // infected D: dropout after second-line switch but before
  // resistance to second-line acquired. For example, Y4NR2S is
  // distinct from Y4DR2S in that ART uptake in the former puts
  // individuals on first-line while ART uptake in the latter puts
  // individuals on second-line.
  Y4DWT, Y5DWT, Y6DWT, // WT
  Y4DR1, Y5DR1, Y6DR1, Y4Dr1, Y5Dr1, Y6Dr1, // R1 and r1
  Y4DC1, Y5DC1, Y6DC1, Y4Dc1, Y5Dc1, Y6Dc1, // C1 and c1
  Y4DQ1, Y5DQ1, Y6DQ1, Y4Dq1, Y5Dq1, Y6Dq1, // Q1 and q1
  Y4DS1, Y5DS1, Y6DS1, Y4Ds1, Y5Ds1, Y6Ds1, // S1 and s1
  Y4DQ2, Y5DQ2, Y6DQ2, Y4Dq2, Y5Dq2, Y6Dq2, // Q2 and q2
  Y4DR2S, Y5DR2S, Y6DR2S, // R2, sensitive to second-line ART
  Y4DC2S, Y5DC2S, Y6DC2S, // C2, sensitive to second-line ART
  Y4DR2R, Y5DR2R, Y6DR2R, // R2, TDR to second-line ART
  Y4DC2R, Y5DC2R, Y6DC2R, // C2, TDR to second-line ART
  Y4Dr2S, Y5Dr2S, Y6Dr2S, // r2, sensitive to second-line ART
  Y4Dc2S, Y5Dc2S, Y6Dc2S, // c2, sensitive to second-line ART
  Y4Dr2R, Y5Dr2R, Y6Dr2R, // r2, TDR to second-line ART
  Y4Dc2R, Y5Dc2R, Y6Dc2R, // c2, TDR to second-line ART
  Y4DS2, Y5DS2, Y6DS2, // S2
  Y4Ds2, Y5Ds2, Y6Ds2, // s2

  // infected U: non-adherent to second-line ART
  Y4UWT, Y5UWT, Y6UWT, // WT
  Y4UR1, Y5UR1, Y6UR1, Y4Ur1, Y5Ur1, Y6Ur1, // R1 and r1
  Y4UC1, Y5UC1, Y6UC1, Y4Uc1, Y5Uc1, Y6Uc1, // C1 and c1
  Y4UQ1, Y5UQ1, Y6UQ1, Y4Uq1, Y5Uq1, Y6Uq1, // Q1 and q1
  Y4US1, Y5US1, Y6US1, Y4Us1, Y5Us1, Y6Us1, // S1 and s1
  Y4UQ2, Y5UQ2, Y6UQ2, Y4Uq2, Y5Uq2, Y6Uq2, // Q2 and q2
  Y4UR2S, Y5UR2S, Y6UR2S, // R2, sensitive to second-line ART
  Y4UC2S, Y5UC2S, Y6UC2S, // C2, sensitive to second-line ART
  Y4UR2R, Y5UR2R, Y6UR2R, // R2, resistant to second-line ART
  Y4UC2R, Y5UC2R, Y6UC2R, // C2, resistant to second-line ART
  Y4Ur2S, Y5Ur2S, Y6Ur2S, // r2, sensitive to second-line ART
  Y4Uc2S, Y5Uc2S, Y6Uc2S, // c2, sensitive to second-line ART
  Y4Ur2R, Y5Ur2R, Y6Ur2R, // r2, resistant to second-line ART
  Y4Uc2R, Y5Uc2R, Y6Uc2R, // c2, resistant to second-line ART

  STATES=504
};

// indices of special compartments
enum Index {
  ACTIVE = SEXCIRC * STATES * LEVELS * AGES * TRACKS, // number of compartments of living people

  // indices of compartments used to store cumulative outcomes from simulation start.

  // Indices prefixed with ACTIVE_ record cumulative outcomes in the
  // population aged 15-54.

  // indices for life-years lived (LYL) by stage of infection
  ACTIVE_LYL_X = SEXCIRC * STATES * LEVELS * AGES * TRACKS, // susceptible, men+women
  ACTIVE_W_LYL_X, ACTIVE_M_LYL_X, // women (W), men (M)
  ACTIVE_LYL_P1, ACTIVE_W_LYL_P1, ACTIVE_M_LYL_P1,
  ACTIVE_LYL_P2, ACTIVE_W_LYL_P2, ACTIVE_M_LYL_P2,
  ACTIVE_LYL_L1, ACTIVE_W_LYL_L1, ACTIVE_M_LYL_L1,
  ACTIVE_LYL_L2, ACTIVE_W_LYL_L2, ACTIVE_M_LYL_L2,
  ACTIVE_LYL_L3, ACTIVE_W_LYL_L3, ACTIVE_M_LYL_L3,
  ACTIVE_LYL_A1, ACTIVE_W_LYL_A1, ACTIVE_M_LYL_A1,

  // indices of compartments used to store cumulative new infections
  ACTIVE_WT, ACTIVE_W_WT, ACTIVE_M_WT, // wild-type
  ACTIVE_R1, ACTIVE_W_R1, ACTIVE_M_R1, // TDR to ART without PrEP cross resistance
  ACTIVE_C1, ACTIVE_W_C1, ACTIVE_M_C1, // TDR to ART with PrEP cross resistance
  ACTIVE_Q1, ACTIVE_W_Q1, ACTIVE_M_Q1, // TDR to PrEP
  ACTIVE_S1, ACTIVE_W_S1, ACTIVE_M_S1, // TDR to 2nd-line ART

  // indices for cumulative mortality
  ACTIVE_HIV_MORT, ACTIVE_W_HIV_MORT, ACTIVE_M_HIV_MORT,
  ACTIVE_BAK_MORT, ACTIVE_W_BAK_MORT, ACTIVE_M_BAK_MORT,

  // indices for person-years of ART by stage of infection. For
  // forwards-compatability, compartments are defined for stages prior
  // to ART eligibility
  ACTIVE_ART_P1, ACTIVE_W_ART_P1, ACTIVE_M_ART_P1,
  ACTIVE_ART_P2, ACTIVE_W_ART_P2, ACTIVE_M_ART_P2,
  ACTIVE_ART_L1, ACTIVE_W_ART_L1, ACTIVE_M_ART_L1,
  ACTIVE_ART_L2, ACTIVE_W_ART_L2, ACTIVE_M_ART_L2,
  ACTIVE_ART_L3, ACTIVE_W_ART_L3, ACTIVE_M_ART_L3,
  ACTIVE_ART_A1, ACTIVE_W_ART_A1, ACTIVE_M_ART_A1,

  // second-line ART
  ACTIVE_2ND_P1, ACTIVE_W_2ND_P1, ACTIVE_M_2ND_P1,
  ACTIVE_2ND_P2, ACTIVE_W_2ND_P2, ACTIVE_M_2ND_P2,
  ACTIVE_2ND_L1, ACTIVE_W_2ND_L1, ACTIVE_M_2ND_L1,
  ACTIVE_2ND_L2, ACTIVE_W_2ND_L2, ACTIVE_M_2ND_L2,
  ACTIVE_2ND_L3, ACTIVE_W_2ND_L3, ACTIVE_M_2ND_L3,
  ACTIVE_2ND_A1, ACTIVE_W_2ND_A1, ACTIVE_M_2ND_A1,

  // Cumulative ART initiation
  ACTIVE_AINIT, ACTIVE_W_AINIT, ACTIVE_M_AINIT,

  // index for person-years of PrEP
  ACTIVE_PRP, ACTIVE_W_PRP, ACTIVE_M_PRP,

  // Cumulative male circumcisions
  ACTIVE_MMC_DEBUT,  // number of men who debut circumcised 
  ACTIVE_MMC_UPTAKE, // number of men circumcised after debut

  // Indices prefixed with COHORT_ record cumulative outcomes among
  // all people who debut before PrEP closure (parameter
  // p.prep_time_closure). Otherwise, these follow the same naming
  // convention as above.

  COHORT_LYL_X,  COHORT_W_LYL_X,  COHORT_M_LYL_X,
  COHORT_LYL_P1, COHORT_W_LYL_P1, COHORT_M_LYL_P1,
  COHORT_LYL_P2, COHORT_W_LYL_P2, COHORT_M_LYL_P2,
  COHORT_LYL_L1, COHORT_W_LYL_L1, COHORT_M_LYL_L1,
  COHORT_LYL_L2, COHORT_W_LYL_L2, COHORT_M_LYL_L2,
  COHORT_LYL_L3, COHORT_W_LYL_L3, COHORT_M_LYL_L3,
  COHORT_LYL_A1, COHORT_W_LYL_A1, COHORT_M_LYL_A1,

  COHORT_WT, COHORT_W_WT, COHORT_M_WT,
  COHORT_R1, COHORT_W_R1, COHORT_M_R1,
  COHORT_C1, COHORT_W_C1, COHORT_M_C1,
  COHORT_Q1, COHORT_W_Q1, COHORT_M_Q1,
  COHORT_S1, COHORT_W_S1, COHORT_M_S1,

  COHORT_HIV_MORT, COHORT_W_HIV_MORT, COHORT_M_HIV_MORT,
  COHORT_BAK_MORT, COHORT_W_BAK_MORT, COHORT_M_BAK_MORT,

  COHORT_ART_P1, COHORT_W_ART_P1, COHORT_M_ART_P1,
  COHORT_ART_P2, COHORT_W_ART_P2, COHORT_M_ART_P2,
  COHORT_ART_L1, COHORT_W_ART_L1, COHORT_M_ART_L1,
  COHORT_ART_L2, COHORT_W_ART_L2, COHORT_M_ART_L2,
  COHORT_ART_L3, COHORT_W_ART_L3, COHORT_M_ART_L3,
  COHORT_ART_A1, COHORT_W_ART_A1, COHORT_M_ART_A1,

  COHORT_2ND_P1, COHORT_W_2ND_P1, COHORT_M_2ND_P1,
  COHORT_2ND_P2, COHORT_W_2ND_P2, COHORT_M_2ND_P2,
  COHORT_2ND_L1, COHORT_W_2ND_L1, COHORT_M_2ND_L1,
  COHORT_2ND_L2, COHORT_W_2ND_L2, COHORT_M_2ND_L2,
  COHORT_2ND_L3, COHORT_W_2ND_L3, COHORT_M_2ND_L3,
  COHORT_2ND_A1, COHORT_W_2ND_A1, COHORT_M_2ND_A1,

  COHORT_AINIT, COHORT_W_AINIT, COHORT_M_AINIT,

  COHORT_PRP, COHORT_W_PRP, COHORT_M_PRP,

  COHORT_MMC_DEBUT,
  COHORT_MMC_UPTAKE,

  // Below are compartments that store DISCOUNTED cumulative outcomes
  // indices for life-years lived (LYL) by stage of infection
  ACTDSC_LYL_X,  ACTDSC_W_LYL_X,  ACTDSC_M_LYL_X, // women (W), men (M)
  ACTDSC_LYL_P1, ACTDSC_W_LYL_P1, ACTDSC_M_LYL_P1,
  ACTDSC_LYL_P2, ACTDSC_W_LYL_P2, ACTDSC_M_LYL_P2,
  ACTDSC_LYL_L1, ACTDSC_W_LYL_L1, ACTDSC_M_LYL_L1,
  ACTDSC_LYL_L2, ACTDSC_W_LYL_L2, ACTDSC_M_LYL_L2,
  ACTDSC_LYL_L3, ACTDSC_W_LYL_L3, ACTDSC_M_LYL_L3,
  ACTDSC_LYL_A1, ACTDSC_W_LYL_A1, ACTDSC_M_LYL_A1,

  // indices of compartments used to store cumulative new infections
  ACTDSC_WT, ACTDSC_W_WT, ACTDSC_M_WT, // wild-type
  ACTDSC_R1, ACTDSC_W_R1, ACTDSC_M_R1, // TDR to ART without PrEP cross resistance
  ACTDSC_C1, ACTDSC_W_C1, ACTDSC_M_C1, // TDR to ART with PrEP cross resistance
  ACTDSC_Q1, ACTDSC_W_Q1, ACTDSC_M_Q1, // TDR to PrEP
  ACTDSC_S1, ACTDSC_W_S1, ACTDSC_M_S1, // TDR to 2nd-line ART

  // indices for cumulative mortality
  ACTDSC_HIV_MORT, ACTDSC_W_HIV_MORT, ACTDSC_M_HIV_MORT,
  ACTDSC_BAK_MORT, ACTDSC_W_BAK_MORT, ACTDSC_M_BAK_MORT,

  // indices for person-years of ART by stage of infection. For
  // forwards-compatability, compartments are defined for stages prior
  // to ART eligibility
  ACTDSC_ART_P1, ACTDSC_W_ART_P1, ACTDSC_M_ART_P1,
  ACTDSC_ART_P2, ACTDSC_W_ART_P2, ACTDSC_M_ART_P2,
  ACTDSC_ART_L1, ACTDSC_W_ART_L1, ACTDSC_M_ART_L1,
  ACTDSC_ART_L2, ACTDSC_W_ART_L2, ACTDSC_M_ART_L2,
  ACTDSC_ART_L3, ACTDSC_W_ART_L3, ACTDSC_M_ART_L3,
  ACTDSC_ART_A1, ACTDSC_W_ART_A1, ACTDSC_M_ART_A1,

  ACTDSC_2ND_P1, ACTDSC_W_2ND_P1, ACTDSC_M_2ND_P1,
  ACTDSC_2ND_P2, ACTDSC_W_2ND_P2, ACTDSC_M_2ND_P2,
  ACTDSC_2ND_L1, ACTDSC_W_2ND_L1, ACTDSC_M_2ND_L1,
  ACTDSC_2ND_L2, ACTDSC_W_2ND_L2, ACTDSC_M_2ND_L2,
  ACTDSC_2ND_L3, ACTDSC_W_2ND_L3, ACTDSC_M_2ND_L3,
  ACTDSC_2ND_A1, ACTDSC_W_2ND_A1, ACTDSC_M_2ND_A1,

  // Cumulative ART initiation
  ACTDSC_AINIT, ACTDSC_W_AINIT, ACTDSC_M_AINIT,

  // index for person-years of PrEP
  ACTDSC_PRP, ACTDSC_W_PRP, ACTDSC_M_PRP,

  // Cumulative male circumcisions
  ACTDSC_MMC_DEBUT,  // number of men who debut circumcised 
  ACTDSC_MMC_UPTAKE, // number of men circumcised after debut

  COHDSC_LYL_X,  COHDSC_W_LYL_X,  COHDSC_M_LYL_X,
  COHDSC_LYL_P1, COHDSC_W_LYL_P1, COHDSC_M_LYL_P1,
  COHDSC_LYL_P2, COHDSC_W_LYL_P2, COHDSC_M_LYL_P2,
  COHDSC_LYL_L1, COHDSC_W_LYL_L1, COHDSC_M_LYL_L1,
  COHDSC_LYL_L2, COHDSC_W_LYL_L2, COHDSC_M_LYL_L2,
  COHDSC_LYL_L3, COHDSC_W_LYL_L3, COHDSC_M_LYL_L3,
  COHDSC_LYL_A1, COHDSC_W_LYL_A1, COHDSC_M_LYL_A1,

  COHDSC_WT, COHDSC_W_WT, COHDSC_M_WT,
  COHDSC_R1, COHDSC_W_R1, COHDSC_M_R1,
  COHDSC_C1, COHDSC_W_C1, COHDSC_M_C1,
  COHDSC_Q1, COHDSC_W_Q1, COHDSC_M_Q1,
  COHDSC_S1, COHDSC_W_S1, COHDSC_M_S1,

  COHDSC_HIV_MORT, COHDSC_W_HIV_MORT, COHDSC_M_HIV_MORT,
  COHDSC_BAK_MORT, COHDSC_W_BAK_MORT, COHDSC_M_BAK_MORT,

  COHDSC_ART_P1, COHDSC_W_ART_P1, COHDSC_M_ART_P1,
  COHDSC_ART_P2, COHDSC_W_ART_P2, COHDSC_M_ART_P2,
  COHDSC_ART_L1, COHDSC_W_ART_L1, COHDSC_M_ART_L1,
  COHDSC_ART_L2, COHDSC_W_ART_L2, COHDSC_M_ART_L2,
  COHDSC_ART_L3, COHDSC_W_ART_L3, COHDSC_M_ART_L3,
  COHDSC_ART_A1, COHDSC_W_ART_A1, COHDSC_M_ART_A1,

  COHDSC_2ND_P1, COHDSC_W_2ND_P1, COHDSC_M_2ND_P1,
  COHDSC_2ND_P2, COHDSC_W_2ND_P2, COHDSC_M_2ND_P2,
  COHDSC_2ND_L1, COHDSC_W_2ND_L1, COHDSC_M_2ND_L1,
  COHDSC_2ND_L2, COHDSC_W_2ND_L2, COHDSC_M_2ND_L2,
  COHDSC_2ND_L3, COHDSC_W_2ND_L3, COHDSC_M_2ND_L3,
  COHDSC_2ND_A1, COHDSC_W_2ND_A1, COHDSC_M_2ND_A1,

  COHDSC_AINIT, COHDSC_W_AINIT, COHDSC_M_AINIT,

  COHDSC_PRP, COHDSC_W_PRP, COHDSC_M_PRP,

  COHDSC_MMC_DEBUT,
  COHDSC_MMC_UPTAKE,

  COMPARTMENTS // number of compartments
};

// used to map stages to cumulative outcomes compartments
const Index active_lyl[] = {ACTIVE_LYL_P1, ACTIVE_LYL_P2, ACTIVE_LYL_L1, ACTIVE_LYL_L2, ACTIVE_LYL_L3, ACTIVE_LYL_A1};
const Index active_art[] = {ACTIVE_ART_P1, ACTIVE_ART_P2, ACTIVE_ART_L1, ACTIVE_ART_L2, ACTIVE_ART_L3, ACTIVE_ART_A1};
const Index active_2nd[] = {ACTIVE_2ND_P1, ACTIVE_2ND_P2, ACTIVE_2ND_L1, ACTIVE_2ND_L2, ACTIVE_2ND_L3, ACTIVE_2ND_A1};
const Index cohort_lyl[] = {COHORT_LYL_P1, COHORT_LYL_P2, COHORT_LYL_L1, COHORT_LYL_L2, COHORT_LYL_L3, COHORT_LYL_A1};
const Index cohort_art[] = {COHORT_ART_P1, COHORT_ART_P2, COHORT_ART_L1, COHORT_ART_L2, COHORT_ART_L3, COHORT_ART_A1};
const Index cohort_2nd[] = {COHORT_2ND_P1, COHORT_2ND_P2, COHORT_2ND_L1, COHORT_2ND_L2, COHORT_2ND_L3, COHORT_2ND_A1};

const Index active_w_lyl[] = {ACTIVE_W_LYL_P1, ACTIVE_W_LYL_P2, ACTIVE_W_LYL_L1, ACTIVE_W_LYL_L2, ACTIVE_W_LYL_L3, ACTIVE_W_LYL_A1};
const Index active_w_art[] = {ACTIVE_W_ART_P1, ACTIVE_W_ART_P2, ACTIVE_W_ART_L1, ACTIVE_W_ART_L2, ACTIVE_W_ART_L3, ACTIVE_W_ART_A1};
const Index active_w_2nd[] = {ACTIVE_W_2ND_P1, ACTIVE_W_2ND_P2, ACTIVE_W_2ND_L1, ACTIVE_W_2ND_L2, ACTIVE_W_2ND_L3, ACTIVE_W_2ND_A1};
const Index cohort_w_lyl[] = {COHORT_W_LYL_P1, COHORT_W_LYL_P2, COHORT_W_LYL_L1, COHORT_W_LYL_L2, COHORT_W_LYL_L3, COHORT_W_LYL_A1};
const Index cohort_w_art[] = {COHORT_W_ART_P1, COHORT_W_ART_P2, COHORT_W_ART_L1, COHORT_W_ART_L2, COHORT_W_ART_L3, COHORT_W_ART_A1};
const Index cohort_w_2nd[] = {COHORT_W_2ND_P1, COHORT_W_2ND_P2, COHORT_W_2ND_L1, COHORT_W_2ND_L2, COHORT_W_2ND_L3, COHORT_W_2ND_A1};

const Index active_m_lyl[] = {ACTIVE_M_LYL_P1, ACTIVE_M_LYL_P2, ACTIVE_M_LYL_L1, ACTIVE_M_LYL_L2, ACTIVE_M_LYL_L3, ACTIVE_M_LYL_A1};
const Index active_m_art[] = {ACTIVE_M_ART_P1, ACTIVE_M_ART_P2, ACTIVE_M_ART_L1, ACTIVE_M_ART_L2, ACTIVE_M_ART_L3, ACTIVE_M_ART_A1};
const Index active_m_2nd[] = {ACTIVE_M_2ND_P1, ACTIVE_M_2ND_P2, ACTIVE_M_2ND_L1, ACTIVE_M_2ND_L2, ACTIVE_M_2ND_L3, ACTIVE_M_2ND_A1};
const Index cohort_m_lyl[] = {COHORT_M_LYL_P1, COHORT_M_LYL_P2, COHORT_M_LYL_L1, COHORT_M_LYL_L2, COHORT_M_LYL_L3, COHORT_M_LYL_A1};
const Index cohort_m_art[] = {COHORT_M_ART_P1, COHORT_M_ART_P2, COHORT_M_ART_L1, COHORT_M_ART_L2, COHORT_M_ART_L3, COHORT_M_ART_A1};
const Index cohort_m_2nd[] = {COHORT_M_2ND_P1, COHORT_M_2ND_P2, COHORT_M_2ND_L1, COHORT_M_2ND_L2, COHORT_M_2ND_L3, COHORT_M_2ND_A1};

const Index actdsc_lyl[] = {ACTDSC_LYL_P1, ACTDSC_LYL_P2, ACTDSC_LYL_L1, ACTDSC_LYL_L2, ACTDSC_LYL_L3, ACTDSC_LYL_A1};
const Index actdsc_art[] = {ACTDSC_ART_P1, ACTDSC_ART_P2, ACTDSC_ART_L1, ACTDSC_ART_L2, ACTDSC_ART_L3, ACTDSC_ART_A1};
const Index actdsc_2nd[] = {ACTDSC_2ND_P1, ACTDSC_2ND_P2, ACTDSC_2ND_L1, ACTDSC_2ND_L2, ACTDSC_2ND_L3, ACTDSC_2ND_A1};
const Index cohdsc_lyl[] = {COHDSC_LYL_P1, COHDSC_LYL_P2, COHDSC_LYL_L1, COHDSC_LYL_L2, COHDSC_LYL_L3, COHDSC_LYL_A1};
const Index cohdsc_art[] = {COHDSC_ART_P1, COHDSC_ART_P2, COHDSC_ART_L1, COHDSC_ART_L2, COHDSC_ART_L3, COHDSC_ART_A1};
const Index cohdsc_2nd[] = {COHDSC_2ND_P1, COHDSC_2ND_P2, COHDSC_2ND_L1, COHDSC_2ND_L2, COHDSC_2ND_L3, COHDSC_2ND_A1};

const Index actdsc_w_lyl[] = {ACTDSC_W_LYL_P1, ACTDSC_W_LYL_P2, ACTDSC_W_LYL_L1, ACTDSC_W_LYL_L2, ACTDSC_W_LYL_L3, ACTDSC_W_LYL_A1};
const Index actdsc_w_art[] = {ACTDSC_W_ART_P1, ACTDSC_W_ART_P2, ACTDSC_W_ART_L1, ACTDSC_W_ART_L2, ACTDSC_W_ART_L3, ACTDSC_W_ART_A1};
const Index actdsc_w_2nd[] = {ACTDSC_W_2ND_P1, ACTDSC_W_2ND_P2, ACTDSC_W_2ND_L1, ACTDSC_W_2ND_L2, ACTDSC_W_2ND_L3, ACTDSC_W_2ND_A1};
const Index cohdsc_w_lyl[] = {COHDSC_W_LYL_P1, COHDSC_W_LYL_P2, COHDSC_W_LYL_L1, COHDSC_W_LYL_L2, COHDSC_W_LYL_L3, COHDSC_W_LYL_A1};
const Index cohdsc_w_art[] = {COHDSC_W_ART_P1, COHDSC_W_ART_P2, COHDSC_W_ART_L1, COHDSC_W_ART_L2, COHDSC_W_ART_L3, COHDSC_W_ART_A1};
const Index cohdsc_w_2nd[] = {COHDSC_W_2ND_P1, COHDSC_W_2ND_P2, COHDSC_W_2ND_L1, COHDSC_W_2ND_L2, COHDSC_W_2ND_L3, COHDSC_W_2ND_A1};

const Index actdsc_m_lyl[] = {ACTDSC_M_LYL_P1, ACTDSC_M_LYL_P2, ACTDSC_M_LYL_L1, ACTDSC_M_LYL_L2, ACTDSC_M_LYL_L3, ACTDSC_M_LYL_A1};
const Index actdsc_m_art[] = {ACTDSC_M_ART_P1, ACTDSC_M_ART_P2, ACTDSC_M_ART_L1, ACTDSC_M_ART_L2, ACTDSC_M_ART_L3, ACTDSC_M_ART_A1};
const Index actdsc_m_2nd[] = {ACTDSC_M_2ND_P1, ACTDSC_M_2ND_P2, ACTDSC_M_2ND_L1, ACTDSC_M_2ND_L2, ACTDSC_M_2ND_L3, ACTDSC_M_2ND_A1};
const Index cohdsc_m_lyl[] = {COHDSC_M_LYL_P1, COHDSC_M_LYL_P2, COHDSC_M_LYL_L1, COHDSC_M_LYL_L2, COHDSC_M_LYL_L3, COHDSC_M_LYL_A1};
const Index cohdsc_m_art[] = {COHDSC_M_ART_P1, COHDSC_M_ART_P2, COHDSC_M_ART_L1, COHDSC_M_ART_L2, COHDSC_M_ART_L3, COHDSC_M_ART_A1};
const Index cohdsc_m_2nd[] = {COHDSC_M_2ND_P1, COHDSC_M_2ND_P2, COHDSC_M_2ND_L1, COHDSC_M_2ND_L2, COHDSC_M_2ND_L3, COHDSC_M_2ND_A1};

// Replacement Schemes (REPLACE_*)
// NONE      Activity levels locked in at sexual debut. The duration
//           of sex work input parameter is ignored.
//
// STATIC    Women exit sex work, and initiate sex work at
//           age-specific, static rates.
//
// ALL       Women exit sex work. Higher activity levels depleted by
//           sex work ceasation and HIV mortality are replenished from
//           lower activity levels so that the risk profile remains
//           static throughout simulation
//
// CSW       Women exit sex work. Sex workers are replenished from
//           lower activity levels so that the proportion of women
//           who are sex workers remains static throughout
//           simulation.
//
// DEMAND    Women exit sex work. Sex workers are replenished from
//           lower activity levels so that unmet demand for sex
//           workers is minimal.
enum Replace {
  REPLACE_NONE = 1,
  REPLACE_STATIC = 2,
  REPLACE_ALL = 3,
  REPLACE_CSW = 4,
  REPLACE_DEMAND = 5,
  REPLACE_SCHEMES
};

// Attributes arrays. Used to associate each infection state above
// with its corresponding infection stage, ARV state, majority virus
// variant, and the transmitted virus variant. default values are used
// for susceptible states
extern Stage attr_stage[STATES]; // stage of infection
extern Drug attr_drug[STATES]; // ARV state
extern Virus attr_virus[STATES]; // majority variant
extern Virus attr_transmit[STATES]; // transmitted variant
extern std::string names[STATES];

struct Parameters {

  // Read demographic, behavioral and epidemiological model inputs
  int read_epidemic(std::string const& filename);

  // Read intervention and drug resistance model inputs
  int read_intervention(std::string const& filename);

  // set the infectivity of each infection state
  void set_infectivity();

  // year that HIV is introduced into the population
  unsigned int init_year;

  double init_size; // initial population size

  double rate_grow_init; // initial population growth rate
  double rate_grow_late; // asymptotic population growth rate
  double rate_grow_decr; // population growth rate decay parameter

  // distribution of sex, age and risk at simulation start
  double init_risk[SEXES][BANDS][LEVELS];

  // distribution of sex and risk at sexual debut
  double prop_risk[SEXES][LEVELS];

  // distribution of sex at debut
  double prop_sex[SEXES];

  // Indicator variable that selects how individuals move between
  // activity levels, if at all. Valid values are defined and
  // explained when enum Replace is declared above
  Replace replace;

  // rate that high risk women (commercial sex workers, or CSWs)
  // transition to lower sexual activity levels. Ignored when
  // replace=REPLACE_NONE
  double rate_csw_exit;

  // In schemes with replacement, prop_risk_exit scales the rate of
  // exit from higher activity levels based on health
  // state. prop_risk_init similarly scales the rate of entry into
  // higher activity levels.
  double prop_risk_exit[STATES];
  double prop_risk_init[STATES];

  // reduction in partner change rate by disease state
  double prop_risk_decr[STATES];

  double rate_age; // rate of aging

  // age-stratified non-HIV mortality rates
  double rate_mort[SEXES][BANDS];

  // Preferential mixing parameters
  double assort_act; // assortativity by sexual activity level
  double assort_age; // assortativity by age
  double assort_dif; // % of age-assortativity between older men and younger women

  // partner change rates
  double partner[SEXES][BANDS][LEVELS];

  // Sexual acts per partnership. The first index is the female
  // partner's activity level, the second index is the male partner's
  // activity level
  double acts[LEVELS][LEVELS];

  double infect_base;            // baseline per-act infectivity
  double infect_stage[STAGES];   // relative infectivity by stage of infection
  double infect_virus[VARIANTS]; // relative infectivity by virus variant

  // infect[c][m][p][s] stores ln(1 - p), where p is the probability
  // of HIV transmission per coital act from a donor in disease state
  // s given condom use (boolean c), MMC status of susceptible partner
  // (boolean m) and PrEP use by susceptible partner (boolean p)
  double infect[2][2][2][STATES];

  // disease progression rates
  double prog[VARIANTS][STAGES];

  // +=+ condom-related parameters +===========================================+
  // Proportion of acts in which condoms are used. The first index is
  // the female partner's activity level, the second index is the male
  // partner's activity level
  double condom_prop[LEVELS][LEVELS];

  // relative per-act infectivity when condoms are used
  double condom_efficacy;

  // changes to condom use
  double condom_change;      // change in use relative to level pre-change
  double condom_change_time; // time that change begins
  double condom_change_span; // time span that change occurs in

  // +=+ circumcision-related parameters +=====================================+
  double mmc_efficacy;     // relative infectivity for circumcised men
  double mmc_prop_init;    // proportion of men who debut circumcised initially
  double mmc_prop_late;    // proportion of men who debut circumcised after scale-up
  double mmc_time_scaleup; // time MMC scale-up begins
  double mmc_span_scaleup; // time MMC scale-up takes to complete

  // +=+ PrEP-related parameters +=============================================+
  double prep_time_rollout;          // first year of PrEP rollout
  double prep_span_rollout;          // duration of PrEP scaleup
  double prep_time_closure;          // time PrEP enrollment stops

  // Rates of PrEP initiation (prep_init), injection (prep_inject) and
  // ceasation (prep_exit) are stratified by sex, age and sexual
  // activity level to support PrEP prioritization.
  double prep_target;                       // Population-level PrEP coverage
  double prep_init[SEXES][BANDS][LEVELS];   // The proportion who start PrEP each year
  double prep_inject[SEXES][BANDS][LEVELS]; // PrEP injection rate
  double prep_exit[SEXES][BANDS][LEVELS];   // rate of exit from PrEP

  // prep_inject and prep_exit above take PrEP eligibility into
  // account. For example, 
  //
  //   prep_exit = prep_base_inject + prep_base_exit 
  //
  // for people who are not eligible for PrEP, as individuals who
  // become ineligible for PrEP are taken off PrEP if they present for
  // a follow-up injection.
  double prep_base_inject;
  double prep_base_exit;

  // HIV testing in the PrEP program. This is expressed as a fraction
  // of injections at which testing is administered (e.g., prep_test=1
  // implies testing before every injection, prep_test=0.5 implies
  // testing before 50% of injections. This only appears in equations
  // for established infection
  double prep_test;

  double prep_prop_high;             // proportion with high drug levels
  double prep_efficacy;              // PrEP efficacy with drug level = 1
  double prep_adh_high;              // PrEP adherence in high concentration stratum
  double prep_adh_poor;              // PrEP adherence in poor concentration stratum
  double prep_adr_high;              // selective pressure with high drug levels
  double prep_adr_poor;              // selective pressure with poor drug levels

  // relative efficacy against transmitted virus variants
  double prep_eff_virus[TRANSMIT];

  // +=+ ART-related parameters +==============================================+

  // For convenience, ART implementation parameters are allocated for
  // every stage of infection. Since ART use compartments are defined
  // only for stages L2, L3 and A1, half of the ART coverage
  // parameters are unnused
  double art_time_start[STAGES]; // time that ART coverage begins
  double art_time_limit[STAGES]; // time that ART coverage stabilizes
  double art_prop_start[STAGES]; // proportion starting ART annually during scale-up
  double art_prop_limit[STAGES]; // proportion starting ART annually during maintainence

  double art_advance;          // inverse of early ART duration

  double art_exitE;            // dropout rate, early ART
  double art_exitL;            // dropout rate, late ART
  double art_exitF;            // dropout rate, non-adherent to ART

  double art_failE;            // non-adherence rate, early ART
  double art_failL;            // non-adherence rate, late ART

  double art_mortE[STAGES];    // HIV-related mortality rate during early ART
  double art_mortL[STAGES];    // HIV-related mortality rate during late ART

  double art_efficacy_wt; // Relative infectivity on ART (wild-type)
  double art_efficacy_dr; // Relative infectivity on ART (pre-existing DR)

  double sec_init_dr; // second-line ART initiation after switching for drug resistance 
  double sec_init_na; // second-line ART initiation after switching for non-adherence
  double sec_exit; // second-line ART dropout
  double sec_mort[STAGES]; // HIV-related mortality on second-line ART
  double sec_fail_dr; // non-adherence rate after switching for drug resistance
  double sec_fail_na; // non-adherence rate after switching for non-adherence

  double xdr_overall; // proportion with PrEP cross resistance upon ART ADR

  // +=+ ARV resistance parameters +===========================================+
  double emerge[ARV_STATES][VARIANTS]; // resistance emergence rates
  double revert[ARV_STATES][VARIANTS]; // resistance reveresion rates

  // +=+ Cost-effectiveness parameters +=======================================+
  // Only the discount factor and horizon that quantification begins
  // are stored by the model. Users are responsible for applying costs
  // and quality weights to the model outputs. This is intended to
  // minimize the need to repeat simulations when cost assumptions
  // change.
  double ce_discount; // discount rate
  double ce_tinit;    // time that outcome accumulation begins
};

// maps a collection of attributes to the corresponding compartment index
inline int index(Sex const sex, Level const level, int const age, State const state, Track const track) {
  return TRACKS * (STATES * (AGES * (sex * LEVELS + level) + age) + state) + track;
}

// maps a collection of attributes to a risk group index
inline int risk(Sex const sex, Level const level, int const age) {
  return AGES * (sex * LEVELS + level) + age;
}

// maps a collection of attributes to a risk group index
inline int risk(int const sex, int const level, int const age) {
  return risk(static_cast<Sex>(sex), static_cast<Level>(level), age);
}

// offset into the vector for states of people with the given
// behavioral and tracking status
inline int offset(Sex const sex, Level const level, int const age, Track const track) {
  return STATES * risk(sex, level, age) + track * (STATES * SEXCIRC * AGES * LEVELS);
}

inline int offset(int const sex, int const level, int const age, int const track) {
  return offset(static_cast<Sex>(sex), static_cast<Level>(level), age, static_cast<Track>(track));
}

// initialize the ODE system state
void initialize(double y[], Parameters const& params);

// evaluate the ODE system
int model(double t, double const y[], double dy[], Parameters const& params);

// overloaded function to evaluate the ODE system. this prototype is
// used by the GSL numerical integration functions
inline int model(double t, double const y[], double dy[], void* params) {
  return model(t, y, dy, *((Parameters*)params));
}

// The first index of lambda should take three values:
// lambda[F ]: force of infection acting on women
// lambda[MU]: force of infection acting on uncircumcised men
// lambda[MC]: force of infection acting on circumcised men
//void force(double const t, double const y[], Parameters const& p, double lambda[][BANDS][LEVELS][TRANSMIT]);
void force(double const t, double const y[], Parameters const& p, double lambda[][BANDS][LEVELS][PREP_STATES][TRANSMIT]);

// Calculates rates of ART initiation. uptake is stratified by stage of
// infection. Uptake is only supported for stages L2 (350<CD4<500), L3
// (200<CD4<350) and A1 (CD4<200)
void rates_art(double const t, Parameters const& p, double uptake[]);

// Calculates rates of PrEP uptake. uptake is stratified by sex, age
// and sexual activity level
void rates_prep(double const t, double const y[], double const dy[], Parameters const& p, double uptake[][BANDS][LEVELS]);

// Helper functions used to define target PrEP coverage (mmc_cover_ft)
// and its time derivative (mmc_cover_dt)
double prep_cover_ft(double const t, Parameters const& p);
double prep_cover_dt(double const t, Parameters const& p);

// Calculates rates of male medical circumcision (MMC) uptake. uptake is stratified by age.
void rates_mmc(double const t, double const y[], double const dy[], Parameters const& p, double uptake[]);

// Helper function used to define target MMC coverage. Returns (by
// reference) the target circumcision coverage ft and the change in
// target coverage dt at time t
void mmc_cover(double const t, Parameters const& p, double& ft, double& dt);

// Calculates rates of commercial sex work initiation. Rates are
// returned in csw_init. These rates are stratified by age band and
// determined by unmet demand among men
void rates_csw_need(double const t, double const y[], double const dy[], Parameters const& p, double csw_init[]);

// Calculates rates of commercial sex work initiation. Rates are
// returned in csw_init. These rates are homogeneous across age bands,
// but csw_init has one cell per age band for compatibility with
// rates_csw_need. Rates are determined so that the proportion of
// women who are sex workers remains constant throughout simulation
void rates_csw_prop(double const t, double const y[], double const dy[], Parameters const& p, double csw_init[]);

// Calculates rates that higher activity levels are replenished from
// lower ones. Rates are stratified by sex and sexual activity level
void rates_replace(double const t, double const y[], double const dy[], Parameters const& p, double rate[][LEVELS]);

#endif // Model_H
