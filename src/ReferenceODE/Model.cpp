#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <vector>
#include <gsl/gsl_errno.h>
#include <Model.H>

using std::min;
using std::max;

int Parameters::read_epidemic(std::string const& filename) {
  std::ifstream in(filename.c_str());
  if (in.fail()) return -1;

  std::istream_iterator<double> start(in), end;
  std::vector<double> input(start, end);
  in.close();

  if (input.size() != 119) return -2;

  init_year = round(input[0]);
  init_size = input[1];
  rate_grow_init = input[2];
  rate_grow_late = input[3];
  rate_grow_decr = input[4];

  prop_sex[M] = input[5];
  prop_sex[F] = 1.0 - input[5];

  prop_risk[F][LEAST] = prop_sex[F] * input[6];
  prop_risk[F][LOW  ] = prop_sex[F] * (1 - (input[6] + input[7] + input[8]));
  prop_risk[F][MID  ] = prop_sex[F] * input[7];
  prop_risk[F][HIGH ] = prop_sex[F] * input[8];

  prop_risk[M][LEAST] = prop_sex[M] * input[9];
  prop_risk[M][LOW  ] = prop_sex[M] * (1 - (input[9] + input[10] + input[11]));
  prop_risk[M][MID  ] = prop_sex[M] * input[10];
  prop_risk[M][HIGH ] = prop_sex[M] * input[11];

  // Initial, the model population just consists of sexually active
  // adults aged 15-54. The oldest age band is assumed empty
  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS - 1; ++b) {
	init_risk[g][b][k] = input[12 + b] * prop_risk[g][k];
      }
      init_risk[g][BANDS-1][k] = 0.0;
    }
  }

  // The user may disable aging by setting the aging rate to
  // zero. This needs to be handled as a special case, since otherwise
  // rate_age is not a number, which will propagate throughout the
  // system state.
  if (input[20] == 0.0) {
    rate_age = 0.0;
  } else {
    rate_age = 1.0 / input[20];
  }

  rate_csw_exit = 1.0 / input[21];

  // Specify relative rates of behavior change based on disease stage
  // and ART status.
  for (int s(0); s < STATES; ++s) {
    switch (attr_stage[s]) {
    case L3: // 200 < CD4 < 350
      prop_risk_exit[s] = ((attr_drug[s] >= ART1_NEW) ? 2.50 : 2.50);
      prop_risk_init[s] = ((attr_drug[s] >= ART1_NEW) ? 0.40 : 0.40);
      prop_risk_decr[s] = ((attr_drug[s] >= ART1_NEW) ? 1.00 : 1.00);
      break;
    case A1: // CD4 < 200
      prop_risk_exit[s] = ((attr_drug[s] >= ART1_NEW) ? 2.50 : 4.00);
      prop_risk_init[s] = ((attr_drug[s] >= ART1_NEW) ? 0.40 : 0.25);
      prop_risk_decr[s] = ((attr_drug[s] >= ART1_NEW) ? 1.00 : 0.65);
      break;
    default:
      prop_risk_exit[s] = ((attr_drug[s] >= ART1_NEW) ? 2.50 : 1.00);
      prop_risk_init[s] = ((attr_drug[s] >= ART1_NEW) ? 0.40 : 1.00);
      prop_risk_decr[s] = ((attr_drug[s] >= ART1_NEW) ? 1.00 : 1.00);
      break;
    }
  }

#warning "Hard-coded survival at ages 55+"
  for (int b(0); b < BANDS - 1; ++b) rate_mort[F][b] = input[22 + b];
  for (int b(0); b < BANDS - 1; ++b) rate_mort[M][b] = input[30 + b];
  rate_mort[F][BANDS-1] = 0.04694; // reciprocal of WHO remaining life expectancy at age 55
  rate_mort[M][BANDS-1] = 0.05828; // reciprocal of WHO remaining life expectancy at age 55

  for (int k(0); k < LEVELS; ++k) {
    for (int b(0); b < BANDS-1; ++b) {
      partner[F][b][k] = input[38] * input[40+b] * input[56+k];
      partner[M][b][k] = input[39] * input[48+b] * input[60+k];
    }
    partner[F][BANDS-1][k] = 0.0; // Assuming sexual activity stops at age 55
    partner[M][BANDS-1][k] = 0.0; // Assuming sexual activity stops at age 55
  }

  for (int kf(0); kf < LEVELS; ++kf) {
    for (int km(0); km < LEVELS; ++km) {
      acts[kf][km] = input[64 + kf * LEVELS + km];
    }
  }

  for (int kf(0); kf < LEVELS; ++kf) {
    for (int km(0); km < LEVELS; ++km) {
      condom_prop[kf][km] = input[80 + kf * LEVELS + km];
    }
  }

  condom_efficacy = 1.0 - input[96];
  condom_change = input[97];
  condom_change_time = input[98];
  condom_change_span = input[99];

  // The parameters of MMC scale-up are set in read_intervention. The
  // values below imply MMC is not scaled up
  mmc_prop_init = input[100];
  mmc_prop_late = input[100];
  mmc_time_scaleup = 1970;
  mmc_span_scaleup = 1.0;

  mmc_efficacy = 1.0 - input[101];

  assort_act = input[102];
  assort_age = input[103];
  assort_dif = input[104];

  infect_base = input[105];
  for (int h(0); h < STAGES; ++h) infect_stage[h] = input[106 + h];

  // Infectivity of non-WT variants is set by read_intervention.
  std::fill(infect_virus, infect_virus + VARIANTS, 0.0);
  infect_virus[WT] = 1.0;
  set_infectivity();

  // fill with invalid values so that it is obvious if some required
  // parameters are not set.
  for (int v(0); v < VARIANTS; ++v) {
    std::fill(prog[v], prog[v] + STAGES, std::numeric_limits<double>::quiet_NaN());
  }

  // Progression with non-WT variants is set by read_intervention
  prog[WT][P1] = 1.0 / input[112];
  prog[WT][P2] = std::min(1.0 / input[113], 365.242);

#ifdef CD4_500
  prog[WT][L1] = 1.0 / input[114];
  prog[WT][L2] = 1.0 / input[115];
#else
  prog[WT][L1] = std::numeric_limits<double>::quiet_NaN();
  prog[WT][L2] = 1.0 / (input[114] + input[115]);
#endif // CD4_500

  prog[WT][L3] = 1.0 / input[116];
  prog[WT][A1] = 1.0 / input[117];

  replace = static_cast<Replace>(input[118]);

  if (replace == 0 || replace >= REPLACE_SCHEMES) {
    fprintf(stderr, "Warning: replacement scheme %d unrecognized, REPLACE_NONE used\n", replace);
    replace = REPLACE_NONE;
  }

  return 0;
}

int Parameters::read_intervention(std::string const& filename) {
  std::ifstream in(filename.c_str());
  if (in.fail()) return -1;
  std::istream_iterator<double> start(in), end;
  std::vector<double> input(start, end);
  in.close();
  double avfailE[VARIANTS], avfailL[VARIANTS];

  if (input.size() != 135) return -2;

  // +-+ PrEP parameters +-----------------------------------------------------+
  prep_time_rollout = input[0];
  prep_span_rollout = input[1];
  prep_time_closure = input[2];
  prep_target = input[3];

  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS-1; ++b) {
	prep_init[g][b][k] = input[4 + (g * LEVELS + k) * (BANDS - 1) + b];
      }
      prep_init[g][BANDS-1][k] = -1.0; // No PrEP uptake at ages 55+
    }
  }

  // PrEP-eligible groups have prep_init >= 0. Individuals who are not
  // eligible for PrEP are removed when they present for another
  // injection. This is implemented by adding the injection rate to
  // the dropout rate and assuming no new injections are received.
  prep_base_inject = input[68];
  //  prep_clear = input[69]; // tail not modeled in this implementation.
  prep_test = input[70];
  prep_base_exit = 1.0 / input[71];
  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS; ++b) {
	if (prep_init[g][b][k] < 0.0) {
	  prep_inject[g][b][k] = 0.0;
	  prep_exit[g][b][k] = prep_base_inject + prep_base_exit;
	} else {
	  prep_inject[g][b][k] = prep_base_inject;
	  prep_exit[g][b][k] = prep_base_exit;
	}
      }
    }
  }

  prep_prop_high = input[72];
  prep_efficacy = 1.0 - input[73];
  prep_adh_high = 1.0; // high drug levels: 100% of prep_efficacy
  prep_adh_poor = 0.0; // poor drug levels:   0% of prep_efficacy
  prep_adr_high = input[74];
  prep_adr_poor = input[75];

  prep_eff_virus[WT] = 1.0;
  prep_eff_virus[R1] = 1.0;
  prep_eff_virus[S1] = 1.0;
  prep_eff_virus[C1] = input[76];
  prep_eff_virus[Q1] = input[77];

  for (int s(0); s < ARV_STATES; ++s) {
    std::fill(emerge[s], emerge[s] + VARIANTS, 0.0);
    std::fill(revert[s], revert[s] + VARIANTS, 0.0);
  }

  const double prep_adr_base_high(-log(1.0 - prep_adr_high * 0.99) / input[78]);
  const double prep_adr_base_poor(-log(1.0 - prep_adr_poor * 0.99) / input[78]);

  emerge[PREP_HIGH][WT ] = prep_adr_base_high;
  emerge[PREP_HIGH][R1 ] = prep_adr_base_high;
  emerge[PREP_HIGH][WR1] = prep_adr_base_high;
  emerge[PREP_HIGH][S1 ] = prep_adr_base_high;
  emerge[PREP_HIGH][WS1] = prep_adr_base_high;
  emerge[PREP_HIGH][WQ1] = 2 * prep_adr_base_high;
  emerge[PREP_HIGH][WC1] = 2 * prep_adr_base_high;

  emerge[PREP_POOR][WT ] = prep_adr_base_poor;
  emerge[PREP_POOR][R1 ] = prep_adr_base_poor;
  emerge[PREP_POOR][WR1] = prep_adr_base_poor;
  emerge[PREP_POOR][S1 ] = prep_adr_base_poor;
  emerge[PREP_POOR][WS1] = prep_adr_base_poor;
  emerge[PREP_POOR][WQ1] = 2 * prep_adr_base_poor;
  emerge[PREP_POOR][WC1] = 2 * prep_adr_base_poor;

  // +-+ ART parameters +------------------------------------------------------+

  // delay ART implementation indefinitely in ART-ineligible stages
  std::fill(art_time_start, art_time_start + STAGES, std::numeric_limits<double>::infinity());
  std::fill(art_time_limit, art_time_limit + STAGES, std::numeric_limits<double>::infinity());
  art_time_start[A1] = 2004; art_time_limit[A1] = 2012;
  art_time_start[L3] = 2010; art_time_limit[L3] = 2017;
  art_prop_start[A1] = 0.1; art_prop_limit[A1] = 0.1;
  art_prop_start[L3] = 0.1; art_prop_limit[L3] = 0.1;
#ifdef CD4_500
    art_time_start[L2] = 2015; art_time_limit[L2] =2020;
    art_prop_start[L2] = 0.1;  art_prop_limit[L2] = 0.1;
#endif // CD4_500

  art_advance = 1.0;

  art_exitE = input[79];
  art_exitL = input[80];
  art_exitF = input[81];

  avfailE[WT] = input[82];
  avfailL[WT] = input[83];
  avfailE[R1] = avfailE[WR1] = avfailE[C1] = avfailE[WC1] = input[84];
  avfailL[R1] = avfailL[WR1] = avfailL[C1] = avfailL[WC1] = input[85];
  avfailE[S1] = avfailE[WS1] = input[86];
  avfailL[S1] = avfailL[WS1] = input[87];
  avfailE[Q1] = avfailE[WQ1] = input[88];
  avfailL[Q1] = avfailL[WQ1] = input[89];
  avfailE[Q2] = avfailE[WQ2] = input[90];
  avfailL[Q2] = avfailL[WQ2] = input[91];

  art_failE = input[92] * -log(1.0 - avfailE[WT]);
  art_failL = input[93] * -log(1.0 - avfailL[WT]);
  for (int v(0); v < VARIANTS; ++v) {
    emerge[ART1_NEW][v] = -log(1.0 - avfailE[v]) - art_failE;
    emerge[ART1_OLD][v] = -log(1.0 - avfailL[v]) - art_failL;
  }

  xdr_overall = input[94];

  std::fill(art_mortE, art_mortE + STAGES, std::numeric_limits<double>::quiet_NaN());
  art_mortE[L2] = input[95] * input[97];
  art_mortE[L3] = input[95] * input[97];
  art_mortE[A1] = input[95];

  std::fill(art_mortL, art_mortL + STAGES, std::numeric_limits<double>::quiet_NaN());
  art_mortL[L2] = input[96] * input[97];
  art_mortL[L3] = input[96] * input[97];
  art_mortL[A1] = input[96];

  art_efficacy_wt = 1.0 - input[98];
  art_efficacy_dr = 1.0 - input[99];

  sec_init_dr = input[100];
  sec_init_na = input[101];
  for (int v(0); v < VARIANTS; ++v) {
    emerge[ART2][v] = -input[104] * log(1.0 - input[102]);
  }
  emerge[ART2][ S1] = emerge[ART2][WS1] = input[105] * emerge[ART2][WT];

  sec_fail_dr = -log(1.0 - input[102]) - emerge[ART2][WT];
  sec_fail_na = -log(1.0 - input[103]) - emerge[ART2][WT];

  std::fill(sec_mort, sec_mort + STAGES, std::numeric_limits<double>::quiet_NaN());
  sec_mort[L2] = input[106] * input[97];
  sec_mort[L3] = input[106] * input[97];
  sec_mort[A1] = input[106];

  sec_exit = input[107];

  // +-+ Drug resistance fitness parameters +----------------------------------+
  revert[NONE][C1] = 1.0 / input[108];
  revert[NONE][C2] = 1.0 / input[109];
  revert[NONE][R1] = 1.0 / input[110];
  revert[NONE][R2] = 1.0 / input[111];
  revert[NONE][S1] = 1.0 / input[112];
  revert[NONE][S2] = 1.0 / input[113];
  revert[NONE][Q1] = 1.0 / input[114];
  revert[NONE][Q2] = 1.0 / input[115];

  revert[PREP_HIGH][R1] = revert[PREP_POOR][R1] = revert[NONE][R1];
  revert[PREP_HIGH][S1] = revert[PREP_POOR][S1] = revert[NONE][S1];

  revert[ART2][R1] = revert[NONE][R1];
  revert[ART2][C1] = revert[NONE][C1];
  revert[ART2][Q1] = revert[NONE][Q1];
  revert[ART2][Q2] = revert[NONE][Q2];

  revert[ART2][R2] = revert[NONE][R2];
  revert[ART2][C2] = revert[NONE][C2];

  memcpy(revert[ART1_NA], revert[NONE], VARIANTS * sizeof(double));
  memcpy(revert[ART2_NA], revert[NONE], VARIANTS * sizeof(double));

  for (int h(0); h < STAGES; ++h) {
    prog[C1 ][h] = prog[WT][h] * input[116];
    prog[C2 ][h] = prog[WT][h] * input[117];
    prog[R1 ][h] = prog[WT][h] * input[118];
    prog[R2 ][h] = prog[WT][h] * input[119];
    prog[S1 ][h] = prog[WT][h] * input[120];
    prog[S2 ][h] = prog[WT][h] * input[121];
    prog[Q1 ][h] = prog[WT][h] * input[122];
    prog[Q2 ][h] = prog[WT][h] * input[123];
    prog[WQ1][h] = prog[WQ2][h] = prog[WT][h];
    prog[WC1][h] = prog[WC2][h] = prog[WT][h];
    prog[WR1][h] = prog[WR2][h] = prog[WT][h];
    prog[WS1][h] = prog[WS2][h] = prog[WT][h];
  }

  std::fill(infect_virus, infect_virus + VARIANTS, 1.0);
  infect_virus[C1] = input[124];
  infect_virus[C2] = input[125];
  infect_virus[R1] = input[126];
  infect_virus[R2] = input[127];
  infect_virus[S1] = input[128];
  infect_virus[S2] = input[129];
  infect_virus[Q1] = input[130];
  infect_virus[Q2] = input[131];

  set_infectivity();

  mmc_prop_late = input[132];
  mmc_time_scaleup = input[133];
  mmc_span_scaleup = input[134] - mmc_time_scaleup;

  return 0;
}

void Parameters::set_infectivity() {
  Drug u;
  Virus v;
  double base, prep;

  for (int s(Y1NWT); s < STATES; ++s) {
    base = infect_base * infect_stage[attr_stage[s]];
    u = attr_drug[s];
    v = attr_virus[s];
    if ((u == ART2) && (v != S2)) {
      // suppressed on second-line ART
      if (v == R1 || v == C1 || v == Q1 || v == S1 || v == R2 || v == C2 || v == Q2) {
	base *= art_efficacy_dr;
      } else {
	base *= art_efficacy_wt;
      }
    } else if ((u == ART1_NEW || u == ART1_OLD) && (v != R2) && (v != C2)) {
      // suppressed on first-line ART
      assert(v != S2); // added assertion; absurd to have 2nd-line ADR on 1st-line ART
      if (v == R1 || v == C1 || v == S1 || v == Q1 || v == Q2) {
	base *= art_efficacy_dr;
      } else {
	base *= art_efficacy_wt;
      }
    } else {
      base *= infect_virus[v];
    }

    // adjust PrEP efficacy for the donor's HIV variant
    prep = 1.0 - (1.0 - prep_efficacy) * prep_eff_virus[attr_transmit[s]];

    // Factor in other preventative interventions and log-transform
    // terms, since it is cheaper to calculate 1-exp(n*sum(log(1-a)))
    // than to calculate 1-product((1-a)^n) after caching the
    // log-transformed terms.
    infect[0][0][0][s] = log(1.0 - base);
    infect[0][0][1][s] = log(1.0 - base * prep);
    infect[0][1][0][s] = log(1.0 - base * mmc_efficacy);
    infect[0][1][1][s] = log(1.0 - base * mmc_efficacy * prep);
    infect[1][0][0][s] = log(1.0 - base * condom_efficacy);
    infect[1][0][1][s] = log(1.0 - base * condom_efficacy * prep);
    infect[1][1][0][s] = log(1.0 - base * condom_efficacy * mmc_efficacy);
    infect[1][1][1][s] = log(1.0 - base * condom_efficacy * mmc_efficacy * prep);
  }
}

int model(double t, double const y[], double dy[], Parameters const& p) {
  int b; // age band

  // discount factor used to accumulate outcomes for
  // cost-effectiveness analysis
  const double discount((t >= p.ce_tinit) * exp(-p.ce_discount * (t - p.ce_tinit)));

  const double prep_propP(1.0 - p.prep_prop_high), prep_propH(p.prep_prop_high);

  double const* z;
  double const* w;
  double *d;
  double debut[SEXES][LEVELS];
  double emergeR2S, emergeR2R, emergeC2S, emergeC2R;
  double emergeS2;
  double rate_age[AGES];
  double mort;

#warning "Second-line ART scale-up time points hard-coded"
  // second-line access
  double sec_access, sec_init_dr, sec_init_na;
  const double as1(0.000), ts1(p.art_time_start[A1]);
  const double as2(0.140), ts2(2015);
  const double as3(1.000), ts3(2020);
  if (t <= ts1) {
    sec_access = as1;
  } else if (t <= ts2) {
    sec_access = as1 + (as2 - as1) * (t - ts1) / (ts2 - ts1);
  } else if (t <= ts3) {
    sec_access = as2 + (as3 - as2) * (t - ts2) / (ts3 - ts2);
  } else {
    sec_access = as3;
  }
  sec_init_dr = sec_access * p.sec_init_dr;
  sec_init_na = sec_access * p.sec_init_na;

  std::fill(rate_age, rate_age + AGES, p.rate_age);
  rate_age[AGES-1] = 0.0; // Oldest age compartments catch all people aged 55+

  // maps sex and circumcision status to sex (F->F, MU->M, MC->M)
  const Sex sex[] = {F, M, M};

  // lambda[F ]: force of infection acting on women
  // lambda[MU]: force of infection acting on uncircumcised men
  // lambda[MC]: force of infection acting on circumcised men
  double lambda[SEXCIRC][BANDS][LEVELS][PREP_STATES][TRANSMIT];
  force(t, y, p, lambda);

  // MMC initiation parameters
  double mmc_prop, mmc_diff; // mmc_diff is required by mmc_cover but unused otherwise
  mmc_cover(t, p, mmc_prop, mmc_diff);
  double cuptake[BANDS];

  // PrEP initiation parameters
  double puptake[SEXES][BANDS][LEVELS];

  double prep_inject[SEXES][BANDS][LEVELS];
  double prep_exit[SEXES][BANDS][LEVELS];

  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS; ++b) {
	if (p.prep_init[g][b][k] < 0.0 || t >= p.prep_time_closure) {
	  prep_inject[g][b][k] = 0.0;
	  prep_exit[g][b][k] = p.prep_base_inject + p.prep_base_exit;
	} else {
	  prep_inject[g][b][k] = p.prep_base_inject;
	  prep_exit[g][b][k] = p.prep_base_exit;
	}
      }
    }
  }

  // ART initiation parameters
  double auptake[STAGES];
  rates_art(t, p, auptake);

  // population growth rate
  const double grow_param[] = {2 * p.rate_grow_init - p.rate_grow_late, p.rate_grow_late, p.rate_grow_decr, 1978};
  const double rate_grow(grow_param[0] + (grow_param[1] - grow_param[0]) / (1 + exp(-grow_param[2] * (t - grow_param[3]))));

  double n(0.0);
  for (int u(0); u < TRACKS; ++u) {
    for (int g(0); g < SEXCIRC; ++g) {
      for (int a(0); a < AGES-1; ++a) {
	for (int k(0); k < LEVELS; ++k) {
	  z = offset(g, k, a, u) + y;
	  n += std::accumulate(z, z + STATES, 0.0);
	}
      }
    }
  }

  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      debut[g][k] = rate_grow * p.prop_risk[g][k] * n;
    }
  }

  std::fill(dy, dy + COMPARTMENTS, 0.0);

  // 1. handle everything except demographic-related inflows 
  for (int u(0); u < TRACKS; ++u) {
    for (int a(0); a < AGES; ++a) {
      b = band[a];
      for (int k(0); k < LEVELS; ++k) {

	// established infection and later (compartments where
	// circumcision status is irrelevant)
	for (int g(0); g < SEXCIRC; ++g) {
	  z = offset(g, k, a, u) + y;
	  d = offset(g, k, a, u) + dy;

	  // +-+ acute window-period infection +---------------------------------+
	  // Note that PrEP concentration dynamics are ignored for
	  // individuals with PrEP ADR (Q2) here. These dynamics are
	  // irrelevant among infected individuals with PrEP ADR,
	  // since high or poor concentrations are assumed sufficient
	  // to maintain resistance. Concentration dynamics are
	  // present in later stages to more explicitly account for
	  // outcomes of follow-up injections (removal via HIV
	  // testing, removal via ineligibility and concentration
	  // change).
	  d[Y1NQ2] = prep_exit[sex[g]][b][k] * (z[Y1HQ2] + z[Y1PQ2]) - (p.prog[Q2][P1] + p.revert[NONE][Q2]) * z[Y1NQ2];
	  d[Y1HQ2] = p.emerge[PREP_HIGH][WT] * z[Y1HWT] + p.emerge[PREP_HIGH][R1] * z[Y1HR1] + p.emerge[PREP_HIGH][WR1] * z[Y1Hr1] + p.emerge[PREP_HIGH][S1] * z[Y1HS1] + p.emerge[PREP_HIGH][WS1] * z[Y1Hs1] + p.emerge[PREP_HIGH][WC1] * z[Y1Hc1] + p.emerge[PREP_HIGH][WQ1] * z[Y1Hq1] - (p.prog[Q2][P1] + prep_exit[sex[g]][b][k]) * z[Y1HQ2];
	  d[Y1PQ2] = p.emerge[PREP_POOR][WT] * z[Y1PWT] + p.emerge[PREP_POOR][R1] * z[Y1PR1] + p.emerge[PREP_POOR][WR1] * z[Y1Pr1] + p.emerge[PREP_POOR][S1] * z[Y1PS1] + p.emerge[PREP_POOR][WS1] * z[Y1Ps1] + p.emerge[PREP_POOR][WC1] * z[Y1Pc1] + p.emerge[PREP_POOR][WQ1] * z[Y1Pq1] - (p.prog[Q2][P1] + prep_exit[sex[g]][b][k]) * z[Y1PQ2];

	  d[Y1Nr1] = p.revert[NONE][R1] * z[Y1NR1] + prep_exit[sex[g]][b][k] * (z[Y1Hr1] + z[Y1Pr1]) - (p.prog[WR1][P1]) * z[Y1Nr1];
	  d[Y1Hr1] = p.revert[PREP_HIGH][R1] * z[Y1HR1] + prep_propH * prep_inject[sex[g]][b][k] * z[Y1Pr1] - (p.prog[WR1][P1] + p.emerge[PREP_HIGH][WR1] + prep_inject[sex[g]][b][k] * prep_propP + prep_exit[sex[g]][b][k]) * z[Y1Hr1];
	  d[Y1Pr1] = p.revert[PREP_POOR][R1] * z[Y1PR1] + prep_propP * prep_inject[sex[g]][b][k] * z[Y1Hr1] - (p.prog[WR1][P1] + p.emerge[PREP_POOR][WR1] + prep_inject[sex[g]][b][k] * prep_propH + prep_exit[sex[g]][b][k]) * z[Y1Pr1];

	  d[Y1Nc1] = p.revert[NONE][C1] * z[Y1NC1] + prep_exit[sex[g]][b][k] * (z[Y1Hc1] + z[Y1Pc1]) - (p.prog[WC1][P1]) * z[Y1Nc1];
	  d[Y1Hc1] = prep_propH * (prep_inject[sex[g]][b][k] * z[Y1Pc1]) - (p.prog[WC1][P1] + p.emerge[PREP_HIGH][WC1] + prep_inject[sex[g]][b][k] * prep_propP + prep_exit[sex[g]][b][k]) * z[Y1Hc1];
	  d[Y1Pc1] = prep_propP * (prep_inject[sex[g]][b][k] * z[Y1Hc1]) - (p.prog[WC1][P1] + p.emerge[PREP_POOR][WC1] + prep_inject[sex[g]][b][k] * prep_propH + prep_exit[sex[g]][b][k]) * z[Y1Pc1];

	  d[Y1Nq1] = p.revert[NONE][Q1] * z[Y1NQ1] + prep_exit[sex[g]][b][k] * (z[Y1Hq1] + z[Y1Pq1]) - (p.prog[WQ1][P1]) * z[Y1Nq1];
	  d[Y1Hq1] = prep_propH * (prep_inject[sex[g]][b][k] * z[Y1Pq1]) - (p.prog[WQ1][P1] + p.emerge[PREP_HIGH][WQ1] + prep_inject[sex[g]][b][k] * prep_propP + prep_exit[sex[g]][b][k]) * z[Y1Hq1];
	  d[Y1Pq1] = prep_propP * (prep_inject[sex[g]][b][k] * z[Y1Hq1]) - (p.prog[WQ1][P1] + p.emerge[PREP_POOR][WQ1] + prep_inject[sex[g]][b][k] * prep_propH + prep_exit[sex[g]][b][k]) * z[Y1Pq1];

 	  d[Y1Ns1] = p.revert[NONE][S1] * z[Y1NS1] + prep_exit[sex[g]][b][k] * (z[Y1Hs1] + z[Y1Ps1]) - (p.prog[WS1][P1]) * z[Y1Ns1];
 	  d[Y1Hs1] = p.revert[PREP_HIGH][S1] * z[Y1HS1] + prep_propH * prep_inject[sex[g]][b][k] * z[Y1Ps1] - (p.prog[WS1][P1] + p.emerge[PREP_HIGH][WS1] + prep_inject[sex[g]][b][k] * prep_propP + prep_exit[sex[g]][b][k]) * z[Y1Hs1];
 	  d[Y1Ps1] = p.revert[PREP_POOR][S1] * z[Y1PS1] + prep_propP * prep_inject[sex[g]][b][k] * z[Y1Hs1] - (p.prog[WS1][P1] + p.emerge[PREP_POOR][WS1] + prep_inject[sex[g]][b][k] * prep_propH + prep_exit[sex[g]][b][k]) * z[Y1Ps1];

	  d[Y1Nq2] = p.revert[NONE][Q2] * z[Y1NQ2] - (p.prog[WQ2][P1]) * z[Y1Nq2];

	  // +-+ acute established infection +-----------------------------------+
	  d[Y2NWT] = p.prog[WT][P1] * z[Y1NWT] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y2HWT] + z[Y2PWT]) - (p.prog[WT][P2]) * z[Y2NWT];
	  d[Y2HWT] = p.prog[WT][P1] * z[Y1HWT] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2PWT] - (p.prog[WT][P2] + p.emerge[PREP_HIGH][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y2HWT];
	  d[Y2PWT] = p.prog[WT][P1] * z[Y1PWT] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2HWT] - (p.prog[WT][P2] + p.emerge[PREP_POOR][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y2PWT];

	  d[Y2NR1] = p.prog[R1][P1] * z[Y1NR1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y2HR1] + z[Y2PR1]) - (p.prog[R1][P2] + p.revert[NONE][R1]) * z[Y2NR1];
	  d[Y2HR1] = p.prog[R1][P1] * z[Y1HR1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2PR1] - (p.prog[R1][P2] + p.revert[PREP_HIGH][R1] + p.emerge[PREP_HIGH][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y2HR1];
	  d[Y2PR1] = p.prog[R1][P1] * z[Y1PR1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2HR1] - (p.prog[R1][P2] + p.revert[PREP_POOR][R1] + p.emerge[PREP_POOR][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y2PR1];

	  d[Y2NC1] = p.prog[C1][P1] * z[Y1NC1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y2HC1] + z[Y2PC1]) - (p.prog[C1][P2] + p.revert[NONE][C1]) * z[Y2NC1];
	  d[Y2HC1] = p.prog[C1][P1] * z[Y1HC1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2PC1] - (p.prog[C1][P2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y2HC1];
	  d[Y2PC1] = p.prog[C1][P1] * z[Y1PC1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2HC1] - (p.prog[C1][P2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y2PC1];

	  d[Y2NQ1] = p.prog[Q1][P1] * z[Y1NQ1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y2HQ1] + z[Y2PQ1]) - (p.prog[Q1][P2] + p.revert[NONE][Q1]) * z[Y2NQ1];
	  d[Y2HQ1] = p.prog[Q1][P1] * z[Y1HQ1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2PQ1] - (p.prog[Q1][P2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y2HQ1];
	  d[Y2PQ1] = p.prog[Q1][P1] * z[Y1PQ1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2HQ1] - (p.prog[Q1][P2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y2PQ1];

 	  d[Y2NS1] = p.prog[S1][P1] * z[Y1NS1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y2HS1] + z[Y2PS1]) - (p.prog[S1][P2] + p.revert[NONE][S1]) * z[Y2NS1];
 	  d[Y2HS1] = p.prog[S1][P1] * z[Y1HS1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2PS1] - (p.prog[S1][P2] + p.revert[PREP_HIGH][S1] + p.emerge[PREP_HIGH][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y2HS1];
 	  d[Y2PS1] = p.prog[S1][P1] * z[Y1PS1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2HS1] - (p.prog[S1][P2] + p.revert[PREP_POOR][S1] + p.emerge[PREP_POOR][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y2PS1];

	  d[Y2NQ2] = p.prog[Q2][P1] * z[Y1NQ2] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y2HQ2] + z[Y2PQ2]) - (p.prog[Q2][P2] + p.revert[NONE][Q2]) * z[Y2NQ2];
	  d[Y2HQ2] = p.prog[Q2][P1] * z[Y1HQ2] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2PQ2] + p.emerge[PREP_HIGH][WT] * z[Y2HWT] + p.emerge[PREP_HIGH][R1] * z[Y2HR1] + p.emerge[PREP_HIGH][WR1] * z[Y2Hr1] + p.emerge[PREP_HIGH][S1] * z[Y2HS1] + p.emerge[PREP_HIGH][WS1] * z[Y2Hs1] + p.emerge[PREP_HIGH][WC1] * z[Y2Hc1] + p.emerge[PREP_HIGH][WQ1] * z[Y2Hq1] - (p.prog[Q2][P2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y2HQ2];
	  d[Y2PQ2] = p.prog[Q2][P1] * z[Y1PQ2] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2HQ2] + p.emerge[PREP_POOR][WT] * z[Y2PWT] + p.emerge[PREP_POOR][R1] * z[Y2PR1] + p.emerge[PREP_POOR][WR1] * z[Y2Pr1] + p.emerge[PREP_POOR][S1] * z[Y2PS1] + p.emerge[PREP_POOR][WS1] * z[Y2Ps1] + p.emerge[PREP_POOR][WC1] * z[Y2Pc1] + p.emerge[PREP_POOR][WQ1] * z[Y2Pq1] - (p.prog[Q2][P2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y2PQ2];

	  d[Y2Nr1] = p.prog[WR1][P1] * z[Y1Nr1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y2Hr1] + z[Y2Pr1]) + p.revert[NONE][R1] * z[Y2NR1] - (p.prog[WR1][P2]) * z[Y2Nr1];
	  d[Y2Hr1] = p.prog[WR1][P1] * z[Y1Hr1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2Pr1] + p.revert[PREP_HIGH][R1] * z[Y2HR1] - (p.prog[WR1][P2] + p.emerge[PREP_HIGH][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y2Hr1];
	  d[Y2Pr1] = p.prog[WR1][P1] * z[Y1Pr1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2Hr1] + p.revert[PREP_POOR][R1] * z[Y2PR1] - (p.prog[WR1][P2] + p.emerge[PREP_POOR][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y2Pr1];

	  d[Y2Nc1] = p.prog[WC1][P1] * z[Y1Nc1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y2Hc1] + z[Y2Pc1]) + p.revert[NONE][C1] * z[Y2NC1] - (p.prog[WC1][P2]) * z[Y2Nc1];
	  d[Y2Hc1] = p.prog[WC1][P1] * z[Y1Hc1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2Pc1] - (p.prog[WC1][P2] + p.emerge[PREP_HIGH][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y2Hc1];
	  d[Y2Pc1] = p.prog[WC1][P1] * z[Y1Pc1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2Hc1] - (p.prog[WC1][P2] + p.emerge[PREP_POOR][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y2Pc1];

	  d[Y2Nq1] = p.prog[WQ1][P1] * z[Y1Nq1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y2Hq1] + z[Y2Pq1]) + p.revert[NONE][Q1] * z[Y2NQ1] - (p.prog[WQ1][P2]) * z[Y2Nq1];
	  d[Y2Hq1] = p.prog[WQ1][P1] * z[Y1Hq1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2Pq1] - (p.prog[WQ1][P2] + p.emerge[PREP_HIGH][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y2Hq1];
	  d[Y2Pq1] = p.prog[WQ1][P1] * z[Y1Pq1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2Hq1] - (p.prog[WQ1][P2] + p.emerge[PREP_POOR][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y2Pq1];

	  d[Y2Ns1] = p.prog[WS1][P1] * z[Y1Ns1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y2Hs1] + z[Y2Ps1]) + p.revert[NONE][S1] * z[Y2NS1] - (p.prog[WS1][P2]) * z[Y2Ns1];
	  d[Y2Hs1] = p.prog[WS1][P1] * z[Y1Hs1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2Ps1] + p.revert[PREP_HIGH][S1] * z[Y2HS1] - (p.prog[WS1][P2] + p.emerge[PREP_HIGH][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y2Hs1];
	  d[Y2Ps1] = p.prog[WS1][P1] * z[Y1Ps1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y2Hs1] + p.revert[PREP_POOR][S1] * z[Y2PS1] - (p.prog[WS1][P2] + p.emerge[PREP_POOR][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y2Ps1];

	  d[Y2Nq2] = p.prog[WQ2][P1] * z[Y1Nq2] + p.revert[NONE][Q2] * z[Y2NQ2] - (p.prog[WQ2][P2]) * z[Y2Nq2];

#ifdef CD4_500
	  // If CD4_500 is defined, chronic infection is stratified into
	  // three stages: CD4>500, 350<CD4<500 and 200<CD4<350. This
	  // block of code implements the stratification between CD4>500
	  // and 350<CD4<500. 200<CD4<350 is included in the model
	  // whether or not CD4_500 is defined.

	  // +-+ chronic infection (CD4>500) +-----------------------------------+
	  d[Y3NWT] = p.prog[WT][P2] * z[Y2NWT] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y3HWT] + z[Y3PWT]) - (p.prog[WT][L1]) * z[Y3NWT];
	  d[Y3HWT] = p.prog[WT][P2] * z[Y2HWT] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3PWT] - (p.prog[WT][L1] + p.emerge[PREP_HIGH][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y3HWT];
	  d[Y3PWT] = p.prog[WT][P2] * z[Y2PWT] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3HWT] - (p.prog[WT][L1] + p.emerge[PREP_POOR][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y3PWT];

	  d[Y3NR1] = p.prog[R1][P2] * z[Y2NR1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y3HR1] + z[Y3PR1]) - (p.prog[R1][L1] + p.revert[NONE][R1]) * z[Y3NR1];
	  d[Y3HR1] = p.prog[R1][P2] * z[Y2HR1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3PR1] - (p.prog[R1][L1] + p.revert[PREP_HIGH][R1] + p.emerge[PREP_HIGH][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y3HR1];
	  d[Y3PR1] = p.prog[R1][P2] * z[Y2PR1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3HR1] - (p.prog[R1][L1] + p.revert[PREP_POOR][R1] + p.emerge[PREP_POOR][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y3PR1];

	  d[Y3NC1] = p.prog[C1][P2] * z[Y2NC1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y3HC1] + z[Y3PC1]) - (p.prog[C1][L1] + p.revert[NONE][C1]) * z[Y3NC1];
	  d[Y3HC1] = p.prog[C1][P2] * z[Y2HC1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3PC1] - (p.prog[C1][L1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y3HC1];
	  d[Y3PC1] = p.prog[C1][P2] * z[Y2PC1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3HC1] - (p.prog[C1][L1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y3PC1];

	  d[Y3NQ1] = p.prog[Q1][P2] * z[Y2NQ1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y3HQ1] + z[Y3PQ1]) - (p.prog[Q1][L1] + p.revert[NONE][Q1]) * z[Y3NQ1];
	  d[Y3HQ1] = p.prog[Q1][P2] * z[Y2HQ1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3PQ1] - (p.prog[Q1][L1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y3HQ1];
	  d[Y3PQ1] = p.prog[Q1][P2] * z[Y2PQ1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3HQ1] - (p.prog[Q1][L1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y3PQ1];

	  d[Y3NS1] = p.prog[S1][P2] * z[Y2NS1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y3HS1] + z[Y3PS1]) - (p.prog[S1][L1] + p.revert[NONE][S1]) * z[Y3NS1];
	  d[Y3HS1] = p.prog[S1][P2] * z[Y2HS1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3PS1] - (p.prog[S1][L1] + p.revert[PREP_HIGH][S1] + p.emerge[PREP_HIGH][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y3HS1];
	  d[Y3PS1] = p.prog[S1][P2] * z[Y2PS1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3HS1] - (p.prog[S1][L1] + p.revert[PREP_POOR][S1] + p.emerge[PREP_POOR][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y3PS1];

	  d[Y3NQ2] = p.prog[Q2][P2] * z[Y2NQ2] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y3HQ2] + z[Y3PQ2]) - (p.prog[Q2][L1] + p.revert[NONE][Q2]) * z[Y3NQ2];
	  d[Y3HQ2] = p.prog[Q2][P2] * z[Y2HQ2] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3PQ2] + p.emerge[PREP_HIGH][WT] * z[Y3HWT] + p.emerge[PREP_HIGH][R1] * z[Y3HR1] + p.emerge[PREP_HIGH][WR1] * z[Y3Hr1] + p.emerge[PREP_HIGH][S1] * z[Y3HS1] + p.emerge[PREP_HIGH][WS1] * z[Y3Hs1] + p.emerge[PREP_HIGH][WC1] * z[Y3Hc1] + p.emerge[PREP_HIGH][WQ1] * z[Y3Hq1] - (p.prog[Q2][L1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y3HQ2];
	  d[Y3PQ2] = p.prog[Q2][P2] * z[Y2PQ2] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3HQ2] + p.emerge[PREP_POOR][WT] * z[Y3PWT] + p.emerge[PREP_POOR][R1] * z[Y3PR1] + p.emerge[PREP_POOR][WR1] * z[Y3Pr1] + p.emerge[PREP_POOR][S1] * z[Y3PS1] + p.emerge[PREP_POOR][WS1] * z[Y3Ps1] + p.emerge[PREP_POOR][WC1] * z[Y3Pc1] + p.emerge[PREP_POOR][WQ1] * z[Y3Pq1] - (p.prog[Q2][L1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y3PQ2];

	  d[Y3Nr1] = p.prog[WR1][P2] * z[Y2Nr1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y3Hr1] + z[Y3Pr1]) + p.revert[NONE][R1] * z[Y3NR1] - (p.prog[WR1][L1]) * z[Y3Nr1];
	  d[Y3Hr1] = p.prog[WR1][P2] * z[Y2Hr1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3Pr1] + p.revert[PREP_HIGH][R1] * z[Y3HR1] - (p.prog[WR1][L1] + p.emerge[PREP_HIGH][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y3Hr1];
	  d[Y3Pr1] = p.prog[WR1][P2] * z[Y2Pr1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3Hr1] + p.revert[PREP_POOR][R1] * z[Y3PR1] - (p.prog[WR1][L1] + p.emerge[PREP_POOR][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y3Pr1];

	  d[Y3Nc1] = p.prog[WC1][P2] * z[Y2Nc1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y3Hc1] + z[Y3Pc1]) + p.revert[NONE][C1] * z[Y3NC1] - (p.prog[WC1][L1]) * z[Y3Nc1];
	  d[Y3Hc1] = p.prog[WC1][P2] * z[Y2Hc1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3Pc1] - (p.prog[WC1][L1] + p.emerge[PREP_HIGH][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y3Hc1];
	  d[Y3Pc1] = p.prog[WC1][P2] * z[Y2Pc1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3Hc1] - (p.prog[WC1][L1] + p.emerge[PREP_POOR][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y3Pc1];

	  d[Y3Nq1] = p.prog[WQ1][P2] * z[Y2Nq1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y3Hq1] + z[Y3Pq1]) + p.revert[NONE][Q1] * z[Y3NQ1] - (p.prog[WQ1][L1]) * z[Y3Nq1];
	  d[Y3Hq1] = p.prog[WQ1][P2] * z[Y2Hq1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3Pq1] - (p.prog[WQ1][L1] + p.emerge[PREP_HIGH][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y3Hq1];
	  d[Y3Pq1] = p.prog[WQ1][P2] * z[Y2Pq1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3Hq1] - (p.prog[WQ1][L1] + p.emerge[PREP_POOR][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y3Pq1];

 	  d[Y3Ns1] = p.prog[WS1][P2] * z[Y2Ns1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y3Hs1] + z[Y3Ps1]) + p.revert[NONE][S1] * z[Y3NS1] - (p.prog[WS1][L1]) * z[Y3Ns1];
 	  d[Y3Hs1] = p.prog[WS1][P2] * z[Y2Hs1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3Ps1] + p.revert[PREP_HIGH][S1] * z[Y3HS1] - (p.prog[WS1][L1] + p.emerge[PREP_HIGH][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y3Hs1];
 	  d[Y3Ps1] = p.prog[WS1][P2] * z[Y2Ps1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y3Hs1] + p.revert[PREP_POOR][S1] * z[Y3PS1] - (p.prog[WS1][L1] + p.emerge[PREP_POOR][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y3Ps1];

	  d[Y3Nq2] = p.prog[WQ2][P2] * z[Y2Nq2] + p.revert[NONE][Q2] * z[Y3NQ2] - (p.prog[WQ2][L1]) * z[Y3Nq2];

	  // +-+ chronic infection (350<CD4<500) +-------------------------------+
	  d[Y4NWT] = p.prog[WT][L1] * z[Y3NWT] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HWT] + z[Y4PWT]) + p.art_exitE * z[Y4EWT] + p.art_exitL * z[Y4LWT] + p.art_exitF * z[Y4FWT] - (p.prog[WT][L2] + auptake[L2]) * z[Y4NWT];
	  d[Y4HWT] = p.prog[WT][L1] * z[Y3HWT] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PWT] - (p.prog[WT][L2] + p.emerge[PREP_HIGH][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HWT];
	  d[Y4PWT] = p.prog[WT][L1] * z[Y3PWT] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HWT] - (p.prog[WT][L2] + p.emerge[PREP_POOR][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PWT];
	  d[Y4EWT] = auptake[L2] * z[Y4NWT] - (p.emerge[ART1_NEW][WT] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4EWT];
	  d[Y4LWT] = p.art_advance * z[Y4EWT] - (p.emerge[ART1_OLD][WT] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4LWT];
	  d[Y4FWT] = p.art_failE * z[Y4EWT] + p.art_failL * z[Y4LWT] - (p.prog[WT][L2] + p.art_exitF + sec_init_na) * z[Y4FWT];
	  d[Y4TWT] = sec_init_na * z[Y4FWT] + auptake[L2] * z[Y4DWT] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WT]) * z[Y4TWT]; // 2014-10-20 DONE
	  d[Y4DWT] = p.sec_exit * (z[Y4TWT] + z[Y4UWT]) - (auptake[L2] + p.prog[WT][L2]) * z[Y4DWT]; // 2014-10-20 DONE
	  d[Y4UWT] = p.sec_fail_na * z[Y4TWT] - (p.prog[WT][L2] + p.sec_exit) * z[Y4UWT]; // 2014-10-20 DONE

	  d[Y4NR1] = p.prog[R1][L1] * z[Y3NR1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HR1] + z[Y4PR1]) + p.art_exitE * z[Y4ER1] + p.art_exitL * z[Y4LR1] + p.art_exitF * z[Y4FR1] - (p.prog[R1][L2] + p.revert[NONE][R1] + auptake[L2]) * z[Y4NR1];
	  d[Y4HR1] = p.prog[R1][L1] * z[Y3HR1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PR1] - (p.prog[R1][L2] + p.revert[PREP_HIGH][R1] + p.emerge[PREP_HIGH][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HR1];
	  d[Y4PR1] = p.prog[R1][L1] * z[Y3PR1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HR1] - (p.prog[R1][L2] + p.revert[PREP_POOR][R1] + p.emerge[PREP_POOR][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PR1];
	  d[Y4ER1] = auptake[L2] * z[Y4NR1] - (p.emerge[ART1_NEW][R1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4ER1];
	  d[Y4LR1] = p.art_advance * z[Y4ER1] - (p.emerge[ART1_OLD][R1] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4LR1];
	  d[Y4FR1] = p.art_failE * z[Y4ER1] + p.art_failL * z[Y4LR1] - (p.prog[R1][L2] + p.revert[ART1_NA][R1] + p.art_exitF + sec_init_na) * z[Y4FR1];
	  d[Y4TR1] = sec_init_na * z[Y4FR1] + auptake[L2] * z[Y4DR1] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.revert[ART2][R1] + p.emerge[ART2][R1]) * z[Y4TR1];
	  d[Y4DR1] = p.sec_exit * (z[Y4TR1] + z[Y4UR1]) - (auptake[L2] + p.prog[R1][L2] + p.revert[NONE][R1]) * z[Y4DR1]; // 2014-10-20 DONE
	  d[Y4UR1] = p.sec_fail_na * z[Y4TR1] - (p.prog[R1][L2] + p.revert[ART2_NA][R1] + p.sec_exit) * z[Y4UR1]; // 2014-10-20 DONE

	  d[Y4NC1] = p.prog[C1][L1] * z[Y3NC1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HC1] + z[Y4PC1]) + p.art_exitE * z[Y4EC1] + p.art_exitL * z[Y4LC1] + p.art_exitF * z[Y4FC1] - (p.prog[C1][L2] + p.revert[NONE][C1] + auptake[L2]) * z[Y4NC1];
	  d[Y4HC1] = p.prog[C1][L1] * z[Y3HC1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PC1] - (p.prog[C1][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HC1];
	  d[Y4PC1] = p.prog[C1][L1] * z[Y3PC1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HC1] - (p.prog[C1][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PC1];
	  d[Y4EC1] = auptake[L2] * z[Y4NC1] - (p.emerge[ART1_NEW][C1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4EC1];
	  d[Y4LC1] = p.art_advance * z[Y4EC1] - (p.emerge[ART1_OLD][C1] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4LC1];
	  d[Y4FC1] = p.art_failE * z[Y4EC1] + p.art_failL * z[Y4LC1] - (p.prog[C1][L2] + p.revert[ART1_NA][C1] + p.art_exitF + sec_init_na) * z[Y4FC1];
	  d[Y4TC1] = sec_init_na * z[Y4FC1] + auptake[L2] * z[Y4DC1] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.revert[ART2][C1] + p.emerge[ART2][C1]) * z[Y4TC1];
	  d[Y4DC1] = p.sec_exit * (z[Y4TC1] + z[Y4UC1]) - (auptake[L2] + p.prog[C1][L2] + p.revert[NONE][C1]) * z[Y4DC1]; // 2014-10-20 DONE
	  d[Y4UC1] = p.sec_fail_na * z[Y4TC1] - (p.prog[C1][L2] + p.revert[ART2_NA][C1] + p.sec_exit) * z[Y4UC1]; // 2014-10-20 DONE

	  d[Y4NQ1] = p.prog[Q1][L1] * z[Y3NQ1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HQ1] + z[Y4PQ1]) + p.art_exitE * z[Y4EQ1] + p.art_exitL * z[Y4LQ1] + p.art_exitF * z[Y4FQ1] - (p.prog[Q1][L2] + p.revert[NONE][Q1] + auptake[L2]) * z[Y4NQ1];
	  d[Y4HQ1] = p.prog[Q1][L1] * z[Y3HQ1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PQ1] - (p.prog[Q1][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HQ1];
	  d[Y4PQ1] = p.prog[Q1][L1] * z[Y3PQ1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HQ1] - (p.prog[Q1][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PQ1];
	  d[Y4EQ1] = auptake[L2] * z[Y4NQ1] - (p.emerge[ART1_NEW][Q1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4EQ1];
	  d[Y4LQ1] = p.art_advance * z[Y4EQ1] - (p.emerge[ART1_OLD][Q1] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4LQ1];
	  d[Y4FQ1] = p.art_failE * z[Y4EQ1] + p.art_failL * z[Y4LQ1] - (p.prog[Q1][L2] + p.revert[ART1_NA][Q1] + p.art_exitF + sec_init_na) * z[Y4FQ1];
	  d[Y4TQ1] = sec_init_na * z[Y4FQ1] + auptake[L2] * z[Y4DQ1] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.revert[ART2][Q1] + p.emerge[ART2][Q1]) * z[Y4TQ1];
	  d[Y4DQ1] = p.sec_exit * (z[Y4TQ1] + z[Y4UQ1]) - (auptake[L2] + p.prog[Q1][L2] + p.revert[NONE][Q1]) * z[Y4DQ1]; // 2014-10-20 DONE
	  d[Y4UQ1] = p.sec_fail_na * z[Y4TQ1] - (p.prog[Q1][L2] + p.revert[ART2_NA][Q1] + p.sec_exit) * z[Y4UQ1]; // 2014-10-20 DONE

 	  d[Y4NS1] = p.prog[S1][L1] * z[Y3NS1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HS1] + z[Y4PS1]) + p.art_exitE * z[Y4ES1] + p.art_exitL * z[Y4LS1] + p.art_exitF * z[Y4FS1] - (p.prog[S1][L2] + p.revert[NONE][S1] + auptake[L2]) * z[Y4NS1];
 	  d[Y4HS1] = p.prog[S1][L1] * z[Y3HS1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PS1] - (p.prog[S1][L2] + p.revert[PREP_HIGH][S1] + p.emerge[PREP_HIGH][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HS1];
 	  d[Y4PS1] = p.prog[S1][L1] * z[Y3PS1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HS1] - (p.prog[S1][L2] + p.revert[PREP_POOR][S1] + p.emerge[PREP_POOR][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PS1];
 	  d[Y4ES1] = auptake[L2] * z[Y4NS1] - (p.emerge[ART1_NEW][S1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4ES1];
 	  d[Y4LS1] = p.art_advance * z[Y4ES1] - (p.emerge[ART1_OLD][S1] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4LS1];
 	  d[Y4FS1] = p.art_failE * z[Y4ES1] + p.art_failL * z[Y4LS1] - (p.prog[S1][L2] + p.revert[ART1_NA][S1] + p.art_exitF + sec_init_na) * z[Y4FS1];
	  d[Y4TS1] = sec_init_na * z[Y4FS1] + auptake[L2] * z[Y4DS1] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][S1]) * z[Y4TS1]; // 2014-10-20 DONE
	  d[Y4DS1] = p.sec_exit * (z[Y4TS1] + z[Y4US1]) - (auptake[L2] + p.prog[S1][L2] + p.revert[NONE][S1]) * z[Y4DS1]; // 2014-10-20 DONE
	  d[Y4US1] = p.sec_fail_na * z[Y4TS1] - (p.prog[S1][L2] + p.revert[ART2_NA][S1] + p.sec_exit) * z[Y4US1]; // 2014-10-20 DONE

	  d[Y4NQ2] = p.prog[Q2][L1] * z[Y3NQ2] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HQ2] + z[Y4PQ2]) + p.art_exitE * z[Y4EQ2] + p.art_exitL * z[Y4LQ2] + p.art_exitF * z[Y4FQ2] - (p.prog[Q2][L2] + p.revert[NONE][Q2] + auptake[L2]) * z[Y4NQ2];
	  d[Y4HQ2] = p.prog[Q2][L1] * z[Y3HQ2] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PQ2] + p.emerge[PREP_HIGH][WT] * z[Y4HWT] + p.emerge[PREP_HIGH][R1] * z[Y4HR1] + p.emerge[PREP_HIGH][WR1] * z[Y4Hr1] + p.emerge[PREP_HIGH][S1] * z[Y4HS1] + p.emerge[PREP_HIGH][WS1] * z[Y4Hs1] + p.emerge[PREP_HIGH][WC1] * z[Y4Hc1] + p.emerge[PREP_HIGH][WQ1] * z[Y4Hq1] - (p.prog[Q2][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HQ2];
	  d[Y4PQ2] = p.prog[Q2][L1] * z[Y3PQ2] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HQ2] + p.emerge[PREP_POOR][WT] * z[Y4PWT] + p.emerge[PREP_POOR][R1] * z[Y4PR1] + p.emerge[PREP_POOR][WR1] * z[Y4Pr1] + p.emerge[PREP_POOR][S1] * z[Y4PS1] + p.emerge[PREP_POOR][WS1] * z[Y4Ps1] + p.emerge[PREP_POOR][WC1] * z[Y4Pc1] + p.emerge[PREP_POOR][WQ1] * z[Y4Pq1] - (p.prog[Q2][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PQ2];
	  d[Y4EQ2] = auptake[L2] * z[Y4NQ2] - (p.emerge[ART1_NEW][Q2] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4EQ2];
	  d[Y4LQ2] = p.art_advance * z[Y4EQ2] - (p.emerge[ART1_OLD][Q2] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4LQ2];
	  d[Y4FQ2] = p.art_failE * z[Y4EQ2] + p.art_failL * z[Y4LQ2] - (p.prog[Q2][L2] + p.revert[ART1_NA][Q2] + p.art_exitF + sec_init_na) * z[Y4FQ2];
	  d[Y4TQ2] = sec_init_na * z[Y4FQ2] + auptake[L2] * z[Y4DQ2] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.revert[ART2][Q2] + p.emerge[ART2][Q2]) * z[Y4TQ2];
	  d[Y4DQ2] = p.sec_exit * (z[Y4TQ2] + z[Y4UQ2]) - (auptake[L2] + p.prog[Q2][L2] + p.revert[NONE][Q2]) * z[Y4DQ2]; // 2014-10-20 DONE
	  d[Y4UQ2] = p.sec_fail_na * z[Y4TQ2] - (p.prog[Q2][L2] + p.revert[ART2_NA][Q2] + p.sec_exit) * z[Y4UQ2]; // 2014-10-20 DONE

	  emergeR2S = (1.0 - p.xdr_overall) *
	    (p.emerge[ART1_NEW][WT ] * z[Y4EWT] + p.emerge[ART1_OLD][WT ] * z[Y4LWT] +
	     p.emerge[ART1_NEW][R1 ] * z[Y4ER1] + p.emerge[ART1_OLD][R1 ] * z[Y4LR1] +
	     p.emerge[ART1_NEW][WR1] * z[Y4Er1] + p.emerge[ART1_OLD][WR1] * z[Y4Lr1]);
	  emergeR2R = (1.0 - p.xdr_overall) *
	    (p.emerge[ART1_NEW][S1 ] * z[Y4ES1] + p.emerge[ART1_OLD][S1 ] * z[Y4LS1] +
	     p.emerge[ART1_NEW][WS1] * z[Y4Es1] + p.emerge[ART1_OLD][WS1] * z[Y4Ls1]);

	  emergeC2S
	    = p.xdr_overall
	    * (p.emerge[ART1_NEW][WT ] * z[Y4EWT] + p.emerge[ART1_OLD][WT ] * z[Y4LWT] +
	       p.emerge[ART1_NEW][R1 ] * z[Y4ER1] + p.emerge[ART1_OLD][R1 ] * z[Y4LR1] +
	       p.emerge[ART1_NEW][WR1] * z[Y4Er1] + p.emerge[ART1_OLD][WR1] * z[Y4Lr1])
	    + (p.emerge[ART1_NEW][C1 ] * z[Y4EC1] + p.emerge[ART1_OLD][C1 ] * z[Y4LC1] +
	       p.emerge[ART1_NEW][WC1] * z[Y4Ec1] + p.emerge[ART1_OLD][WC1] * z[Y4Lc1] +
	       p.emerge[ART1_NEW][Q1 ] * z[Y4EQ1] + p.emerge[ART1_OLD][Q1 ] * z[Y4LQ1] +
	       p.emerge[ART1_NEW][WQ1] * z[Y4Eq1] + p.emerge[ART1_OLD][WQ1] * z[Y4Lq1] +
	       p.emerge[ART1_NEW][Q2 ] * z[Y4EQ2] + p.emerge[ART1_OLD][Q2 ] * z[Y4LQ2] +
	       p.emerge[ART1_NEW][WQ2] * z[Y4Eq2] + p.emerge[ART1_OLD][WQ2] * z[Y4Lq2]);
	  emergeC2R
	    = p.xdr_overall
	    * (p.emerge[ART1_NEW][S1 ] * z[Y4ES1] + p.emerge[ART1_OLD][S1 ] * z[Y4LS1] +
	       p.emerge[ART1_NEW][WS1] * z[Y4Es1] + p.emerge[ART1_OLD][WS1] * z[Y4Ls1]);

	  d[Y4LR2S] = auptake[L2] * (z[Y4NR2S] + z[Y4Nr2S]) + emergeR2S - (p.prog[R2][L2] + p.art_exitL + sec_init_dr) * z[Y4LR2S];
	  d[Y4NR2S] = p.art_exitL * z[Y4LR2S] - (p.prog[R2][L2] + p.revert[NONE][R2] + auptake[L2]) * z[Y4NR2S];
	  d[Y4Nr2S] = p.revert[NONE][R2] * z[Y4NR2S] - (p.prog[WR2][L2] + auptake[L2]) * z[Y4Nr2S];
	  d[Y4LR2R] = auptake[L2] * (z[Y4NR2R] + z[Y4Nr2R]) + emergeR2R - (p.prog[R2][L2] + p.art_exitL + sec_init_dr) * z[Y4LR2R];
	  d[Y4NR2R] = p.art_exitL * z[Y4LR2R] - (p.prog[R2][L2] + p.revert[NONE][R2] + auptake[L2]) * z[Y4NR2R];
	  d[Y4Nr2R] = p.revert[NONE][R2] * z[Y4NR2R] - (p.prog[WR2][L2] + auptake[L2]) * z[Y4Nr2R];

	  d[Y4NC2S] = p.art_exitL * z[Y4LC2S] - (p.prog[C2][L2] + p.revert[NONE][C2] + auptake[L2]) * z[Y4NC2S];
	  d[Y4LC2S] = auptake[L2] * (z[Y4NC2S] + z[Y4Nc2S]) + emergeC2S - (p.prog[C2][L2] + p.art_exitL + sec_init_dr) * z[Y4LC2S];
	  d[Y4Nc2S] = p.revert[NONE][C2] * z[Y4NC2S] - (p.prog[WC2][L2] + auptake[L2]) * z[Y4Nc2S];
	  d[Y4NC2R] = p.art_exitL * z[Y4LC2R] - (p.prog[C2][L2] + p.revert[NONE][C2] + auptake[L2]) * z[Y4NC2R];
	  d[Y4LC2R] = auptake[L2] * (z[Y4NC2R] + z[Y4Nc2R]) + emergeC2R - (p.prog[C2][L2] + p.art_exitL + sec_init_dr) * z[Y4LC2R];
	  d[Y4Nc2R] = p.revert[NONE][C2] * z[Y4NC2R] - (p.prog[WC2][L2] + auptake[L2]) * z[Y4Nc2R];

	  d[Y4Nr1] = p.prog[WR1][L1] * z[Y3Nr1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4Hr1] + z[Y4Pr1]) + p.revert[NONE][R1] * z[Y4NR1] + p.art_exitE * z[Y4Er1] + p.art_exitL * z[Y4Lr1] + p.art_exitF * z[Y4Fr1] - (p.prog[WR1][L2] + auptake[L2]) * z[Y4Nr1];
	  d[Y4Hr1] = p.prog[WR1][L1] * z[Y3Hr1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Pr1] + p.revert[PREP_HIGH][R1] * z[Y4HR1] - (p.prog[WR1][L2] + p.emerge[PREP_HIGH][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4Hr1];
	  d[Y4Pr1] = p.prog[WR1][L1] * z[Y3Pr1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Hr1] + p.revert[PREP_POOR][R1] * z[Y4PR1] - (p.prog[WR1][L2] + p.emerge[PREP_POOR][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4Pr1];
	  d[Y4Er1] = auptake[L2] * z[Y4Nr1] - (p.emerge[ART1_NEW][WR1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4Er1];
	  d[Y4Lr1] = p.art_advance * z[Y4Er1] - (p.emerge[ART1_OLD][WR1] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4Lr1];
	  d[Y4Fr1] = p.revert[ART1_NA][R1] * z[Y4FR1] + p.art_failE * z[Y4Er1] + p.art_failL * z[Y4Lr1] - (p.prog[WR1][L2] + p.art_exitF + sec_init_na) * z[Y4Fr1];
	  d[Y4Tr1] = sec_init_na * z[Y4Fr1] + auptake[L2] * z[Y4Dr1] + p.revert[ART2][R1] * z[Y4TR1] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WR1]) * z[Y4Tr1];
	  d[Y4Dr1] = p.sec_exit * (z[Y4Tr1] + z[Y4Ur1]) + p.revert[NONE][R1] * z[Y4DR1] - (auptake[L2] + p.prog[WR1][L2]) * z[Y4Dr1]; // 2014-10-20 DONE
	  d[Y4Ur1] = p.sec_fail_na * z[Y4Tr1] + p.revert[ART2_NA][R1] * z[Y4UR1] - (p.prog[WR1][L2] + p.sec_exit) * z[Y4Ur1]; // 2014-10-20 DONE

	  d[Y4Nc1] = p.prog[WC1][L1] * z[Y3Nc1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4Hc1] + z[Y4Pc1]) + p.revert[NONE][C1] * z[Y4NC1] + p.art_exitE * z[Y4Ec1] + p.art_exitL * z[Y4Lc1] + p.art_exitF * z[Y4Fc1] - (p.prog[WC1][L2] + auptake[L2]) * z[Y4Nc1];
	  d[Y4Hc1] = p.prog[WC1][L1] * z[Y3Hc1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Pc1] - (p.prog[WC1][L2] + p.emerge[PREP_HIGH][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4Hc1];
	  d[Y4Pc1] = p.prog[WC1][L1] * z[Y3Pc1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Hc1] - (p.prog[WC1][L2] + p.emerge[PREP_POOR][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4Pc1];
	  d[Y4Ec1] = auptake[L2] * z[Y4Nc1] - (p.emerge[ART1_NEW][WC1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4Ec1];
	  d[Y4Lc1] = p.art_advance * z[Y4Ec1] - (p.emerge[ART1_OLD][WC1] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4Lc1];
	  d[Y4Fc1] = p.revert[ART1_NA][C1] * z[Y4FC1] + p.art_failE * z[Y4Ec1] + p.art_failL * z[Y4Lc1] - (p.prog[WC1][L2] + p.art_exitF + sec_init_na) * z[Y4Fc1];
	  d[Y4Tc1] = sec_init_na * z[Y4Fc1] + auptake[L2] * z[Y4Dc1] + p.revert[ART2][C1] * z[Y4TC1] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WC1]) * z[Y4Tc1];
	  d[Y4Dc1] = p.sec_exit * (z[Y4Tc1] + z[Y4Uc1]) + p.revert[NONE][C1] * z[Y4DC1] - (auptake[L2] + p.prog[WC1][L2]) * z[Y4Dc1]; // 2014-10-20 DONE
	  d[Y4Uc1] = p.sec_fail_na * z[Y4Tc1] + p.revert[ART2_NA][C1] * z[Y4UC1] - (p.prog[WC1][L2] + p.sec_exit) * z[Y4Uc1]; // 2014-10-20 DONE

	  d[Y4Nq1] = p.prog[WQ1][L1] * z[Y3Nq1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4Hq1] + z[Y4Pq1]) + p.revert[NONE][Q1] * z[Y4NQ1] + p.art_exitE * z[Y4Eq1] + p.art_exitL * z[Y4Lq1] + p.art_exitF * z[Y4Fq1] - (p.prog[WQ1][L2] + auptake[L2]) * z[Y4Nq1];
	  d[Y4Hq1] = p.prog[WQ1][L1] * z[Y3Hq1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Pq1] - (p.prog[WQ1][L2] + p.emerge[PREP_HIGH][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4Hq1];
	  d[Y4Pq1] = p.prog[WQ1][L1] * z[Y3Pq1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Hq1] - (p.prog[WQ1][L2] + p.emerge[PREP_POOR][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4Pq1];
	  d[Y4Eq1] = auptake[L2] * z[Y4Nq1] - (p.emerge[ART1_NEW][WQ1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4Eq1];
	  d[Y4Lq1] = p.art_advance * z[Y4Eq1] - (p.emerge[ART1_OLD][WQ1] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4Lq1];
	  d[Y4Fq1] = p.revert[ART1_NA][Q1] * z[Y4FQ1] + p.art_failE * z[Y4Eq1] + p.art_failL * z[Y4Lq1] - (p.prog[WQ1][L2] + p.art_exitF + sec_init_na) * z[Y4Fq1];
	  d[Y4Tq1] = sec_init_na * z[Y4Fq1] + auptake[L2] * z[Y4Dq1] + p.revert[ART2][Q1] * z[Y4TQ1] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WQ1]) * z[Y4Tq1];
	  d[Y4Dq1] = p.sec_exit * (z[Y4Tq1] + z[Y4Uq1]) + p.revert[NONE][Q1] * z[Y4DQ1] - (auptake[L2] + p.prog[WQ1][L2]) * z[Y4Dq1]; // 2014-10-20 DONE
	  d[Y4Uq1] = p.sec_fail_na * z[Y4Tq1] + p.revert[ART2_NA][Q1] * z[Y4UQ1] - (p.prog[WQ1][L2] + p.sec_exit) * z[Y4Uq1]; // 2014-10-20 DONE

 	  d[Y4Ns1] = p.prog[WS1][L1] * z[Y3Ns1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4Hs1] + z[Y4Ps1]) + p.revert[NONE][S1] * z[Y4NS1] + p.art_exitE * z[Y4Es1] + p.art_exitL * z[Y4Ls1] + p.art_exitF * z[Y4Fs1] - (p.prog[WS1][L2] + auptake[L2]) * z[Y4Ns1];
 	  d[Y4Hs1] = p.prog[WS1][L1] * z[Y3Hs1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Ps1] + p.revert[PREP_HIGH][S1] * z[Y4HS1] - (p.prog[WS1][L2] + p.emerge[PREP_HIGH][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4Hs1];
 	  d[Y4Ps1] = p.prog[WS1][L1] * z[Y3Ps1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Hs1] + p.revert[PREP_POOR][S1] * z[Y4PS1] - (p.prog[WS1][L2] + p.emerge[PREP_POOR][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4Ps1];
 	  d[Y4Es1] = auptake[L2] * z[Y4Ns1] - (p.emerge[ART1_NEW][WS1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4Es1];
 	  d[Y4Ls1] = p.art_advance * z[Y4Es1] - (p.emerge[ART1_OLD][WS1] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4Ls1];
 	  d[Y4Fs1] = p.revert[ART1_NA][S1] * z[Y4FS1] + p.art_failE * z[Y4Es1] + p.art_failL * z[Y4Ls1] - (p.prog[WS1][L2] + p.art_exitF + sec_init_na) * z[Y4Fs1];
	  d[Y4Ts1] = sec_init_na * z[Y4Fs1] + auptake[L2] * z[Y4Ds1] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WS1]) * z[Y4Ts1]; // 2014-10-20 DONE
	  d[Y4Ds1] = p.sec_exit * (z[Y4Ts1] + z[Y4Us1]) + p.revert[NONE][S1] * z[Y4DS1] - (auptake[L2] + p.prog[WS1][L2]) * z[Y4Ds1]; // 2014-10-20 DONE
	  d[Y4Us1] = p.sec_fail_na * z[Y4Ts1] + p.revert[ART2_NA][S1] * z[Y4US1] - (p.prog[WS1][L2] + p.sec_exit) * z[Y4Us1]; // 2014-10-20 DONE

	  d[Y4Nq2] = p.prog[WQ2][L1] * z[Y3Nq2] + p.revert[NONE][Q2] * z[Y4NQ2] + p.art_exitE * z[Y4Eq2] + p.art_exitL * z[Y4Lq2] + p.art_exitF * z[Y4Fq2] - (p.prog[WQ2][L2] + auptake[L2]) * z[Y4Nq2];
	  d[Y4Eq2] = auptake[L2] * z[Y4Nq2] - (p.emerge[ART1_NEW][WQ2] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L2]) * z[Y4Eq2];
	  d[Y4Lq2] = p.art_advance * z[Y4Eq2] - (p.emerge[ART1_OLD][WQ2] + p.art_failL + p.art_exitL + p.art_mortL[L2]) * z[Y4Lq2];
	  d[Y4Fq2] = p.revert[ART1_NA][Q2] * z[Y4FQ2] + p.art_failE * z[Y4Eq2] + p.art_failL * z[Y4Lq2] - (p.prog[WQ2][L2] + p.art_exitF + sec_init_na) * z[Y4Fq2];
	  d[Y4Tq2] = sec_init_na * z[Y4Fq2] + auptake[L2] * z[Y4Dq2] + p.revert[ART2][Q2] * z[Y4TQ2] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WQ2]) * z[Y4Tq2];
	  d[Y4Dq2] = p.sec_exit * (z[Y4Tq2] + z[Y4Uq2]) + p.revert[NONE][Q2] * z[Y4DQ2] - (auptake[L2] + p.prog[WQ2][L2]) * z[Y4Dq2]; // 2014-10-20 DONE
	  d[Y4Uq2] = p.sec_fail_na * z[Y4Tq2] + p.revert[ART2_NA][Q2] * z[Y4UQ2] - (p.prog[WQ2][L2] + p.sec_exit) * z[Y4Uq2]; // 2014-10-20 DONE

 	  d[Y4TR2S] = sec_init_dr * z[Y4LR2S] + auptake[L2] * z[Y4DR2S] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][R2 ] + p.revert[ART2][R2]) * z[Y4TR2S];
 	  d[Y4TR2R] = sec_init_dr * z[Y4LR2R] + auptake[L2] * z[Y4DR2R] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][WS1] + p.revert[ART2][R2]) * z[Y4TR2R];
 	  d[Y4TC2S] = sec_init_dr * z[Y4LC2S] + auptake[L2] * z[Y4DC2S] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][C2 ] + p.revert[ART2][C2]) * z[Y4TC2S];
 	  d[Y4TC2R] = sec_init_dr * z[Y4LC2R] + auptake[L2] * z[Y4DC2R] - (p.sec_mort[L2] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][WS1] + p.revert[ART2][C2]) * z[Y4TC2R];

 	  d[Y4Tr2S] = auptake[L2] * z[Y4Dr2S] + p.revert[ART2][R2] * z[Y4TR2S] - (p.sec_mort[L2] + p.sec_fail_dr + p.emerge[ART2][WR2] + p.sec_exit) * z[Y4Tr2S];
	  d[Y4Tr2R] = auptake[L2] * z[Y4Dr2R] + p.revert[ART2][R2] * z[Y4TR2R] - (p.sec_mort[L2] + p.sec_fail_dr + p.emerge[ART2][WS1] + p.sec_exit) * z[Y4Tr2R];
 	  d[Y4Tc2S] = auptake[L2] * z[Y4Dc2S] + p.revert[ART2][C2] * z[Y4TC2S] - (p.sec_mort[L2] + p.sec_fail_dr + p.emerge[ART2][WC2] + p.sec_exit) * z[Y4Tc2S];
 	  d[Y4Tc2R] = auptake[L2] * z[Y4Dc2R] + p.revert[ART2][C2] * z[Y4TC2R] - (p.sec_mort[L2] + p.sec_fail_dr + p.emerge[ART2][WS1] + p.sec_exit) * z[Y4Tc2R];

 	  d[Y4DR2S] = p.sec_exit * (z[Y4TR2S] + z[Y4UR2S]) - (auptake[L2] + p.prog[R2][L2] + p.revert[NONE][R2]) * z[Y4DR2S];
 	  d[Y4DR2R] = p.sec_exit * (z[Y4TR2R] + z[Y4UR2R]) - (auptake[L2] + p.prog[R2][L2] + p.revert[NONE][R2]) * z[Y4DR2R];
 	  d[Y4DC2S] = p.sec_exit * (z[Y4TC2S] + z[Y4UC2S]) - (auptake[L2] + p.prog[C2][L2] + p.revert[NONE][C2]) * z[Y4DC2S];
 	  d[Y4DC2R] = p.sec_exit * (z[Y4TC2R] + z[Y4UC2R]) - (auptake[L2] + p.prog[C2][L2] + p.revert[NONE][C2]) * z[Y4DC2R];

	  d[Y4Dr2S] = p.sec_exit * (z[Y4Tr2S] + z[Y4Ur2S]) + p.revert[NONE][R2] * z[Y4DR2S] - (auptake[L2] + p.prog[WR2][L2]) * z[Y4Dr2S];
	  d[Y4Dr2R] = p.sec_exit * (z[Y4Tr2R] + z[Y4Ur2R]) + p.revert[NONE][R2] * z[Y4DR2R] - (auptake[L2] + p.prog[WR2][L2]) * z[Y4Dr2R];
	  d[Y4Dc2S] = p.sec_exit * (z[Y4Tc2S] + z[Y4Uc2S]) + p.revert[NONE][C2] * z[Y4DC2S] - (auptake[L2] + p.prog[WC2][L2]) * z[Y4Dc2S];
	  d[Y4Dc2R] = p.sec_exit * (z[Y4Tc2R] + z[Y4Uc2R]) + p.revert[NONE][C2] * z[Y4DC2R] - (auptake[L2] + p.prog[WC2][L2]) * z[Y4Dc2R];

 	  d[Y4UR2S] = p.sec_fail_dr * z[Y4TR2S] - (p.prog[R2][L2] + p.sec_exit + p.revert[ART2_NA][R2]) * z[Y4UR2S]; // 2014-10-20 DONE
 	  d[Y4UR2R] = p.sec_fail_dr * z[Y4TR2R] - (p.prog[R2][L2] + p.sec_exit + p.revert[ART2_NA][R2]) * z[Y4UR2R]; // 2014-10-20 DONE
 	  d[Y4UC2S] = p.sec_fail_dr * z[Y4TC2S] - (p.prog[C2][L2] + p.sec_exit + p.revert[ART2_NA][C2]) * z[Y4UC2S]; // 2014-10-20 DONE
 	  d[Y4UC2R] = p.sec_fail_dr * z[Y4TC2R] - (p.prog[C2][L2] + p.sec_exit + p.revert[ART2_NA][C2]) * z[Y4UC2R]; // 2014-10-20 DONE

	  d[Y4Ur2S] = p.sec_fail_dr * z[Y4Tr2S] + p.revert[ART2_NA][R2] * z[Y4UR2S] - (p.prog[WR2][L2] + p.sec_exit) * z[Y4Ur2S]; // 2014-10-20 DONE
	  d[Y4Ur2R] = p.sec_fail_dr * z[Y4Tr2R] + p.revert[ART2_NA][R2] * z[Y4UR2R] - (p.prog[WR2][L2] + p.sec_exit) * z[Y4Ur2R]; // 2014-10-20 DONE
	  d[Y4Uc2S] = p.sec_fail_dr * z[Y4Tc2S] + p.revert[ART2_NA][C2] * z[Y4UC2S] - (p.prog[WC2][L2] + p.sec_exit) * z[Y4Uc2S]; // 2014-10-20 DONE
	  d[Y4Uc2R] = p.sec_fail_dr * z[Y4Tc2R] + p.revert[ART2_NA][C2] * z[Y4UC2R] - (p.prog[WC2][L2] + p.sec_exit) * z[Y4Uc2R]; // 2014-10-20 DONE

	  // 2014-10-20 DONE
	  emergeS2 = p.emerge[ART2][WT] * z[Y4TWT]
	    + p.emerge[ART2][R1] * z[Y4TR1] + p.emerge[ART2][WR1] * z[Y4Tr1]
	    + p.emerge[ART2][C1] * z[Y4TC1] + p.emerge[ART2][WC1] * z[Y4Tc1]
	    + p.emerge[ART2][S1] * z[Y4TS1] + p.emerge[ART2][WS1] * z[Y4Ts1]
	    + p.emerge[ART2][Q1] * z[Y4TQ1] + p.emerge[ART2][WQ1] * z[Y4Tq1]
	    + p.emerge[ART2][Q2] * z[Y4TQ2] + p.emerge[ART2][WQ2] * z[Y4Tq2]
	    + p.emerge[ART2][ R2] * z[Y4TR2S] + p.emerge[ART2][ C2] * z[Y4TC2S]
	    + p.emerge[ART2][WR2] * z[Y4Tr2S] + p.emerge[ART2][WC2] * z[Y4Tc2S]
	    + p.emerge[ART2][WS1] * z[Y4TR2R] + p.emerge[ART2][WS1] * z[Y4TC2R]
	    + p.emerge[ART2][WS1] * z[Y4Tr2R] + p.emerge[ART2][WS1] * z[Y4Tc2R];

 	  d[Y4TS2] = auptake[L2] * (z[Y4DS2] + z[Y4Ds2]) + emergeS2 - (p.prog[S2][L2] + p.sec_exit) * z[Y4TS2];
	  d[Y4DS2] = p.sec_exit * z[Y4TS2] - (auptake[L2] + p.prog[S2][L2] + p.revert[NONE][S2]) * z[Y4DS2];
	  d[Y4Ds2] = p.revert[NONE][S2] * z[Y4DS2] - (auptake[L2] + p.prog[WS2][L2]) * z[Y4Ds2];

	  // update cumulative HIV-related deaths with those occuring with 350<CD4<500 (on ART; individuals off ART progress to later stages)
	  mort
	    = p.art_mortE[L2] * (z[Y4EWT] + z[Y4ER1] + z[Y4EC1] + z[Y4ES1] + z[Y4EQ1] + z[Y4EQ2] + z[Y4Er1] + z[Y4Ec1] + z[Y4Es1] + z[Y4Eq1] + z[Y4Eq2])
	    + p.art_mortL[L2] * (z[Y4LWT] + z[Y4LR1] + z[Y4LC1] + z[Y4LS1] + z[Y4LQ1] + z[Y4LQ2] + z[Y4Lr1] + z[Y4Lc1] + z[Y4Ls1] + z[Y4Lq1] + z[Y4Lq2])
	    + p.sec_mort[L2] * (z[Y4TR2S] + z[Y4TR2R] + z[Y4TC2S] + z[Y4TC2R] + z[Y4Tr2S] + z[Y4Tr2R] + z[Y4Tc2S] + z[Y4Tc2R])
	    + p.sec_mort[L2] * (z[Y4TWT] + z[Y4TR1] + z[Y4TC1] + z[Y4TS1] + z[Y4TQ1] + z[Y4TQ2] + z[Y4Tr1] + z[Y4Tc1] + z[Y4Ts1] + z[Y4Tq1] + z[Y4Tq2]); // 2014-10-20 DONE

	  if (b < BANDS - 1) {
	    dy[ACTIVE_HIV_MORT] += mort;
	    dy[ACTIVE_W_HIV_MORT] += mort * (g == F);
	    dy[ACTIVE_M_HIV_MORT] += mort * (g != F);

	    dy[ACTDSC_HIV_MORT] += discount * mort;
	    dy[ACTDSC_W_HIV_MORT] += discount * mort * (g == F);
	    dy[ACTDSC_M_HIV_MORT] += discount * mort * (g != F);
	  }

	  if (u == ANTE) {
	    dy[COHORT_HIV_MORT] += mort;
	    dy[COHORT_W_HIV_MORT] += mort * (g == F);
	    dy[COHORT_M_HIV_MORT] += mort * (g != F);

	    dy[COHDSC_HIV_MORT] += discount * mort;
	    dy[COHDSC_W_HIV_MORT] += discount * mort * (g == F);
	    dy[COHDSC_M_HIV_MORT] += discount * mort * (g != F);
	  }
#else
	  // This code is used if CD4_500 is undefined. In this case,
	  // chronic infection is stratified into just two stages:
	  // CD4>350 and 200<CD4<350. This block of code implements this
	  // by bypassing compartments for CD4>500. Stage durations are
	  // added together elsewhere to keep survival consistent whether
	  // or not CD4_500 is defined.

	  // +-+ chronic infection (350<CD4) +-----------------------------------+
	  d[Y4NWT] = p.prog[WT][P2] * z[Y2NWT] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HWT] + z[Y4PWT]) - (p.prog[WT][L2]) * z[Y4NWT];
	  d[Y4HWT] = p.prog[WT][P2] * z[Y2HWT] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PWT] - (p.prog[WT][L2] + p.emerge[PREP_HIGH][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HWT];
	  d[Y4PWT] = p.prog[WT][P2] * z[Y2PWT] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HWT] - (p.prog[WT][L2] + p.emerge[PREP_POOR][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PWT];

	  d[Y4NR1] = p.prog[R1][P2] * z[Y2NR1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HR1] + z[Y4PR1]) - (p.prog[R1][L2] + p.revert[NONE][R1]) * z[Y4NR1];
	  d[Y4HR1] = p.prog[R1][P2] * z[Y2HR1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PR1] - (p.prog[R1][L2] + p.revert[PREP_HIGH][R1] + p.emerge[PREP_HIGH][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HR1];
	  d[Y4PR1] = p.prog[R1][P2] * z[Y2PR1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HR1] - (p.prog[R1][L2] + p.revert[PREP_POOR][R1] + p.emerge[PREP_POOR][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PR1];

	  d[Y4NC1] = p.prog[C1][P2] * z[Y2NC1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HC1] + z[Y4PC1]) - (p.prog[C1][L2] + p.revert[NONE][C1]) * z[Y4NC1];
	  d[Y4HC1] = p.prog[C1][P2] * z[Y2HC1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PC1] - (p.prog[C1][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HC1];
	  d[Y4PC1] = p.prog[C1][P2] * z[Y2PC1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HC1] - (p.prog[C1][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PC1];

	  d[Y4NQ1] = p.prog[Q1][P2] * z[Y2NQ1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HQ1] + z[Y4PQ1]) - (p.prog[Q1][L2] + p.revert[NONE][Q1]) * z[Y4NQ1];
	  d[Y4HQ1] = p.prog[Q1][P2] * z[Y2HQ1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PQ1] - (p.prog[Q1][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HQ1];
	  d[Y4PQ1] = p.prog[Q1][P2] * z[Y2PQ1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HQ1] - (p.prog[Q1][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PQ1];

 	  d[Y4NS1] = p.prog[S1][P2] * z[Y2NS1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HS1] + z[Y4PS1]) - (p.prog[S1][L2] + p.revert[NONE][S1]) * z[Y4NS1];
 	  d[Y4HS1] = p.prog[S1][P2] * z[Y2HS1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PS1] - (p.prog[S1][L2] + p.revert[PREP_HIGH][S1] + p.emerge[PREP_HIGH][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HS1];
 	  d[Y4PS1] = p.prog[S1][P2] * z[Y2PS1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HS1] - (p.prog[S1][L2] + p.revert[PREP_POOR][S1] + p.emerge[PREP_POOR][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PS1];

	  d[Y4NQ2] = p.prog[Q2][P2] * z[Y2NQ2] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4HQ2] + z[Y4PQ2]) - (p.prog[Q2][L2] + p.revert[NONE][Q2]) * z[Y4NQ2];
	  d[Y4HQ2] = p.prog[Q2][P2] * z[Y2HQ2] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4PQ2] + p.emerge[PREP_HIGH][WT] * z[Y4HWT] + p.emerge[PREP_HIGH][R1] * z[Y4HR1] + p.emerge[PREP_HIGH][WR1] * z[Y4Hr1] + p.emerge[PREP_HIGH][S1] * z[Y4HS1] + p.emerge[PREP_HIGH][WS1] * z[Y4Hs1] + p.emerge[PREP_HIGH][WC1] * z[Y4Hc1] + p.emerge[PREP_HIGH][WQ1] * z[Y4Hq1] - (p.prog[Q2][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4HQ2];
	  d[Y4PQ2] = p.prog[Q2][P2] * z[Y2PQ2] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4HQ2] + p.emerge[PREP_POOR][WT] * z[Y4PWT] + p.emerge[PREP_POOR][R1] * z[Y4PR1] + p.emerge[PREP_POOR][WR1] * z[Y4Pr1] + p.emerge[PREP_POOR][S1] * z[Y4PS1] + p.emerge[PREP_POOR][WS1] * z[Y4Ps1] + p.emerge[PREP_POOR][WC1] * z[Y4Pc1] + p.emerge[PREP_POOR][WQ1] * z[Y4Pq1] - (p.prog[Q2][L2] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4PQ2];

	  d[Y4Nr1] = p.prog[WR1][P2] * z[Y2Nr1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4Hr1] + z[Y4Pr1]) + p.revert[NONE][R1] * z[Y4NR1] - (p.prog[WR1][L2]) * z[Y4Nr1];
	  d[Y4Hr1] = p.prog[WR1][P2] * z[Y2Hr1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Pr1] + p.revert[PREP_HIGH][R1] * z[Y4HR1] - (p.prog[WR1][L2] + p.emerge[PREP_HIGH][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4Hr1];
	  d[Y4Pr1] = p.prog[WR1][P2] * z[Y2Pr1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Hr1] + p.revert[PREP_POOR][R1] * z[Y4PR1] - (p.prog[WR1][L2] + p.emerge[PREP_POOR][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4Pr1];

	  d[Y4Nc1] = p.prog[WC1][P2] * z[Y2Nc1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4Hc1] + z[Y4Pc1]) + p.revert[NONE][C1] * z[Y4NC1] - (p.prog[WC1][L2]) * z[Y4Nc1];
	  d[Y4Hc1] = p.prog[WC1][P2] * z[Y2Hc1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Pc1] - (p.prog[WC1][L2] + p.emerge[PREP_HIGH][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4Hc1];
	  d[Y4Pc1] = p.prog[WC1][P2] * z[Y2Pc1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Hc1] - (p.prog[WC1][L2] + p.emerge[PREP_POOR][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4Pc1];

	  d[Y4Nq1] = p.prog[WQ1][P2] * z[Y2Nq1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4Hq1] + z[Y4Pq1]) + p.revert[NONE][Q1] * z[Y4NQ1] - (p.prog[WQ1][L2]) * z[Y4Nq1];
	  d[Y4Hq1] = p.prog[WQ1][P2] * z[Y2Hq1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Pq1] - (p.prog[WQ1][L2] + p.emerge[PREP_HIGH][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4Hq1];
	  d[Y4Pq1] = p.prog[WQ1][P2] * z[Y2Pq1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Hq1] - (p.prog[WQ1][L2] + p.emerge[PREP_POOR][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4Pq1];

	  d[Y4Ns1] = p.prog[WS1][P2] * z[Y2Ns1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y4Hs1] + z[Y4Ps1]) + p.revert[NONE][S1] * z[Y4NS1] - (p.prog[WS1][L2]) * z[Y4Ns1];
	  d[Y4Hs1] = p.prog[WS1][P2] * z[Y2Hs1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Ps1] + p.revert[PREP_HIGH][S1] * z[Y4HS1] - (p.prog[WS1][L2] + p.emerge[PREP_HIGH][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y4Hs1];
	  d[Y4Ps1] = p.prog[WS1][P2] * z[Y2Ps1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y4Hs1] + p.revert[PREP_POOR][S1] * z[Y4PS1] - (p.prog[WS1][L2] + p.emerge[PREP_POOR][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y4Ps1];

	  d[Y4Nq2] = p.prog[WQ2][P2] * z[Y2Nq2] + p.revert[NONE][Q2] * z[Y4NQ2] - (p.prog[WQ2][L2]) * z[Y4Nq2];
#endif // CD4_500
	  // +-+ chronic infection (200<CD4<350) +-------------------------------+
	  d[Y5NWT] = p.prog[WT][L2] * z[Y4NWT] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y5HWT] + z[Y5PWT]) + p.art_exitE * z[Y5EWT] + p.art_exitL * z[Y5LWT] + p.art_exitF * z[Y5FWT] - (p.prog[WT][L3] + auptake[L3]) * z[Y5NWT];
	  d[Y5HWT] = p.prog[WT][L2] * z[Y4HWT] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5PWT] - (p.prog[WT][L3] + p.emerge[PREP_HIGH][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y5HWT];
	  d[Y5PWT] = p.prog[WT][L2] * z[Y4PWT] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5HWT] - (p.prog[WT][L3] + p.emerge[PREP_POOR][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y5PWT];
	  d[Y5EWT] = auptake[L3] * z[Y5NWT] - (p.emerge[ART1_NEW][WT] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5EWT];
	  d[Y5LWT] = p.art_advance * z[Y5EWT] - (p.emerge[ART1_OLD][WT] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5LWT];
	  d[Y5FWT] = p.prog[WT][L2] * z[Y4FWT] + p.art_failE * z[Y5EWT] + p.art_failL * z[Y5LWT] - (p.prog[WT][L3] + p.art_exitF + sec_init_na) * z[Y5FWT];
	  d[Y5TWT] = sec_init_na * z[Y5FWT] + auptake[L3] * z[Y5DWT] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WT]) * z[Y5TWT]; // 2014-10-20 DONE
 	  d[Y5DWT] = p.prog[WT][L2] * z[Y4DWT] + p.sec_exit * (z[Y5TWT] + z[Y5UWT]) - (auptake[L3] + p.prog[WT][L3]) * z[Y5DWT]; // 2014-10-20 DONE
 	  d[Y5UWT] = p.prog[WT][L2] * z[Y4UWT] + p.sec_fail_na * z[Y5TWT] - (p.prog[WT][L3] + p.sec_exit) * z[Y5UWT]; // 2014-10-20 DONE

	  d[Y5NR1] = p.prog[R1][L2] * z[Y4NR1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y5HR1] + z[Y5PR1]) + p.art_exitE * z[Y5ER1] + p.art_exitL * z[Y5LR1] + p.art_exitF * z[Y5FR1] - (p.prog[R1][L3] + p.revert[NONE][R1] + auptake[L3]) * z[Y5NR1];
	  d[Y5HR1] = p.prog[R1][L2] * z[Y4HR1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5PR1] - (p.prog[R1][L3] + p.revert[PREP_HIGH][R1] + p.emerge[PREP_HIGH][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y5HR1];
	  d[Y5PR1] = p.prog[R1][L2] * z[Y4PR1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5HR1] - (p.prog[R1][L3] + p.revert[PREP_POOR][R1] + p.emerge[PREP_POOR][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y5PR1];
	  d[Y5ER1] = auptake[L3] * z[Y5NR1] - (p.emerge[ART1_NEW][R1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5ER1];
	  d[Y5LR1] = p.art_advance * z[Y5ER1] - (p.emerge[ART1_OLD][R1] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5LR1];
	  d[Y5FR1] = p.prog[R1][L2] * z[Y4FR1] + p.art_failE * z[Y5ER1] + p.art_failL * z[Y5LR1] - (p.prog[R1][L3] + p.revert[ART1_NA][R1] + p.art_exitF + sec_init_na) * z[Y5FR1];
	  d[Y5TR1] = sec_init_na * z[Y5FR1] + auptake[L3] * z[Y5DR1] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.revert[ART2][R1] + p.emerge[ART2][R1]) * z[Y5TR1];
 	  d[Y5DR1] = p.prog[R1][L2] * z[Y4DR1] + p.sec_exit * (z[Y5TR1] + z[Y5UR1]) - (auptake[L3] + p.prog[R1][L3] + p.revert[NONE][R1]) * z[Y5DR1]; // 2014-10-20 DONE
 	  d[Y5UR1] = p.prog[R1][L2] * z[Y4UR1] + p.sec_fail_na * z[Y5TR1] - (p.prog[R1][L3] + p.revert[ART2_NA][R1] + p.sec_exit) * z[Y5UR1]; // 2014-10-20 DONE

	  d[Y5NC1] = p.prog[C1][L2] * z[Y4NC1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y5HC1] + z[Y5PC1]) + p.art_exitE * z[Y5EC1] + p.art_exitL * z[Y5LC1] + p.art_exitF * z[Y5FC1] - (p.prog[C1][L3] + p.revert[NONE][C1] + auptake[L3]) * z[Y5NC1];
	  d[Y5HC1] = p.prog[C1][L2] * z[Y4HC1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5PC1] - (p.prog[C1][L3] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y5HC1];
	  d[Y5PC1] = p.prog[C1][L2] * z[Y4PC1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5HC1] - (p.prog[C1][L3] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y5PC1];
	  d[Y5EC1] = auptake[L3] * z[Y5NC1] - (p.emerge[ART1_NEW][C1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5EC1];
	  d[Y5LC1] = p.art_advance * z[Y5EC1] - (p.emerge[ART1_OLD][C1] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5LC1];
	  d[Y5FC1] = p.prog[C1][L2] * z[Y4FC1] + p.art_failE * z[Y5EC1] + p.art_failL * z[Y5LC1] - (p.prog[C1][L3] + p.revert[ART1_NA][C1] + p.art_exitF + sec_init_na) * z[Y5FC1];
	  d[Y5TC1] = sec_init_na * z[Y5FC1] + auptake[L3] * z[Y5DC1] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.revert[ART2][C1] + p.emerge[ART2][C1]) * z[Y5TC1];
 	  d[Y5DC1] = p.prog[C1][L2] * z[Y4DC1] + p.sec_exit * (z[Y5TC1] + z[Y5UC1]) - (auptake[L3] + p.prog[C1][L3] + p.revert[NONE][C1]) * z[Y5DC1]; // 2014-10-20 DONE
 	  d[Y5UC1] = p.prog[C1][L2] * z[Y4UC1] + p.sec_fail_na * z[Y5TC1] - (p.prog[C1][L3] + p.revert[ART2_NA][C1] + p.sec_exit) * z[Y5UC1]; // 2014-10-20 DONE

	  d[Y5NQ1] = p.prog[Q1][L2] * z[Y4NQ1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y5HQ1] + z[Y5PQ1]) + p.art_exitE * z[Y5EQ1] + p.art_exitL * z[Y5LQ1] + p.art_exitF * z[Y5FQ1] - (p.prog[Q1][L3] + p.revert[NONE][Q1] + auptake[L3]) * z[Y5NQ1];
	  d[Y5HQ1] = p.prog[Q1][L2] * z[Y4HQ1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5PQ1] - (p.prog[Q1][L3] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y5HQ1];
	  d[Y5PQ1] = p.prog[Q1][L2] * z[Y4PQ1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5HQ1] - (p.prog[Q1][L3] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y5PQ1];
	  d[Y5EQ1] = auptake[L3] * z[Y5NQ1] - (p.emerge[ART1_NEW][Q1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5EQ1];
	  d[Y5LQ1] = p.art_advance * z[Y5EQ1] - (p.emerge[ART1_OLD][Q1] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5LQ1];
	  d[Y5FQ1] = p.prog[Q1][L2] * z[Y4FQ1] + p.art_failE * z[Y5EQ1] + p.art_failL * z[Y5LQ1] - (p.prog[Q1][L3] + p.revert[ART1_NA][Q1] + p.art_exitF + sec_init_na) * z[Y5FQ1];
	  d[Y5TQ1] = sec_init_na * z[Y5FQ1] + auptake[L3] * z[Y5DQ1] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.revert[ART2][Q1] + p.emerge[ART2][Q1]) * z[Y5TQ1];
 	  d[Y5DQ1] = p.prog[Q1][L2] * z[Y4DQ1] + p.sec_exit * (z[Y5TQ1] + z[Y5UQ1]) - (auptake[L3] + p.prog[Q1][L3] + p.revert[NONE][Q1]) * z[Y5DQ1]; // 2014-10-20 DONE
 	  d[Y5UQ1] = p.prog[Q1][L2] * z[Y4UQ1] + p.sec_fail_na * z[Y5TQ1] - (p.prog[Q1][L3] + p.revert[ART2_NA][Q1] + p.sec_exit) * z[Y5UQ1]; // 2014-10-20 DONE

	  d[Y5NS1] = p.prog[S1][L2] * z[Y4NS1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y5HS1] + z[Y5PS1]) + p.art_exitE * z[Y5ES1] + p.art_exitL * z[Y5LS1] + p.art_exitF * z[Y5FS1] - (p.prog[S1][L3] + p.revert[NONE][S1] + auptake[L3]) * z[Y5NS1];
	  d[Y5HS1] = p.prog[S1][L2] * z[Y4HS1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5PS1] - (p.prog[S1][L3] + p.revert[PREP_HIGH][S1] + p.emerge[PREP_HIGH][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y5HS1];
	  d[Y5PS1] = p.prog[S1][L2] * z[Y4PS1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5HS1] - (p.prog[S1][L3] + p.revert[PREP_POOR][S1] + p.emerge[PREP_POOR][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y5PS1];
	  d[Y5ES1] = auptake[L3] * z[Y5NS1] - (p.emerge[ART1_NEW][S1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5ES1];
	  d[Y5LS1] = p.art_advance * z[Y5ES1] - (p.emerge[ART1_OLD][S1] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5LS1];
	  d[Y5FS1] = p.prog[S1][L2] * z[Y4FS1] + p.art_failE * z[Y5ES1] + p.art_failL * z[Y5LS1] - (p.prog[S1][L3] + p.revert[ART1_NA][S1] + p.art_exitF + sec_init_na) * z[Y5FS1];
 	  d[Y5TS1] = sec_init_na * z[Y5FS1] + auptake[L3] * z[Y5DS1] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][S1]) * z[Y5TS1]; // 2014-10-20 DONE
 	  d[Y5DS1] = p.prog[S1][L2] * z[Y4DS1] + p.sec_exit * (z[Y5TS1] + z[Y5US1]) - (auptake[L3] + p.prog[S1][L3] + p.revert[NONE][S1]) * z[Y5DS1]; // 2014-10-20 DONE
 	  d[Y5US1] = p.prog[S1][L2] * z[Y4US1] + p.sec_fail_na * z[Y5TS1] - (p.prog[S1][L3] + p.revert[ART2_NA][S1] + p.sec_exit) * z[Y5US1]; // 2014-10-20 DONE

	  d[Y5NQ2] = p.prog[Q2][L2] * z[Y4NQ2] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y5HQ2] + z[Y5PQ2]) + p.art_exitE * z[Y5EQ2] + p.art_exitL * z[Y5LQ2] + p.art_exitF * z[Y5FQ2] - (p.prog[Q2][L3] + p.revert[NONE][Q2] + auptake[L3]) * z[Y5NQ2];
	  d[Y5HQ2] = p.prog[Q2][L2] * z[Y4HQ2] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5PQ2] + p.emerge[PREP_HIGH][WT] * z[Y5HWT] + p.emerge[PREP_HIGH][R1] * z[Y5HR1] + p.emerge[PREP_HIGH][WR1] * z[Y5Hr1] + p.emerge[PREP_HIGH][S1] * z[Y5HS1] + p.emerge[PREP_HIGH][WS1] * z[Y5Hs1] + p.emerge[PREP_HIGH][WC1] * z[Y5Hc1] + p.emerge[PREP_HIGH][WQ1] * z[Y5Hq1] - (p.prog[Q2][L3] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y5HQ2];
	  d[Y5PQ2] = p.prog[Q2][L2] * z[Y4PQ2] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5HQ2] + p.emerge[PREP_POOR][WT] * z[Y5PWT] + p.emerge[PREP_POOR][R1] * z[Y5PR1] + p.emerge[PREP_POOR][WR1] * z[Y5Pr1] + p.emerge[PREP_POOR][S1] * z[Y5PS1] + p.emerge[PREP_POOR][WS1] * z[Y5Ps1] + p.emerge[PREP_POOR][WC1] * z[Y5Pc1] + p.emerge[PREP_POOR][WQ1] * z[Y5Pq1] - (p.prog[Q2][L3] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y5PQ2];
	  d[Y5EQ2] = auptake[L3] * z[Y5NQ2] - (p.emerge[ART1_NEW][Q2] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5EQ2];
	  d[Y5LQ2] = p.art_advance * z[Y5EQ2] - (p.emerge[ART1_OLD][Q2] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5LQ2];
	  d[Y5FQ2] = p.prog[Q2][L2] * z[Y4FQ2] + p.art_failE * z[Y5EQ2] + p.art_failL * z[Y5LQ2] - (p.prog[Q2][L3] + p.revert[ART1_NA][Q2] + p.art_exitF + sec_init_na) * z[Y5FQ2];
	  d[Y5TQ2] = sec_init_na * z[Y5FQ2] + auptake[L3] * z[Y5DQ2] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.revert[ART2][Q2] + p.emerge[ART2][Q2]) * z[Y5TQ2];
 	  d[Y5DQ2] = p.prog[Q2][L2] * z[Y4DQ2] + p.sec_exit * (z[Y5TQ2] + z[Y5UQ2]) - (auptake[L3] + p.prog[Q2][L3] + p.revert[NONE][Q2]) * z[Y5DQ2]; // 2014-10-20 DONE
 	  d[Y5UQ2] = p.prog[Q2][L2] * z[Y4UQ2] + p.sec_fail_na * z[Y5TQ2] - (p.prog[Q2][L3] + p.revert[ART2_NA][Q2] + p.sec_exit) * z[Y5UQ2]; // 2014-10-20 DONE

	  emergeR2S = (1.0 - p.xdr_overall) *
	    (p.emerge[ART1_NEW][WT ] * z[Y5EWT] + p.emerge[ART1_OLD][WT ] * z[Y5LWT] +
	     p.emerge[ART1_NEW][R1 ] * z[Y5ER1] + p.emerge[ART1_OLD][R1 ] * z[Y5LR1] +
	     p.emerge[ART1_NEW][WR1] * z[Y5Er1] + p.emerge[ART1_OLD][WR1] * z[Y5Lr1]);
	  emergeR2R = (1.0 - p.xdr_overall) *
	    (p.emerge[ART1_NEW][S1 ] * z[Y5ES1] + p.emerge[ART1_OLD][S1 ] * z[Y5LS1] +
	     p.emerge[ART1_NEW][WS1] * z[Y5Es1] + p.emerge[ART1_OLD][WS1] * z[Y5Ls1]);

	  emergeC2S
	    = p.xdr_overall
	    * (p.emerge[ART1_NEW][WT ] * z[Y5EWT] + p.emerge[ART1_OLD][WT ] * z[Y5LWT] +
	       p.emerge[ART1_NEW][R1 ] * z[Y5ER1] + p.emerge[ART1_OLD][R1 ] * z[Y5LR1] +
	       p.emerge[ART1_NEW][WR1] * z[Y5Er1] + p.emerge[ART1_OLD][WR1] * z[Y5Lr1])
	    + (p.emerge[ART1_NEW][C1 ] * z[Y5EC1] + p.emerge[ART1_OLD][C1 ] * z[Y5LC1] +
	       p.emerge[ART1_NEW][WC1] * z[Y5Ec1] + p.emerge[ART1_OLD][WC1] * z[Y5Lc1] +
	       p.emerge[ART1_NEW][Q1 ] * z[Y5EQ1] + p.emerge[ART1_OLD][Q1 ] * z[Y5LQ1] +
	       p.emerge[ART1_NEW][WQ1] * z[Y5Eq1] + p.emerge[ART1_OLD][WQ1] * z[Y5Lq1] +
	       p.emerge[ART1_NEW][Q2 ] * z[Y5EQ2] + p.emerge[ART1_OLD][Q2 ] * z[Y5LQ2] +
	       p.emerge[ART1_NEW][WQ2] * z[Y5Eq2] + p.emerge[ART1_OLD][WQ2] * z[Y5Lq2]);
	  emergeC2R
	    = p.xdr_overall
	    * (p.emerge[ART1_NEW][S1 ] * z[Y5ES1] + p.emerge[ART1_OLD][S1 ] * z[Y5LS1] +
	       p.emerge[ART1_NEW][WS1] * z[Y5Es1] + p.emerge[ART1_OLD][WS1] * z[Y5Ls1]);

	  d[Y5LR2S] = p.prog[ R2][L2] * z[Y4LR2S] + auptake[L3] * (z[Y5NR2S] + z[Y5Nr2S]) + emergeR2S - (p.prog[R2][L3] + p.art_exitL + sec_init_dr) * z[Y5LR2S];
	  d[Y5NR2S] = p.prog[ R2][L2] * z[Y4NR2S] + p.art_exitL * z[Y5LR2S] - (p.prog[R2][L3] + p.revert[NONE][R2] + auptake[L3]) * z[Y5NR2S];
	  d[Y5Nr2S] = p.prog[WR2][L2] * z[Y4Nr2S] + p.revert[NONE][R2] * z[Y5NR2S] - (p.prog[WR2][L3] + auptake[L3]) * z[Y5Nr2S];
	  d[Y5LR2R] = p.prog[ R2][L2] * z[Y4LR2R] + auptake[L3] * (z[Y5NR2R] + z[Y5Nr2R]) + emergeR2R - (p.prog[R2][L3] + p.art_exitL + sec_init_dr) * z[Y5LR2R];
	  d[Y5NR2R] = p.prog[ R2][L2] * z[Y4NR2R] + p.art_exitL * z[Y5LR2R] - (p.prog[R2][L3] + p.revert[NONE][R2] + auptake[L3]) * z[Y5NR2R];
	  d[Y5Nr2R] = p.prog[WR2][L2] * z[Y4Nr2R] + p.revert[NONE][R2] * z[Y5NR2R] - (p.prog[WR2][L3] + auptake[L3]) * z[Y5Nr2R];

	  d[Y5LC2S] = p.prog[ C2][L2] * z[Y4LC2S] + auptake[L3] * (z[Y5NC2S] + z[Y5Nc2S]) + emergeC2S - (p.prog[C2][L3] + p.art_exitL + sec_init_dr) * z[Y5LC2S];
	  d[Y5NC2S] = p.prog[ C2][L2] * z[Y4NC2S] + p.art_exitL * z[Y5LC2S] - (p.prog[C2][L3] + p.revert[NONE][C2] + auptake[L3]) * z[Y5NC2S];
	  d[Y5Nc2S] = p.prog[WC2][L2] * z[Y4Nc2S] + p.revert[NONE][C2] * z[Y5NC2S] - (p.prog[WC2][L3] + auptake[L3]) * z[Y5Nc2S];
	  d[Y5LC2R] = p.prog[ C2][L2] * z[Y4LC2R] + auptake[L3] * (z[Y5NC2R] + z[Y5Nc2R]) + emergeC2R - (p.prog[C2][L3] + p.art_exitL + sec_init_dr) * z[Y5LC2R];
	  d[Y5NC2R] = p.prog[ C2][L2] * z[Y4NC2R] + p.art_exitL * z[Y5LC2R] - (p.prog[C2][L3] + p.revert[NONE][C2] + auptake[L3]) * z[Y5NC2R];
	  d[Y5Nc2R] = p.prog[WC2][L2] * z[Y4Nc2R] + p.revert[NONE][C2] * z[Y5NC2R] - (p.prog[WC2][L3] + auptake[L3]) * z[Y5Nc2R];

	  d[Y5Nr1] = p.prog[WR1][L2] * z[Y4Nr1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y5Hr1] + z[Y5Pr1]) + p.revert[NONE][R1] * z[Y5NR1] + p.art_exitE * z[Y5Er1] + p.art_exitL * z[Y5Lr1] + p.art_exitF * z[Y5Fr1] - (p.prog[WR1][L3] + auptake[L3]) * z[Y5Nr1];
	  d[Y5Hr1] = p.prog[WR1][L2] * z[Y4Hr1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5Pr1] + p.revert[PREP_HIGH][R1] * z[Y5HR1] - (p.prog[WR1][L3] + p.emerge[PREP_HIGH][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y5Hr1];
	  d[Y5Pr1] = p.prog[WR1][L2] * z[Y4Pr1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5Hr1] + p.revert[PREP_POOR][R1] * z[Y5PR1] - (p.prog[WR1][L3] + p.emerge[PREP_POOR][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y5Pr1];
	  d[Y5Er1] = auptake[L3] * z[Y5Nr1] - (p.emerge[ART1_NEW][WR1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5Er1];
	  d[Y5Lr1] = p.art_advance * z[Y5Er1] - (p.emerge[ART1_OLD][WR1] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5Lr1];
	  d[Y5Fr1] = p.revert[ART1_NA][R1] * z[Y5FR1] + p.prog[WR1][L2] * z[Y4Fr1] + p.art_failE * z[Y5Er1] + p.art_failL * z[Y5Lr1] - (p.prog[WR1][L3] + p.art_exitF + sec_init_na) * z[Y5Fr1];
	  d[Y5Tr1] = sec_init_na * z[Y5Fr1] + auptake[L3] * z[Y5Dr1] + p.revert[ART2][R1] * z[Y5TR1] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WR1]) * z[Y5Tr1];
 	  d[Y5Dr1] = p.prog[WR1][L2] * z[Y4Dr1] + p.sec_exit * (z[Y5Tr1] + z[Y5Ur1]) + p.revert[NONE][R1] * z[Y5DR1] - (auptake[L3] + p.prog[WR1][L3]) * z[Y5Dr1]; // 2014-10-20 DONE
 	  d[Y5Ur1] = p.prog[WR1][L2] * z[Y4Ur1] + p.sec_fail_na * z[Y5Tr1] + p.revert[ART2_NA][R1] * z[Y5UR1] - (p.prog[WR1][L3] + p.sec_exit) * z[Y5Ur1]; // 2014-10-20 DONE

	  d[Y5Nc1] = p.prog[WC1][L2] * z[Y4Nc1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y5Hc1] + z[Y5Pc1]) + p.revert[NONE][C1] * z[Y5NC1] + p.art_exitE * z[Y5Ec1] + p.art_exitL * z[Y5Lc1] + p.art_exitF * z[Y5Fc1] - (p.prog[WC1][L3] + auptake[L3]) * z[Y5Nc1];
	  d[Y5Hc1] = p.prog[WC1][L2] * z[Y4Hc1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5Pc1] - (p.prog[WC1][L3] + p.emerge[PREP_HIGH][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y5Hc1];
	  d[Y5Pc1] = p.prog[WC1][L2] * z[Y4Pc1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5Hc1] - (p.prog[WC1][L3] + p.emerge[PREP_POOR][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y5Pc1];
	  d[Y5Ec1] = auptake[L3] * z[Y5Nc1] - (p.emerge[ART1_NEW][WC1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5Ec1];
	  d[Y5Lc1] = p.art_advance * z[Y5Ec1] - (p.emerge[ART1_OLD][WC1] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5Lc1];
	  d[Y5Fc1] = p.revert[ART1_NA][C1] * z[Y5FC1] + p.prog[WC1][L2] * z[Y4Fc1] + p.art_failE * z[Y5Ec1] + p.art_failL * z[Y5Lc1] - (p.prog[WC1][L3] + p.art_exitF + sec_init_na) * z[Y5Fc1];
	  d[Y5Tc1] = sec_init_na * z[Y5Fc1] + auptake[L3] * z[Y5Dc1] + p.revert[ART2][C1] * z[Y5TC1] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WC1]) * z[Y5Tc1];
	  d[Y5Dc1] = p.prog[WC1][L2] * z[Y4Dc1] + p.sec_exit * (z[Y5Tc1] + z[Y5Uc1]) + p.revert[NONE][C1] * z[Y5DC1] - (auptake[L3] + p.prog[WC1][L3]) * z[Y5Dc1]; // 2014-10-20 DONE
	  d[Y5Uc1] = p.prog[WC1][L2] * z[Y4Uc1] + p.sec_fail_na * z[Y5Tc1] + p.revert[ART2_NA][C1] * z[Y5UC1] - (p.prog[WC1][L3] + p.sec_exit) * z[Y5Uc1]; // 2014-10-20 DONE

	  d[Y5Nq1] = p.prog[WQ1][L2] * z[Y4Nq1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y5Hq1] + z[Y5Pq1]) + p.revert[NONE][Q1] * z[Y5NQ1] + p.art_exitE * z[Y5Eq1] + p.art_exitL * z[Y5Lq1] + p.art_exitF * z[Y5Fq1] - (p.prog[WQ1][L3] + auptake[L3]) * z[Y5Nq1];
	  d[Y5Hq1] = p.prog[WQ1][L2] * z[Y4Hq1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5Pq1] - (p.prog[WQ1][L3] + p.emerge[PREP_HIGH][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y5Hq1];
	  d[Y5Pq1] = p.prog[WQ1][L2] * z[Y4Pq1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5Hq1] - (p.prog[WQ1][L3] + p.emerge[PREP_POOR][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y5Pq1];
	  d[Y5Eq1] = auptake[L3] * z[Y5Nq1] - (p.emerge[ART1_NEW][WQ1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5Eq1];
	  d[Y5Lq1] = p.art_advance * z[Y5Eq1] - (p.emerge[ART1_OLD][WQ1] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5Lq1];
	  d[Y5Fq1] = p.revert[ART1_NA][Q1] * z[Y5FQ1] + p.prog[WQ1][L2] * z[Y4Fq1] + p.art_failE * z[Y5Eq1] + p.art_failL * z[Y5Lq1] - (p.prog[WQ1][L3] + p.art_exitF + sec_init_na) * z[Y5Fq1];
	  d[Y5Tq1] = sec_init_na * z[Y5Fq1] + auptake[L3] * z[Y5Dq1] + p.revert[ART2][Q1] * z[Y5TQ1] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WQ1]) * z[Y5Tq1];
	  d[Y5Dq1] = p.prog[WQ1][L2] * z[Y4Dq1] + p.sec_exit * (z[Y5Tq1] + z[Y5Uq1]) + p.revert[NONE][Q1] * z[Y5DQ1] - (auptake[L3] + p.prog[WQ1][L3]) * z[Y5Dq1]; // 2014-10-20 DONE
	  d[Y5Uq1] = p.prog[WQ1][L2] * z[Y4Uq1] + p.sec_fail_na * z[Y5Tq1] + p.revert[ART2_NA][Q1] * z[Y5UQ1] - (p.prog[WQ1][L3] + p.sec_exit) * z[Y5Uq1]; // 2014-10-20 DONE

	  d[Y5Ns1] = p.prog[WS1][L2] * z[Y4Ns1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y5Hs1] + z[Y5Ps1]) + p.revert[NONE][S1] * z[Y5NS1] + p.art_exitE * z[Y5Es1] + p.art_exitL * z[Y5Ls1] + p.art_exitF * z[Y5Fs1] - (p.prog[WS1][L3] + auptake[L3]) * z[Y5Ns1];
	  d[Y5Hs1] = p.prog[WS1][L2] * z[Y4Hs1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5Ps1] + p.revert[PREP_HIGH][S1] * z[Y5HS1] - (p.prog[WS1][L3] + p.emerge[PREP_HIGH][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y5Hs1];
	  d[Y5Ps1] = p.prog[WS1][L2] * z[Y4Ps1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y5Hs1] + p.revert[PREP_POOR][S1] * z[Y5PS1] - (p.prog[WS1][L3] + p.emerge[PREP_POOR][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y5Ps1];
	  d[Y5Es1] = auptake[L3] * z[Y5Ns1] - (p.emerge[ART1_NEW][WS1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5Es1];
	  d[Y5Ls1] = p.art_advance * z[Y5Es1] - (p.emerge[ART1_OLD][WS1] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5Ls1];
	  d[Y5Fs1] = p.revert[ART1_NA][S1] * z[Y5FS1] + p.prog[WS1][L2] * z[Y4Fs1] + p.art_failE * z[Y5Es1] + p.art_failL * z[Y5Ls1] - (p.prog[WS1][L3] + p.art_exitF + sec_init_na) * z[Y5Fs1];
     	  d[Y5Ts1] = sec_init_na * z[Y5Fs1] + auptake[L3] * z[Y5Ds1] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WS1]) * z[Y5Ts1]; // 2014-10-20 DONE
	  d[Y5Ds1] = p.prog[WS1][L2] * z[Y4Ds1] + p.sec_exit * (z[Y5Ts1] + z[Y5Us1]) + p.revert[NONE][S1] * z[Y5DS1] - (auptake[L3] + p.prog[WS1][L3]) * z[Y5Ds1]; // 2014-10-20 DONE
	  d[Y5Us1] = p.prog[WS1][L2] * z[Y4Us1] + p.sec_fail_na * z[Y5Ts1] + p.revert[ART2_NA][S1] * z[Y5US1] - (p.prog[WS1][L3] + p.sec_exit) * z[Y5Us1]; // 2014-10-20 DONE

	  d[Y5Nq2] = p.prog[WQ2][L2] * z[Y4Nq2] + p.revert[NONE][Q2] * z[Y5NQ2] + p.art_exitE * z[Y5Eq2] + p.art_exitL * z[Y5Lq2] + p.art_exitF * z[Y5Fq2] - (p.prog[WQ2][L3] + auptake[L3]) * z[Y5Nq2];
	  d[Y5Eq2] = auptake[L3] * z[Y5Nq2] - (p.emerge[ART1_NEW][WQ2] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[L3]) * z[Y5Eq2];
	  d[Y5Lq2] = p.art_advance * z[Y5Eq2] - (p.emerge[ART1_OLD][WQ2] + p.art_failL + p.art_exitL + p.art_mortL[L3]) * z[Y5Lq2];
	  d[Y5Fq2] = p.revert[ART1_NA][Q2] * z[Y5FQ2] + p.prog[WQ2][L2] * z[Y4Fq2] + p.art_failE * z[Y5Eq2] + p.art_failL * z[Y5Lq2] - (p.prog[WQ2][L3] + p.art_exitF + sec_init_na) * z[Y5Fq2];
	  d[Y5Tq2] = sec_init_na * z[Y5Fq2] + auptake[L3] * z[Y5Dq2] + p.revert[ART2][Q2] * z[Y5TQ2] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WQ2]) * z[Y5Tq2];
	  d[Y5Dq2] = p.prog[WQ2][L2] * z[Y4Dq2] + p.sec_exit * (z[Y5Tq2] + z[Y5Uq2]) + p.revert[NONE][Q2] * z[Y5DQ2] - (auptake[L3] + p.prog[WQ2][L3]) * z[Y5Dq2]; // 2014-10-20 DONE
	  d[Y5Uq2] = p.prog[WQ2][L2] * z[Y4Uq2] + p.sec_fail_na * z[Y5Tq2] + p.revert[ART2_NA][Q2] * z[Y5UQ2] - (p.prog[WQ2][L3] + p.sec_exit) * z[Y5Uq2]; // 2014-10-20 DONE

	  d[Y5TR2S] = sec_init_dr * z[Y5LR2S] + auptake[L3] * z[Y5DR2S] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][R2 ] + p.revert[ART2][R2]) * z[Y5TR2S];
	  d[Y5TR2R] = sec_init_dr * z[Y5LR2R] + auptake[L3] * z[Y5DR2R] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][WS1] + p.revert[ART2][R2]) * z[Y5TR2R];
	  d[Y5TC2S] = sec_init_dr * z[Y5LC2S] + auptake[L3] * z[Y5DC2S] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][C2 ] + p.revert[ART2][C2]) * z[Y5TC2S];
	  d[Y5TC2R] = sec_init_dr * z[Y5LC2R] + auptake[L3] * z[Y5DC2R] - (p.sec_mort[L3] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][WS1] + p.revert[ART2][C2]) * z[Y5TC2R];

	  d[Y5Tr2S] = auptake[L3] * z[Y5Dr2S] + p.revert[ART2][R2] * z[Y5TR2S] - (p.sec_mort[L3] + p.sec_fail_dr + p.emerge[ART2][WR2] + p.sec_exit) * z[Y5Tr2S];
	  d[Y5Tr2R] = auptake[L3] * z[Y5Dr2R] + p.revert[ART2][R2] * z[Y5TR2R] - (p.sec_mort[L3] + p.sec_fail_dr + p.emerge[ART2][WS1] + p.sec_exit) * z[Y5Tr2R];
	  d[Y5Tc2S] = auptake[L3] * z[Y5Dc2S] + p.revert[ART2][C2] * z[Y5TC2S] - (p.sec_mort[L3] + p.sec_fail_dr + p.emerge[ART2][WC2] + p.sec_exit) * z[Y5Tc2S];
	  d[Y5Tc2R] = auptake[L3] * z[Y5Dc2R] + p.revert[ART2][C2] * z[Y5TC2R] - (p.sec_mort[L3] + p.sec_fail_dr + p.emerge[ART2][WS1] + p.sec_exit) * z[Y5Tc2R];

	  d[Y5DR2S] = p.sec_exit * (z[Y5TR2S] + z[Y5UR2S]) - (auptake[L3] + p.prog[R2][L3] + p.revert[NONE][R2]) * z[Y5DR2S] + p.prog[R2][L2] * z[Y4DR2S];
	  d[Y5DR2R] = p.sec_exit * (z[Y5TR2R] + z[Y5UR2R]) - (auptake[L3] + p.prog[R2][L3] + p.revert[NONE][R2]) * z[Y5DR2R] + p.prog[R2][L2] * z[Y4DR2R];
	  d[Y5DC2S] = p.sec_exit * (z[Y5TC2S] + z[Y5UC2S]) - (auptake[L3] + p.prog[C2][L3] + p.revert[NONE][C2]) * z[Y5DC2S] + p.prog[C2][L2] * z[Y4DC2S];
	  d[Y5DC2R] = p.sec_exit * (z[Y5TC2R] + z[Y5UC2R]) - (auptake[L3] + p.prog[C2][L3] + p.revert[NONE][C2]) * z[Y5DC2R] + p.prog[C2][L2] * z[Y4DC2R];

	  d[Y5Dr2S] = p.sec_exit * (z[Y5Tr2S] + z[Y5Ur2S]) + p.revert[NONE][R2] * z[Y5DR2S] - (auptake[L3] + p.prog[WR2][L3]) * z[Y5Dr2S] + p.prog[WR2][L2] * z[Y4Dr2S];
	  d[Y5Dr2R] = p.sec_exit * (z[Y5Tr2R] + z[Y5Ur2R]) + p.revert[NONE][R2] * z[Y5DR2R] - (auptake[L3] + p.prog[WR2][L3]) * z[Y5Dr2R] + p.prog[WR2][L2] * z[Y4Dr2R];
	  d[Y5Dc2S] = p.sec_exit * (z[Y5Tc2S] + z[Y5Uc2S]) + p.revert[NONE][C2] * z[Y5DC2S] - (auptake[L3] + p.prog[WC2][L3]) * z[Y5Dc2S] + p.prog[WC2][L2] * z[Y4Dc2S];
	  d[Y5Dc2R] = p.sec_exit * (z[Y5Tc2R] + z[Y5Uc2R]) + p.revert[NONE][C2] * z[Y5DC2R] - (auptake[L3] + p.prog[WC2][L3]) * z[Y5Dc2R] + p.prog[WC2][L2] * z[Y4Dc2R];

 	  d[Y5UR2S] = p.sec_fail_dr * z[Y5TR2S] - (p.prog[R2][L3] + p.sec_exit + p.revert[ART2_NA][R2]) * z[Y5UR2S] + p.prog[R2][L2] * z[Y4UR2S]; // 2014-10-20 DONE
 	  d[Y5UR2R] = p.sec_fail_dr * z[Y5TR2R] - (p.prog[R2][L3] + p.sec_exit + p.revert[ART2_NA][R2]) * z[Y5UR2R] + p.prog[R2][L2] * z[Y4UR2R]; // 2014-10-20 DONE
 	  d[Y5UC2S] = p.sec_fail_dr * z[Y5TC2S] - (p.prog[C2][L3] + p.sec_exit + p.revert[ART2_NA][C2]) * z[Y5UC2S] + p.prog[C2][L2] * z[Y4UC2S]; // 2014-10-20 DONE
 	  d[Y5UC2R] = p.sec_fail_dr * z[Y5TC2R] - (p.prog[C2][L3] + p.sec_exit + p.revert[ART2_NA][C2]) * z[Y5UC2R] + p.prog[C2][L2] * z[Y4UC2R]; // 2014-10-20 DONE

	  d[Y5Ur2S] = p.sec_fail_dr * z[Y5Tr2S] + p.revert[ART2_NA][R2] * z[Y5UR2S] - (p.prog[WR2][L3] + p.sec_exit) * z[Y5Ur2S] + p.prog[WR2][L2] * z[Y4Ur2S]; // 2014-10-20 DONE
	  d[Y5Ur2R] = p.sec_fail_dr * z[Y5Tr2R] + p.revert[ART2_NA][R2] * z[Y5UR2R] - (p.prog[WR2][L3] + p.sec_exit) * z[Y5Ur2R] + p.prog[WR2][L2] * z[Y4Ur2R]; // 2014-10-20 DONE
	  d[Y5Uc2S] = p.sec_fail_dr * z[Y5Tc2S] + p.revert[ART2_NA][C2] * z[Y5UC2S] - (p.prog[WC2][L3] + p.sec_exit) * z[Y5Uc2S] + p.prog[WC2][L2] * z[Y4Uc2S]; // 2014-10-20 DONE
	  d[Y5Uc2R] = p.sec_fail_dr * z[Y5Tc2R] + p.revert[ART2_NA][C2] * z[Y5UC2R] - (p.prog[WC2][L3] + p.sec_exit) * z[Y5Uc2R] + p.prog[WC2][L2] * z[Y4Uc2R]; // 2014-10-20 DONE

	  // 2014-10-20 DONE
	  emergeS2 = p.emerge[ART2][WT] * z[Y5TWT]
	    + p.emerge[ART2][R1] * z[Y5TR1] + p.emerge[ART2][WR1] * z[Y5Tr1]
	    + p.emerge[ART2][C1] * z[Y5TC1] + p.emerge[ART2][WC1] * z[Y5Tc1]
	    + p.emerge[ART2][S1] * z[Y5TS1] + p.emerge[ART2][WS1] * z[Y5Ts1]
	    + p.emerge[ART2][Q1] * z[Y5TQ1] + p.emerge[ART2][WQ1] * z[Y5Tq1]
	    + p.emerge[ART2][Q2] * z[Y5TQ2] + p.emerge[ART2][WQ2] * z[Y5Tq2]
	    + p.emerge[ART2][ R2] * z[Y5TR2S] + p.emerge[ART2][ C2] * z[Y5TC2S]
	    + p.emerge[ART2][WR2] * z[Y5Tr2S] + p.emerge[ART2][WC2] * z[Y5Tc2S]
	    + p.emerge[ART2][WS1] * z[Y5TR2R] + p.emerge[ART2][WS1] * z[Y5TC2R]
	    + p.emerge[ART2][WS1] * z[Y5Tr2R] + p.emerge[ART2][WS1] * z[Y5Tc2R];

	  d[Y5TS2] = auptake[L3] * (z[Y5DS2] + z[Y5Ds2]) + p.prog[S2][L2] * z[Y4TS2] + emergeS2 - (p.prog[S2][L3] + p.sec_exit) * z[Y5TS2];
	  d[Y5DS2] = p.sec_exit * z[Y5TS2] - (auptake[L3] + p.prog[S2][L3] + p.revert[NONE][S2]) * z[Y5DS2] + p.prog[S2][L2] * z[Y4DS2];
	  d[Y5Ds2] = p.revert[NONE][S2] * z[Y5DS2] - (auptake[L3] + p.prog[WS2][L3]) * z[Y5Ds2] + p.prog[WS2][L2] * z[Y4Ds2];

	  // update cumulative HIV-related deaths with those occuring with CD4<350 (on ART; individuals off ART progress to later stages)
	  mort
	    = p.art_mortE[L3] * (z[Y5EWT] + z[Y5ER1] + z[Y5EC1] + z[Y5ES1] + z[Y5EQ1] + z[Y5EQ2] + z[Y5Er1] + z[Y5Ec1] + z[Y5Es1] + z[Y5Eq1] + z[Y5Eq2])
	    + p.art_mortL[L3] * (z[Y5LWT] + z[Y5LR1] + z[Y5LC1] + z[Y5LS1] + z[Y5LQ1] + z[Y5LQ2] + z[Y5Lr1] + z[Y5Lc1] + z[Y5Ls1] + z[Y5Lq1] + z[Y5Lq2])
	    + p.sec_mort[L3] * (z[Y5TR2S] + z[Y5TR2R] + z[Y5TC2S] + z[Y5TC2R] + z[Y5Tr2S] + z[Y5Tr2R] + z[Y5Tc2S] + z[Y5Tc2R])
	    + p.sec_mort[L3] * (z[Y5TWT] + z[Y5TR1] + z[Y5TC1] + z[Y5TS1] + z[Y5TQ1] + z[Y5TQ2] + z[Y5Tr1] + z[Y5Tc1] + z[Y5Ts1] + z[Y5Tq1] + z[Y5Tq2]); // 2014-10-20 DONE

	  if (b < BANDS - 1) {
	    dy[ACTIVE_HIV_MORT] += mort;
	    dy[ACTIVE_W_HIV_MORT] += mort * (g == F);
	    dy[ACTIVE_M_HIV_MORT] += mort * (g != F);

	    dy[ACTDSC_HIV_MORT] += discount * mort;
	    dy[ACTDSC_W_HIV_MORT] += discount * mort * (g == F);
	    dy[ACTDSC_M_HIV_MORT] += discount * mort * (g != F);
	  }

	  if (u == ANTE) {
	    dy[COHORT_HIV_MORT] += mort;
	    dy[COHORT_W_HIV_MORT] += mort * (g == F);
	    dy[COHORT_M_HIV_MORT] += mort * (g != F);

	    dy[COHDSC_HIV_MORT] += discount * mort;
	    dy[COHDSC_W_HIV_MORT] += discount * mort * (g == F);
	    dy[COHDSC_M_HIV_MORT] += discount * mort * (g != F);
	  }

	  // +-+ AIDS (CD4<200) +------------------------------------------------+
	  d[Y6NWT] = p.prog[WT][L3] * z[Y5NWT] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y6HWT] + z[Y6PWT]) + p.art_exitE * z[Y6EWT] + p.art_exitL * z[Y6LWT] + p.art_exitF * z[Y6FWT] - (p.prog[WT][A1] + auptake[A1]) * z[Y6NWT];
	  d[Y6HWT] = p.prog[WT][L3] * z[Y5HWT] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6PWT] - (p.prog[WT][A1] + p.emerge[PREP_HIGH][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y6HWT];
	  d[Y6PWT] = p.prog[WT][L3] * z[Y5PWT] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6HWT] - (p.prog[WT][A1] + p.emerge[PREP_POOR][WT] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y6PWT];
	  d[Y6EWT] = auptake[A1] * z[Y6NWT] - (p.emerge[ART1_NEW][WT] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6EWT];
	  d[Y6LWT] = p.art_advance * z[Y6EWT] - (p.emerge[ART1_OLD][WT] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6LWT];
	  d[Y6FWT] = p.prog[WT][L3] * z[Y5FWT] + p.art_failE * z[Y6EWT] + p.art_failL * z[Y6LWT] - (p.prog[WT][A1] + p.art_exitF + sec_init_na) * z[Y6FWT];
	  d[Y6TWT] = sec_init_na * z[Y6FWT] + auptake[A1] * z[Y6DWT] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WT]) * z[Y6TWT]; // 2014-10-20 DONE
	  d[Y6DWT] = p.prog[WT][L3] * z[Y5DWT] + p.sec_exit * (z[Y6TWT] + z[Y6UWT]) - (auptake[A1] + p.prog[WT][A1]) * z[Y6DWT]; // 2014-10-20 DONE
	  d[Y6UWT] = p.prog[WT][L3] * z[Y5UWT] + p.sec_fail_na * z[Y6TWT] - (p.prog[WT][A1] + p.sec_exit) * z[Y6UWT]; // 2014-10-20 DONE

	  d[Y6NR1] = p.prog[R1][L3] * z[Y5NR1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y6HR1] + z[Y6PR1]) + p.art_exitE * z[Y6ER1] + p.art_exitL * z[Y6LR1] + p.art_exitF * z[Y6FR1] - (p.prog[R1][A1] + p.revert[NONE][R1] + auptake[A1]) * z[Y6NR1];
	  d[Y6HR1] = p.prog[R1][L3] * z[Y5HR1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6PR1] - (p.prog[R1][A1] + p.revert[PREP_HIGH][R1] + p.emerge[PREP_HIGH][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y6HR1];
	  d[Y6PR1] = p.prog[R1][L3] * z[Y5PR1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6HR1] - (p.prog[R1][A1] + p.revert[PREP_POOR][R1] + p.emerge[PREP_POOR][R1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y6PR1];
	  d[Y6ER1] = auptake[A1] * z[Y6NR1] - (p.emerge[ART1_NEW][R1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6ER1];
	  d[Y6LR1] = p.art_advance * z[Y6ER1] - (p.emerge[ART1_OLD][R1] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6LR1];
	  d[Y6FR1] = p.prog[R1][L3] * z[Y5FR1] + p.art_failE * z[Y6ER1] + p.art_failL * z[Y6LR1] - (p.prog[R1][A1] + p.revert[ART1_NA][R1] + p.art_exitF + sec_init_na) * z[Y6FR1];
	  d[Y6TR1] = sec_init_na * z[Y6FR1] + auptake[A1] * z[Y6DR1] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.revert[ART2][R1] + p.emerge[ART2][R1]) * z[Y6TR1];
	  d[Y6DR1] = p.prog[R1][L3] * z[Y5DR1] + p.sec_exit * (z[Y6TR1] + z[Y6UR1]) - (auptake[A1] + p.prog[R1][A1] + p.revert[NONE][R1]) * z[Y6DR1]; // 2014-10-20 DONE
	  d[Y6UR1] = p.prog[R1][L3] * z[Y5UR1] + p.sec_fail_na * z[Y6TR1] - (p.prog[R1][A1] + p.revert[ART2_NA][R1] + p.sec_exit) * z[Y6UR1]; // 2014-10-20 DONE

	  d[Y6NC1] = p.prog[C1][L3] * z[Y5NC1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y6HC1] + z[Y6PC1]) + p.art_exitE * z[Y6EC1] + p.art_exitL * z[Y6LC1] + p.art_exitF * z[Y6FC1] - (p.prog[C1][A1] + p.revert[NONE][C1] + auptake[A1]) * z[Y6NC1];
	  d[Y6HC1] = p.prog[C1][L3] * z[Y5HC1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6PC1] - (p.prog[C1][A1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y6HC1];
	  d[Y6PC1] = p.prog[C1][L3] * z[Y5PC1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6HC1] - (p.prog[C1][A1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y6PC1];
	  d[Y6EC1] = auptake[A1] * z[Y6NC1] - (p.emerge[ART1_NEW][C1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6EC1];
	  d[Y6LC1] = p.art_advance * z[Y6EC1] - (p.emerge[ART1_OLD][C1] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6LC1];
	  d[Y6FC1] = p.prog[C1][L3] * z[Y5FC1] + p.art_failE * z[Y6EC1] + p.art_failL * z[Y6LC1] - (p.prog[C1][A1] + p.revert[ART1_NA][C1] + p.art_exitF + sec_init_na) * z[Y6FC1];
	  d[Y6TC1] = sec_init_na * z[Y6FC1] + auptake[A1] * z[Y6DC1] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.revert[ART2][C1] + p.emerge[ART2][C1]) * z[Y6TC1];
	  d[Y6DC1] = p.prog[C1][L3] * z[Y5DC1] + p.sec_exit * (z[Y6TC1] + z[Y6UC1]) - (auptake[A1] + p.prog[C1][A1] + p.revert[NONE][C1]) * z[Y6DC1]; // 2014-10-20 DONE
	  d[Y6UC1] = p.prog[C1][L3] * z[Y5UC1] + p.sec_fail_na * z[Y6TC1] - (p.prog[C1][A1] + p.revert[ART2_NA][C1] + p.sec_exit) * z[Y6UC1]; // 2014-10-20 DONE

	  d[Y6NQ1] = p.prog[Q1][L3] * z[Y5NQ1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y6HQ1] + z[Y6PQ1]) + p.art_exitE * z[Y6EQ1] + p.art_exitL * z[Y6LQ1] + p.art_exitF * z[Y6FQ1] - (p.prog[Q1][A1] + p.revert[NONE][Q1] + auptake[A1]) * z[Y6NQ1];
	  d[Y6HQ1] = p.prog[Q1][L3] * z[Y5HQ1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6PQ1] - (p.prog[Q1][A1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y6HQ1];
	  d[Y6PQ1] = p.prog[Q1][L3] * z[Y5PQ1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6HQ1] - (p.prog[Q1][A1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y6PQ1];
	  d[Y6EQ1] = auptake[A1] * z[Y6NQ1] - (p.emerge[ART1_NEW][Q1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6EQ1];
	  d[Y6LQ1] = p.art_advance * z[Y6EQ1] - (p.emerge[ART1_OLD][Q1] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6LQ1];
	  d[Y6FQ1] = p.prog[Q1][L3] * z[Y5FQ1] + p.art_failE * z[Y6EQ1] + p.art_failL * z[Y6LQ1] - (p.prog[Q1][A1] + p.revert[ART1_NA][Q1] + p.art_exitF + sec_init_na) * z[Y6FQ1];
	  d[Y6TQ1] = sec_init_na * z[Y6FQ1] + auptake[A1] * z[Y6DQ1] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.revert[ART2][Q1] + p.emerge[ART2][Q1]) * z[Y6TQ1];
	  d[Y6DQ1] = p.prog[Q1][L3] * z[Y5DQ1] + p.sec_exit * (z[Y6TQ1] + z[Y6UQ1]) - (auptake[A1] + p.prog[Q1][A1] + p.revert[NONE][Q1]) * z[Y6DQ1]; // 2014-10-20 DONE
	  d[Y6UQ1] = p.prog[Q1][L3] * z[Y5UQ1] + p.sec_fail_na * z[Y6TQ1] - (p.prog[Q1][A1] + p.revert[ART2_NA][Q1] + p.sec_exit) * z[Y6UQ1]; // 2014-10-20 DONE

	  d[Y6NS1] = p.prog[S1][L3] * z[Y5NS1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y6HS1] + z[Y6PS1]) + p.art_exitE * z[Y6ES1] + p.art_exitL * z[Y6LS1] + p.art_exitF * z[Y6FS1] - (p.prog[S1][A1] + p.revert[NONE][S1] + auptake[A1]) * z[Y6NS1];
	  d[Y6HS1] = p.prog[S1][L3] * z[Y5HS1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6PS1] - (p.prog[S1][A1] + p.revert[PREP_HIGH][S1] + p.emerge[PREP_HIGH][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y6HS1];
	  d[Y6PS1] = p.prog[S1][L3] * z[Y5PS1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6HS1] - (p.prog[S1][A1] + p.revert[PREP_POOR][S1] + p.emerge[PREP_POOR][S1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y6PS1];
	  d[Y6ES1] = auptake[A1] * z[Y6NS1] - (p.emerge[ART1_NEW][S1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6ES1];
	  d[Y6LS1] = p.art_advance * z[Y6ES1] - (p.emerge[ART1_OLD][S1] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6LS1];
	  d[Y6FS1] = p.prog[S1][L3] * z[Y5FS1] + p.art_failE * z[Y6ES1] + p.art_failL * z[Y6LS1] - (p.prog[S1][A1] + p.revert[ART1_NA][S1] + p.art_exitF + sec_init_na) * z[Y6FS1];
	  d[Y6TS1] = sec_init_na * z[Y6FS1] + auptake[A1] * z[Y6DS1] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][S1]) * z[Y6TS1]; // 2014-10-20 DONE
	  d[Y6DS1] = p.prog[S1][L3] * z[Y5DS1] + p.sec_exit * (z[Y6TS1] + z[Y6US1]) - (auptake[A1] + p.prog[S1][A1] + p.revert[NONE][S1]) * z[Y6DS1]; // 2014-10-20 DONE
	  d[Y6US1] = p.prog[S1][L3] * z[Y5US1] + p.sec_fail_na * z[Y6TS1] - (p.prog[S1][A1] + p.revert[ART2_NA][S1] + p.sec_exit) * z[Y6US1]; // 2014-10-20 DONE

	  d[Y6NQ2] = p.prog[Q2][L3] * z[Y5NQ2] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y6HQ2] + z[Y6PQ2]) + p.art_exitE * z[Y6EQ2] + p.art_exitL * z[Y6LQ2] + p.art_exitF * z[Y6FQ2] - (p.prog[Q2][A1] + p.revert[NONE][Q2] + auptake[A1]) * z[Y6NQ2];
	  d[Y6HQ2] = p.prog[Q2][L3] * z[Y5HQ2] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6PQ2] + p.emerge[PREP_HIGH][WT] * z[Y6HWT] + p.emerge[PREP_HIGH][R1] * z[Y6HR1] + p.emerge[PREP_HIGH][WR1] * z[Y6Hr1] + p.emerge[PREP_HIGH][S1] * z[Y6HS1] + p.emerge[PREP_HIGH][WS1] * z[Y6Hs1] + p.emerge[PREP_HIGH][WC1] * z[Y6Hc1] + p.emerge[PREP_HIGH][WQ1] * z[Y6Hq1] - (p.prog[Q2][A1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y6HQ2];
	  d[Y6PQ2] = p.prog[Q2][L3] * z[Y5PQ2] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6HQ2] + p.emerge[PREP_POOR][WT] * z[Y6PWT] + p.emerge[PREP_POOR][R1] * z[Y6PR1] + p.emerge[PREP_POOR][WR1] * z[Y6Pr1] + p.emerge[PREP_POOR][S1] * z[Y6PS1] + p.emerge[PREP_POOR][WS1] * z[Y6Ps1] + p.emerge[PREP_POOR][WC1] * z[Y6Pc1] + p.emerge[PREP_POOR][WQ1] * z[Y6Pq1] - (p.prog[Q2][A1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y6PQ2];
	  d[Y6EQ2] = auptake[A1] * z[Y6NQ2] - (p.emerge[ART1_NEW][Q2] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6EQ2];
	  d[Y6LQ2] = p.art_advance * z[Y6EQ2] - (p.emerge[ART1_OLD][Q2] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6LQ2];
	  d[Y6FQ2] = p.prog[Q2][L3] * z[Y5FQ2] + p.art_failE * z[Y6EQ2] + p.art_failL * z[Y6LQ2] - (p.prog[Q2][A1] + p.revert[ART1_NA][Q2] + p.art_exitF + sec_init_na) * z[Y6FQ2];
	  d[Y6TQ2] = sec_init_na * z[Y6FQ2] + auptake[A1] * z[Y6DQ2] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.revert[ART2][Q2] + p.emerge[ART2][Q2]) * z[Y6TQ2];
	  d[Y6DQ2] = p.prog[Q2][L3] * z[Y5DQ2] + p.sec_exit * (z[Y6TQ2] + z[Y6UQ2]) - (auptake[A1] + p.prog[Q2][A1] + p.revert[NONE][Q2]) * z[Y6DQ2]; // 2014-10-20 DONE
	  d[Y6UQ2] = p.prog[Q2][L3] * z[Y5UQ2] + p.sec_fail_na * z[Y6TQ2] - (p.prog[Q2][A1] + p.revert[ART2_NA][Q2] + p.sec_exit) * z[Y6UQ2]; // 2014-10-20 DONE

	  emergeR2S = (1.0 - p.xdr_overall) *
	    (p.emerge[ART1_NEW][WT ] * z[Y6EWT] + p.emerge[ART1_OLD][WT ] * z[Y6LWT] +
	     p.emerge[ART1_NEW][R1 ] * z[Y6ER1] + p.emerge[ART1_OLD][R1 ] * z[Y6LR1] +
	     p.emerge[ART1_NEW][WR1] * z[Y6Er1] + p.emerge[ART1_OLD][WR1] * z[Y6Lr1]);
	  emergeR2R = (1.0 - p.xdr_overall) *
	    (p.emerge[ART1_NEW][S1 ] * z[Y6ES1] + p.emerge[ART1_OLD][S1 ] * z[Y6LS1] +
	     p.emerge[ART1_NEW][WS1] * z[Y6Es1] + p.emerge[ART1_OLD][WS1] * z[Y6Ls1]);

	  emergeC2S
	    = p.xdr_overall
	    * (p.emerge[ART1_NEW][WT ] * z[Y6EWT] + p.emerge[ART1_OLD][WT ] * z[Y6LWT] +
	       p.emerge[ART1_NEW][R1 ] * z[Y6ER1] + p.emerge[ART1_OLD][R1 ] * z[Y6LR1] +
	       p.emerge[ART1_NEW][WR1] * z[Y6Er1] + p.emerge[ART1_OLD][WR1] * z[Y6Lr1])
	    + (p.emerge[ART1_NEW][C1 ] * z[Y6EC1] + p.emerge[ART1_OLD][C1 ] * z[Y6LC1] +
	       p.emerge[ART1_NEW][WC1] * z[Y6Ec1] + p.emerge[ART1_OLD][WC1] * z[Y6Lc1] +
	       p.emerge[ART1_NEW][Q1 ] * z[Y6EQ1] + p.emerge[ART1_OLD][Q1 ] * z[Y6LQ1] +
	       p.emerge[ART1_NEW][WQ1] * z[Y6Eq1] + p.emerge[ART1_OLD][WQ1] * z[Y6Lq1] +
	       p.emerge[ART1_NEW][Q2 ] * z[Y6EQ2] + p.emerge[ART1_OLD][Q2 ] * z[Y6LQ2] +
	       p.emerge[ART1_NEW][WQ2] * z[Y6Eq2] + p.emerge[ART1_OLD][WQ2] * z[Y6Lq2]);
	  emergeC2R
	    = p.xdr_overall
	    * (p.emerge[ART1_NEW][S1 ] * z[Y6ES1] + p.emerge[ART1_OLD][S1 ] * z[Y6LS1] +
	       p.emerge[ART1_NEW][WS1] * z[Y6Es1] + p.emerge[ART1_OLD][WS1] * z[Y6Ls1]);

	  d[Y6LR2S] = p.prog[ R2][L3] * z[Y5LR2S] + auptake[A1] * (z[Y6NR2S] + z[Y6Nr2S]) + emergeR2S - (p.prog[R2][A1] + p.art_exitL + sec_init_dr) * z[Y6LR2S];
	  d[Y6NR2S] = p.prog[ R2][L3] * z[Y5NR2S] + p.art_exitL * z[Y6LR2S] - (p.prog[R2][A1] + p.revert[NONE][R2] + auptake[A1]) * z[Y6NR2S];
	  d[Y6Nr2S] = p.prog[WR2][L3] * z[Y5Nr2S] + p.revert[NONE][R2] * z[Y6NR2S] - (p.prog[WR2][A1] + auptake[A1]) * z[Y6Nr2S];
	  d[Y6LR2R] = p.prog[ R2][L3] * z[Y5LR2R] + auptake[A1] * (z[Y6NR2R] + z[Y6Nr2R]) + emergeR2R - (p.prog[R2][A1] + p.art_exitL + sec_init_dr) * z[Y6LR2R];
	  d[Y6NR2R] = p.prog[ R2][L3] * z[Y5NR2R] + p.art_exitL * z[Y6LR2R] - (p.prog[R2][A1] + p.revert[NONE][R2] + auptake[A1]) * z[Y6NR2R];
	  d[Y6Nr2R] = p.prog[WR2][L3] * z[Y5Nr2R] + p.revert[NONE][R2] * z[Y6NR2R] - (p.prog[WR2][A1] + auptake[A1]) * z[Y6Nr2R];

	  d[Y6LC2S] = p.prog[ C2][L3] * z[Y5LC2S] + auptake[A1] * (z[Y6NC2S] + z[Y6Nc2S]) + emergeC2S - (p.prog[C2][A1] + p.art_exitL + sec_init_dr) * z[Y6LC2S];
	  d[Y6NC2S] = p.prog[ C2][L3] * z[Y5NC2S] + p.art_exitL * z[Y6LC2S] - (p.prog[C2][A1] + p.revert[NONE][C2] + auptake[A1]) * z[Y6NC2S];
	  d[Y6Nc2S] = p.prog[WC2][L3] * z[Y5Nc2S] + p.revert[NONE][C2] * z[Y6NC2S] - (p.prog[WC2][A1] + auptake[A1]) * z[Y6Nc2S];
	  d[Y6LC2R] = p.prog[ C2][L3] * z[Y5LC2R] + auptake[A1] * (z[Y6NC2R] + z[Y6Nc2R]) + emergeC2R - (p.prog[C2][A1] + p.art_exitL + sec_init_dr) * z[Y6LC2R];
	  d[Y6NC2R] = p.prog[ C2][L3] * z[Y5NC2R] + p.art_exitL * z[Y6LC2R] - (p.prog[C2][A1] + p.revert[NONE][C2] + auptake[A1]) * z[Y6NC2R];
	  d[Y6Nc2R] = p.prog[WC2][L3] * z[Y5Nc2R] + p.revert[NONE][C2] * z[Y6NC2R] - (p.prog[WC2][A1] + auptake[A1]) * z[Y6Nc2R];

	  d[Y6Nr1] = p.prog[WR1][L3] * z[Y5Nr1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y6Hr1] + z[Y6Pr1]) + p.revert[NONE][R1] * z[Y6NR1] + p.art_exitE * z[Y6Er1] + p.art_exitL * z[Y6Lr1] + p.art_exitF * z[Y6Fr1] - (p.prog[WR1][A1] + auptake[A1]) * z[Y6Nr1];
	  d[Y6Hr1] = p.prog[WR1][L3] * z[Y5Hr1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6Pr1] + p.revert[PREP_HIGH][R1] * z[Y6HR1] - (p.prog[WR1][A1] + p.emerge[PREP_HIGH][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y6Hr1];
	  d[Y6Pr1] = p.prog[WR1][L3] * z[Y5Pr1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6Hr1] + p.revert[PREP_POOR][R1] * z[Y6PR1] - (p.prog[WR1][A1] + p.emerge[PREP_POOR][WR1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y6Pr1];
	  d[Y6Er1] = auptake[A1] * z[Y6Nr1] - (p.emerge[ART1_NEW][WR1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6Er1];
	  d[Y6Lr1] = p.art_advance * z[Y6Er1] - (p.emerge[ART1_OLD][WR1] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6Lr1];
	  d[Y6Fr1] = p.revert[ART1_NA][R1] * z[Y6FR1] + p.prog[WR1][L3] * z[Y5Fr1] + p.art_failE * z[Y6Er1] + p.art_failL * z[Y6Lr1] - (p.prog[WR1][A1] + p.art_exitF + sec_init_na) * z[Y6Fr1];
	  d[Y6Tr1] = sec_init_na * z[Y6Fr1] + auptake[A1] * z[Y6Dr1] + p.revert[ART2][R1] * z[Y6TR1] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WR1]) * z[Y6Tr1];
	  d[Y6Dr1] = p.prog[WR1][L3] * z[Y5Dr1] + p.sec_exit * (z[Y6Tr1] + z[Y6Ur1]) + p.revert[NONE][R1] * z[Y6DR1] - (auptake[A1] + p.prog[WR1][A1]) * z[Y6Dr1]; // 2014-10-20 DONE
	  d[Y6Ur1] = p.prog[WR1][L3] * z[Y5Ur1] + p.sec_fail_na * z[Y6Tr1] + p.revert[ART2_NA][R1] * z[Y6UR1] - (p.prog[WR1][A1] + p.sec_exit) * z[Y6Ur1]; // 2014-10-20 DONE

	  d[Y6Nc1] = p.prog[WC1][L3] * z[Y5Nc1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k]) * (z[Y6Hc1] + z[Y6Pc1]) + p.revert[NONE][C1] * z[Y6NC1] + p.art_exitE * z[Y6Ec1] + p.art_exitL * z[Y6Lc1] + p.art_exitF * z[Y6Fc1] - (p.prog[WC1][A1] + auptake[A1]) * z[Y6Nc1];
	  d[Y6Hc1] = p.prog[WC1][L3] * z[Y5Hc1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6Pc1] - (p.prog[WC1][A1] + p.emerge[PREP_HIGH][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y6Hc1];
	  d[Y6Pc1] = p.prog[WC1][L3] * z[Y5Pc1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6Hc1] - (p.prog[WC1][A1] + p.emerge[PREP_POOR][WC1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y6Pc1];
	  d[Y6Ec1] = auptake[A1] * z[Y6Nc1] - (p.emerge[ART1_NEW][WC1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6Ec1];
	  d[Y6Lc1] = p.art_advance * z[Y6Ec1] - (p.emerge[ART1_OLD][WC1] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6Lc1];
	  d[Y6Fc1] = p.revert[ART1_NA][C1] * z[Y6FC1] + p.prog[WC1][L3] * z[Y5Fc1] + p.art_failE * z[Y6Ec1] + p.art_failL * z[Y6Lc1] - (p.prog[WC1][A1] + p.art_exitF + sec_init_na) * z[Y6Fc1];
	  d[Y6Tc1] = sec_init_na * z[Y6Fc1] + auptake[A1] * z[Y6Dc1] + p.revert[ART2][C1] * z[Y6TC1] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WC1]) * z[Y6Tc1];
	  d[Y6Dc1] = p.prog[WC1][L3] * z[Y5Dc1] + p.sec_exit * (z[Y6Tc1] + z[Y6Uc1]) + p.revert[NONE][C1] * z[Y6DC1] - (auptake[A1] + p.prog[WC1][A1]) * z[Y6Dc1]; // 2014-10-20 DONE
	  d[Y6Uc1] = p.prog[WC1][L3] * z[Y5Uc1] + p.sec_fail_na * z[Y6Tc1] + p.revert[ART2_NA][C1] * z[Y6UC1] - (p.prog[WC1][A1] + p.sec_exit) * z[Y6Uc1]; // 2014-10-20 DONE

	  d[Y6Nq1] = p.prog[WQ1][L3] * z[Y5Nq1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k])  * (z[Y6Hq1] + z[Y6Pq1]) + p.revert[NONE][Q1] * z[Y6NQ1] + p.art_exitE * z[Y6Eq1] + p.art_exitL * z[Y6Lq1] + p.art_exitF * z[Y6Fq1] - (p.prog[WQ1][A1] + auptake[A1]) * z[Y6Nq1];
	  d[Y6Hq1] = p.prog[WQ1][L3] * z[Y5Hq1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6Pq1] - (p.prog[WQ1][A1] + p.emerge[PREP_HIGH][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y6Hq1];
	  d[Y6Pq1] = p.prog[WQ1][L3] * z[Y5Pq1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6Hq1] - (p.prog[WQ1][A1] + p.emerge[PREP_POOR][WQ1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y6Pq1];
	  d[Y6Eq1] = auptake[A1] * z[Y6Nq1] - (p.emerge[ART1_NEW][WQ1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6Eq1];
	  d[Y6Lq1] = p.art_advance * z[Y6Eq1] - (p.emerge[ART1_OLD][WQ1] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6Lq1];
	  d[Y6Fq1] = p.revert[ART1_NA][Q1] * z[Y6FQ1] + p.prog[WQ1][L3] * z[Y5Fq1] + p.art_failE * z[Y6Eq1] + p.art_failL * z[Y6Lq1] - (p.prog[WQ1][A1] + p.art_exitF + sec_init_na) * z[Y6Fq1];
	  d[Y6Tq1] = sec_init_na * z[Y6Fq1] + auptake[A1] * z[Y6Dq1] + p.revert[ART2][Q1] * z[Y6TQ1] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WQ1]) * z[Y6Tq1];
	  d[Y6Dq1] = p.prog[WQ1][L3] * z[Y5Dq1] + p.sec_exit * (z[Y6Tq1] + z[Y6Uq1]) + p.revert[NONE][Q1] * z[Y6DQ1] - (auptake[A1] + p.prog[WQ1][A1]) * z[Y6Dq1]; // 2014-10-20 DONE
	  d[Y6Uq1] = p.prog[WQ1][L3] * z[Y5Uq1] + p.sec_fail_na * z[Y6Tq1] + p.revert[ART2_NA][Q1] * z[Y6UQ1] - (p.prog[WQ1][A1] + p.sec_exit) * z[Y6Uq1]; // 2014-10-20 DONE

	  d[Y6Ns1] = p.prog[WS1][L3] * z[Y5Ns1] + (prep_exit[sex[g]][b][k] + p.prep_test * prep_inject[sex[g]][b][k])  * (z[Y6Hs1] + z[Y6Ps1]) + p.revert[NONE][S1] * z[Y6NS1] + p.art_exitE * z[Y6Es1] + p.art_exitL * z[Y6Ls1] + p.art_exitF * z[Y6Fs1] - (p.prog[WS1][A1] + auptake[A1]) * z[Y6Ns1];
	  d[Y6Hs1] = p.prog[WS1][L3] * z[Y5Hs1] + prep_propH * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6Ps1] + p.revert[PREP_HIGH][S1] * z[Y6HS1] - (p.prog[WS1][A1] + p.emerge[PREP_HIGH][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propP) + prep_exit[sex[g]][b][k]) * z[Y6Hs1];
	  d[Y6Ps1] = p.prog[WS1][L3] * z[Y5Ps1] + prep_propP * (1 - p.prep_test) * prep_inject[sex[g]][b][k] * z[Y6Hs1] + p.revert[PREP_POOR][S1] * z[Y6PS1] - (p.prog[WS1][A1] + p.emerge[PREP_POOR][WS1] + prep_inject[sex[g]][b][k] * (p.prep_test + (1 - p.prep_test) * prep_propH) + prep_exit[sex[g]][b][k]) * z[Y6Ps1];
	  d[Y6Es1] = auptake[A1] * z[Y6Ns1] - (p.emerge[ART1_NEW][WS1] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6Es1];
	  d[Y6Ls1] = p.art_advance * z[Y6Es1] - (p.emerge[ART1_OLD][WS1] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6Ls1];
	  d[Y6Fs1] = p.revert[ART1_NA][S1] * z[Y6FS1] + p.prog[WS1][L3] * z[Y5Fs1] + p.art_failE * z[Y6Es1] + p.art_failL * z[Y6Ls1] - (p.prog[WS1][A1] + p.art_exitF + sec_init_na) * z[Y6Fs1];
	  d[Y6Ts1] = sec_init_na * z[Y6Fs1] + auptake[A1] * z[Y6Ds1] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WS1]) * z[Y6Ts1]; // 2014-10-20 DONE
	  d[Y6Ds1] = p.prog[WS1][L3] * z[Y5Ds1] + p.sec_exit * (z[Y6Ts1] + z[Y6Us1]) + p.revert[NONE][S1] * z[Y6DS1] - (auptake[A1] + p.prog[WS1][A1]) * z[Y6Ds1]; // 2014-10-20 DONE
	  d[Y6Us1] = p.prog[WS1][L3] * z[Y5Us1] + p.sec_fail_na * z[Y6Ts1] + p.revert[ART2_NA][S1] * z[Y6US1] - (p.prog[WS1][A1] + p.sec_exit) * z[Y6Us1]; // 2014-10-20 DONE

	  d[Y6Nq2] = p.prog[WQ2][L3] * z[Y5Nq2] + p.revert[NONE][Q2] * z[Y6NQ2] + p.art_exitE * z[Y6Eq2] + p.art_exitL * z[Y6Lq2] + p.art_exitF * z[Y6Fq2] - (p.prog[WQ2][A1] + auptake[A1]) * z[Y6Nq2];
	  d[Y6Eq2] = auptake[A1] * z[Y6Nq2] - (p.emerge[ART1_NEW][WQ2] + p.art_failE + p.art_exitE + p.art_advance + p.art_mortE[A1]) * z[Y6Eq2];
	  d[Y6Lq2] = p.art_advance * z[Y6Eq2] - (p.emerge[ART1_OLD][WQ2] + p.art_failL + p.art_exitL + p.art_mortL[A1]) * z[Y6Lq2];
	  d[Y6Fq2] = p.revert[ART1_NA][Q2] * z[Y6FQ2] + p.prog[WQ2][L3] * z[Y5Fq2] + p.art_failE * z[Y6Eq2] + p.art_failL * z[Y6Lq2] - (p.prog[WQ2][A1] + p.art_exitF + sec_init_na) * z[Y6Fq2];
	  d[Y6Tq2] = sec_init_na * z[Y6Fq2] + auptake[A1] * z[Y6Dq2] + p.revert[ART2][Q2] * z[Y6TQ2] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_na + p.emerge[ART2][WQ2]) * z[Y6Tq2];
	  d[Y6Dq2] = p.prog[WQ2][L3] * z[Y5Dq2] + p.sec_exit * (z[Y6Tq2] + z[Y6Uq2]) + p.revert[NONE][Q2] * z[Y6DQ2] - (auptake[A1] + p.prog[WQ2][A1]) * z[Y6Dq2]; // 2014-10-20 DONE
	  d[Y6Uq2] = p.prog[WQ2][L3] * z[Y5Uq2] + p.sec_fail_na * z[Y6Tq2] + p.revert[ART2_NA][Q2] * z[Y6UQ2] - (p.prog[WQ2][A1] + p.sec_exit) * z[Y6Uq2]; // 2014-10-20 DONE

	  d[Y6TR2S] = sec_init_dr * z[Y6LR2S] + auptake[A1] * z[Y6DR2S] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][R2 ] + p.revert[ART2][R2]) * z[Y6TR2S];
	  d[Y6TR2R] = sec_init_dr * z[Y6LR2R] + auptake[A1] * z[Y6DR2R] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][WS1] + p.revert[ART2][R2]) * z[Y6TR2R];
	  d[Y6TC2S] = sec_init_dr * z[Y6LC2S] + auptake[A1] * z[Y6DC2S] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][C2 ] + p.revert[ART2][C2]) * z[Y6TC2S];
	  d[Y6TC2R] = sec_init_dr * z[Y6LC2R] + auptake[A1] * z[Y6DC2R] - (p.sec_mort[A1] + p.sec_exit + p.sec_fail_dr + p.emerge[ART2][WS1] + p.revert[ART2][C2]) * z[Y6TC2R];

	  d[Y6Tr2S] = auptake[A1] * z[Y6Dr2S] + p.revert[ART2][R2] * z[Y6TR2S] - (p.sec_mort[A1] + p.sec_fail_dr + p.emerge[ART2][WR2] + p.sec_exit) * z[Y6Tr2S];
	  d[Y6Tr2R] = auptake[A1] * z[Y6Dr2R] + p.revert[ART2][R2] * z[Y6TR2R] - (p.sec_mort[A1] + p.sec_fail_dr + p.emerge[ART2][WS1] + p.sec_exit) * z[Y6Tr2R];
	  d[Y6Tc2S] = auptake[A1] * z[Y6Dc2S] + p.revert[ART2][C2] * z[Y6TC2S] - (p.sec_mort[A1] + p.sec_fail_dr + p.emerge[ART2][WC2] + p.sec_exit) * z[Y6Tc2S];
	  d[Y6Tc2R] = auptake[A1] * z[Y6Dc2R] + p.revert[ART2][C2] * z[Y6TC2R] - (p.sec_mort[A1] + p.sec_fail_dr + p.emerge[ART2][WS1] + p.sec_exit) * z[Y6Tc2R];

	  d[Y6DR2S] = p.sec_exit * (z[Y6TR2S] + z[Y6UR2S]) - (auptake[A1] + p.prog[R2][A1] + p.revert[NONE][R2]) * z[Y6DR2S] + p.prog[R2][L3] * z[Y5DR2S];
	  d[Y6DR2R] = p.sec_exit * (z[Y6TR2R] + z[Y6UR2R]) - (auptake[A1] + p.prog[R2][A1] + p.revert[NONE][R2]) * z[Y6DR2R] + p.prog[R2][L3] * z[Y5DR2R];
	  d[Y6DC2S] = p.sec_exit * (z[Y6TC2S] + z[Y6UC2S]) - (auptake[A1] + p.prog[C2][A1] + p.revert[NONE][C2]) * z[Y6DC2S] + p.prog[C2][L3] * z[Y5DC2S];
	  d[Y6DC2R] = p.sec_exit * (z[Y6TC2R] + z[Y6UC2R]) - (auptake[A1] + p.prog[C2][A1] + p.revert[NONE][C2]) * z[Y6DC2R] + p.prog[C2][L3] * z[Y5DC2R];

	  d[Y6Dr2S] = p.sec_exit * (z[Y6Tr2S] + z[Y6Ur2S]) + p.revert[NONE][R2] * z[Y6DR2S] - (auptake[A1] + p.prog[WR2][A1]) * z[Y6Dr2S] + p.prog[WR2][L3] * z[Y5Dr2S];
	  d[Y6Dr2R] = p.sec_exit * (z[Y6Tr2R] + z[Y6Ur2R]) + p.revert[NONE][R2] * z[Y6DR2R] - (auptake[A1] + p.prog[WR2][A1]) * z[Y6Dr2R] + p.prog[WR2][L3] * z[Y5Dr2R];
	  d[Y6Dc2S] = p.sec_exit * (z[Y6Tc2S] + z[Y6Uc2S]) + p.revert[NONE][C2] * z[Y6DC2S] - (auptake[A1] + p.prog[WC2][A1]) * z[Y6Dc2S] + p.prog[WC2][L3] * z[Y5Dc2S];
	  d[Y6Dc2R] = p.sec_exit * (z[Y6Tc2R] + z[Y6Uc2R]) + p.revert[NONE][C2] * z[Y6DC2R] - (auptake[A1] + p.prog[WC2][A1]) * z[Y6Dc2R] + p.prog[WC2][L3] * z[Y5Dc2R];

 	  d[Y6UR2S] = p.sec_fail_dr * z[Y6TR2S] - (p.prog[R2][A1] + p.sec_exit + p.revert[ART2_NA][R2]) * z[Y6UR2S] + p.prog[R2][L3] * z[Y5UR2S]; // 2014-10-20 DONE
 	  d[Y6UR2R] = p.sec_fail_dr * z[Y6TR2R] - (p.prog[R2][A1] + p.sec_exit + p.revert[ART2_NA][R2]) * z[Y6UR2R] + p.prog[R2][L3] * z[Y5UR2R]; // 2014-10-20 DONE
 	  d[Y6UC2S] = p.sec_fail_dr * z[Y6TC2S] - (p.prog[C2][A1] + p.sec_exit + p.revert[ART2_NA][C2]) * z[Y6UC2S] + p.prog[C2][L3] * z[Y5UC2S]; // 2014-10-20 DONE
 	  d[Y6UC2R] = p.sec_fail_dr * z[Y6TC2R] - (p.prog[C2][A1] + p.sec_exit + p.revert[ART2_NA][C2]) * z[Y6UC2R] + p.prog[C2][L3] * z[Y5UC2R]; // 2014-10-20 DONE

	  d[Y6Ur2S] = p.sec_fail_dr * z[Y6Tr2S] + p.revert[ART2_NA][R2] * z[Y6UR2S] - (p.prog[WR2][A1] + p.sec_exit) * z[Y6Ur2S] + p.prog[WR2][L3] * z[Y5Ur2S]; // 2014-10-20 DONE
	  d[Y6Ur2R] = p.sec_fail_dr * z[Y6Tr2R] + p.revert[ART2_NA][R2] * z[Y6UR2R] - (p.prog[WR2][A1] + p.sec_exit) * z[Y6Ur2R] + p.prog[WR2][L3] * z[Y5Ur2R]; // 2014-10-20 DONE
	  d[Y6Uc2S] = p.sec_fail_dr * z[Y6Tc2S] + p.revert[ART2_NA][C2] * z[Y6UC2S] - (p.prog[WC2][A1] + p.sec_exit) * z[Y6Uc2S] + p.prog[WC2][L3] * z[Y5Uc2S]; // 2014-10-20 DONE
	  d[Y6Uc2R] = p.sec_fail_dr * z[Y6Tc2R] + p.revert[ART2_NA][C2] * z[Y6UC2R] - (p.prog[WC2][A1] + p.sec_exit) * z[Y6Uc2R] + p.prog[WC2][L3] * z[Y5Uc2R]; // 2014-10-20 DONE

	  // 2014-10-20 DONE
	  emergeS2 = p.emerge[ART2][WT] * z[Y6TWT]
	    + p.emerge[ART2][R1] * z[Y6TR1] + p.emerge[ART2][WR1] * z[Y6Tr1]
	    + p.emerge[ART2][C1] * z[Y6TC1] + p.emerge[ART2][WC1] * z[Y6Tc1]
	    + p.emerge[ART2][S1] * z[Y6TS1] + p.emerge[ART2][WS1] * z[Y6Ts1]
	    + p.emerge[ART2][Q1] * z[Y6TQ1] + p.emerge[ART2][WQ1] * z[Y6Tq1]
	    + p.emerge[ART2][Q2] * z[Y6TQ2] + p.emerge[ART2][WQ2] * z[Y6Tq2]
	    + p.emerge[ART2][ R2] * z[Y6TR2S] + p.emerge[ART2][ C2] * z[Y6TC2S]
	    + p.emerge[ART2][WR2] * z[Y6Tr2S] + p.emerge[ART2][WC2] * z[Y6Tc2S]
	    + p.emerge[ART2][WS1] * z[Y6TR2R] + p.emerge[ART2][WS1] * z[Y6TC2R]
	    + p.emerge[ART2][WS1] * z[Y6Tr2R] + p.emerge[ART2][WS1] * z[Y6Tc2R];

	  d[Y6TS2] = auptake[A1] * (z[Y6DS2] + z[Y6Ds2]) + p.prog[S2][L3] * z[Y5TS2] + emergeS2 - (p.prog[S2][A1] + p.sec_exit) * z[Y6TS2];
	  d[Y6DS2] = p.sec_exit * z[Y6TS2] - (auptake[A1] + p.prog[S2][A1] + p.revert[NONE][S2]) * z[Y6DS2] + p.prog[S2][L3] * z[Y5DS2];
	  d[Y6Ds2] = p.revert[NONE][S2] * z[Y6DS2] - (auptake[A1] + p.prog[WS2][A1]) * z[Y6Ds2] + p.prog[WS2][L3] * z[Y5Ds2];

	  // update cumulative HIV-related deaths with those occuring with CD4<200
	  mort
	    = p.prog[WT][A1] * (z[Y6NWT] + z[Y6HWT] + z[Y6PWT] + z[Y6FWT] + z[Y6DWT] + z[Y6UWT])
	    + p.prog[R1][A1] * (z[Y6NR1] + z[Y6HR1] + z[Y6PR1] + z[Y6FR1] + z[Y6DR1] + z[Y6UR1])
	    + p.prog[C1][A1] * (z[Y6NC1] + z[Y6HC1] + z[Y6PC1] + z[Y6FC1] + z[Y6DC1] + z[Y6UC1])
	    + p.prog[S1][A1] * (z[Y6NS1] + z[Y6HS1] + z[Y6PS1] + z[Y6FS1] + z[Y6DS1] + z[Y6US1])
	    + p.prog[Q1][A1] * (z[Y6NQ1] + z[Y6HQ1] + z[Y6PQ1] + z[Y6FQ1] + z[Y6DQ1] + z[Y6UQ1])
	    + p.prog[Q2][A1] * (z[Y6NQ2] + z[Y6HQ2] + z[Y6PQ2] + z[Y6FQ2] + z[Y6DQ2] + z[Y6UQ2])
	    + p.prog[R2][A1] * (z[Y6NR2S] + z[Y6LR2S] + z[Y6DR2S] + z[Y6UR2S] + z[Y6NR2R] + z[Y6LR2R] + z[Y6DR2R] + z[Y6UR2R])
	    + p.prog[C2][A1] * (z[Y6NC2S] + z[Y6LC2S] + z[Y6DC2S] + z[Y6UC2S] + z[Y6NC2R] + z[Y6LC2R] + z[Y6DC2R] + z[Y6UC2R])
	    + p.prog[WR1][A1] * (z[Y6Nr1] + z[Y6Hr1] + z[Y6Pr1] + z[Y6Fr1] + z[Y6Dr1] + z[Y6Ur1])
	    + p.prog[WC1][A1] * (z[Y6Nc1] + z[Y6Hc1] + z[Y6Pc1] + z[Y6Fc1] + z[Y6Dc1] + z[Y6Uc1])
	    + p.prog[WS1][A1] * (z[Y6Ns1] + z[Y6Hs1] + z[Y6Ps1] + z[Y6Fs1] + z[Y6Ds1] + z[Y6Us1])
	    + p.prog[WQ1][A1] * (z[Y6Nq1] + z[Y6Hq1] + z[Y6Pq1] + z[Y6Fq1] + z[Y6Dq1] + z[Y6Uq1])
	    + p.prog[WQ2][A1] * (z[Y6Nq2] + z[Y6Fq2] + z[Y6Dq2] + z[Y6Uq2])
	    + p.prog[WR2][A1] * (z[Y6Nr2S] + z[Y6Dr2S] + z[Y6Ur2S] + z[Y6Nr2R] + z[Y6Dr2R] + z[Y6Ur2R])
	    + p.prog[WC2][A1] * (z[Y6Nc2S] + z[Y6Dc2S] + z[Y6Uc2S] + z[Y6Nc2R] + z[Y6Dc2R] + z[Y6Uc2R])
	    + p.prog[S2][A1] * (z[Y6TS2] + z[Y6DS2]) // 2014-09-24
	    + p.prog[WS2][A1] * z[Y6Ds2] // 2014-09-24
	    + p.art_mortE[A1] * (z[Y6EWT] + z[Y6ER1] + z[Y6EC1] + z[Y6ES1] + z[Y6EQ1] + z[Y6EQ2] + z[Y6Er1] + z[Y6Ec1] + z[Y6Es1] + z[Y6Eq1] + z[Y6Eq2])
	    + p.art_mortL[A1] * (z[Y6LWT] + z[Y6LR1] + z[Y6LC1] + z[Y6LS1] + z[Y6LQ1] + z[Y6LQ2] + z[Y6Lr1] + z[Y6Lc1] + z[Y6Ls1] + z[Y6Lq1] + z[Y6Lq2])
	    + p.sec_mort[A1] * (z[Y6TR2S] + z[Y6TR2R] + z[Y6TC2S] + z[Y6TC2R] + z[Y6Tr2S] + z[Y6Tr2R] + z[Y6Tc2S] + z[Y6Tc2R])
	    + p.sec_mort[A1] * (z[Y6TWT] + z[Y6TR1] + z[Y6TC1] + z[Y6TS1] + z[Y6TQ1] + z[Y6TQ2] + z[Y6Tr1] + z[Y6Tc1] + z[Y6Ts1] + z[Y6Tq1] + z[Y6Tq2]); // 2014-10-20 DONE

	  if (b < BANDS - 1) {
	    dy[ACTIVE_HIV_MORT] += mort;
	    dy[ACTIVE_W_HIV_MORT] += mort * (g == F);
	    dy[ACTIVE_M_HIV_MORT] += mort * (g != F);

	    dy[ACTDSC_HIV_MORT] += discount * mort;
	    dy[ACTDSC_W_HIV_MORT] += discount * mort * (g == F);
	    dy[ACTDSC_M_HIV_MORT] += discount * mort * (g != F);
	  }

	  if (u == ANTE) {
	    dy[COHORT_HIV_MORT] += mort;
	    dy[COHORT_W_HIV_MORT] += mort * (g == F);
	    dy[COHORT_M_HIV_MORT] += mort * (g != F);

	    dy[COHDSC_HIV_MORT] += discount * mort;
	    dy[COHDSC_W_HIV_MORT] += discount * mort * (g == F);
	    dy[COHDSC_M_HIV_MORT] += discount * mort * (g != F);
	  }
	}

	// +=+ susceptible and window-period infection +=========================+
	for (int g(0); g < SEXCIRC; ++g) {
	  z = offset(g, k, a, u) + y;
	  d = offset(g, k, a, u) + dy;

	  d[XN] = prep_exit[sex[g]][b][k] * (z[XH] + z[XP]) - (lambda[g][b][k][NONE][WT] + lambda[g][b][k][NONE][Q1] + lambda[g][b][k][NONE][R1] + lambda[g][b][k][NONE][C1] + lambda[g][b][k][NONE][S1]) * z[XN];
	  d[XH] = prep_propH * (prep_inject[sex[g]][b][k] * z[XP]) - (lambda[g][b][k][PREP_HIGH][WT] + lambda[g][b][k][PREP_HIGH][Q1] + lambda[g][b][k][PREP_HIGH][R1] + lambda[g][b][k][PREP_HIGH][C1] + lambda[g][b][k][PREP_HIGH][S1] + prep_inject[sex[g]][b][k] * prep_propP + prep_exit[sex[g]][b][k]) * z[XH];
	  d[XP] = prep_propP * (prep_inject[sex[g]][b][k] * z[XH]) - (lambda[g][b][k][PREP_POOR][WT] + lambda[g][b][k][PREP_POOR][Q1] + lambda[g][b][k][PREP_POOR][R1] + lambda[g][b][k][PREP_POOR][C1] + lambda[g][b][k][PREP_POOR][S1] + prep_inject[sex[g]][b][k] * prep_propH + prep_exit[sex[g]][b][k]) * z[XP];

	  d[Y1NWT] = lambda[g][b][k][NONE][WT] * z[XN] + prep_exit[sex[g]][b][k] * (z[Y1HWT] + z[Y1PWT]) - (p.prog[WT][P1]) * z[Y1NWT];
	  d[Y1HWT] = lambda[g][b][k][PREP_HIGH][WT] * z[XH] + prep_propH * (prep_inject[sex[g]][b][k] * z[Y1PWT]) - (p.prog[WT][P1] + p.emerge[PREP_HIGH][WT] + prep_inject[sex[g]][b][k] * prep_propP + prep_exit[sex[g]][b][k]) * z[Y1HWT];
	  d[Y1PWT] = lambda[g][b][k][PREP_POOR][WT] * z[XP] + prep_propP * (prep_inject[sex[g]][b][k] * z[Y1HWT]) - (p.prog[WT][P1] + p.emerge[PREP_POOR][WT] + prep_inject[sex[g]][b][k] * prep_propH + prep_exit[sex[g]][b][k]) * z[Y1PWT];

	  d[Y1NR1] = lambda[g][b][k][NONE][R1] * z[XN] + prep_exit[sex[g]][b][k] * (z[Y1HR1] + z[Y1PR1]) - (p.prog[R1][P1] + p.revert[NONE][R1]) * z[Y1NR1];
	  d[Y1HR1] = lambda[g][b][k][PREP_HIGH][R1] * z[XH] + prep_propH * prep_inject[sex[g]][b][k] * z[Y1PR1] - (p.prog[R1][P1] + p.revert[PREP_HIGH][R1] + p.emerge[PREP_HIGH][R1] + prep_inject[sex[g]][b][k] * prep_propP + prep_exit[sex[g]][b][k]) * z[Y1HR1];
	  d[Y1PR1] = lambda[g][b][k][PREP_POOR][R1] * z[XP] + prep_propP * prep_inject[sex[g]][b][k] * z[Y1HR1] - (p.prog[R1][P1] + p.revert[PREP_POOR][R1] + p.emerge[PREP_POOR][R1] + prep_inject[sex[g]][b][k] * prep_propH + prep_exit[sex[g]][b][k]) * z[Y1PR1];      

	  d[Y1NC1] = lambda[g][b][k][NONE][C1] * z[XN] + prep_exit[sex[g]][b][k] * (z[Y1HC1] + z[Y1PC1]) - (p.prog[C1][P1] + p.revert[NONE][C1]) * z[Y1NC1];
	  d[Y1HC1] = lambda[g][b][k][PREP_HIGH][C1] * z[XH] + prep_propH * prep_inject[sex[g]][b][k] * z[Y1PC1] - (p.prog[C1][P1] + prep_inject[sex[g]][b][k] * prep_propP + prep_exit[sex[g]][b][k]) * z[Y1HC1];
	  d[Y1PC1] = lambda[g][b][k][PREP_POOR][C1] * z[XP] + prep_propP * prep_inject[sex[g]][b][k] * z[Y1HC1] - (p.prog[C1][P1] + prep_inject[sex[g]][b][k] * prep_propH + prep_exit[sex[g]][b][k]) * z[Y1PC1];      

	  d[Y1NQ1] = lambda[g][b][k][NONE][Q1] * z[XN] + prep_exit[sex[g]][b][k] * (z[Y1HQ1] + z[Y1PQ1]) - (p.prog[Q1][P1] + p.revert[NONE][Q1]) * z[Y1NQ1];
	  d[Y1HQ1] = lambda[g][b][k][PREP_HIGH][Q1] * z[XH] + prep_propH * prep_inject[sex[g]][b][k] * z[Y1PQ1] - (p.prog[Q1][P1] + prep_inject[sex[g]][b][k] * prep_propP + prep_exit[sex[g]][b][k]) * z[Y1HQ1];
	  d[Y1PQ1] = lambda[g][b][k][PREP_POOR][Q1] * z[XP] + prep_propP * prep_inject[sex[g]][b][k] * z[Y1HQ1] - (p.prog[Q1][P1] + prep_inject[sex[g]][b][k] * prep_propH + prep_exit[sex[g]][b][k]) * z[Y1PQ1];

	  d[Y1NS1] = lambda[g][b][k][NONE][S1] * z[XN] + prep_exit[sex[g]][b][k] * (z[Y1HS1] + z[Y1PS1]) - (p.prog[S1][P1] + p.revert[NONE][S1]) * z[Y1NS1];
	  d[Y1HS1] = lambda[g][b][k][PREP_HIGH][S1] * z[XH] + prep_propH * prep_inject[sex[g]][b][k] * z[Y1PS1] - (p.prog[S1][P1] + p.revert[PREP_HIGH][S1] + p.emerge[PREP_HIGH][S1] + prep_inject[sex[g]][b][k] * prep_propP + prep_exit[sex[g]][b][k]) * z[Y1HS1];
	  d[Y1PS1] = lambda[g][b][k][PREP_POOR][S1] * z[XP] + prep_propP * prep_inject[sex[g]][b][k] * z[Y1HS1] - (p.prog[S1][P1] + p.revert[PREP_POOR][S1] + p.emerge[PREP_POOR][S1] + prep_inject[sex[g]][b][k] * prep_propH + prep_exit[sex[g]][b][k]) * z[Y1PS1];

	  const double inci[] = {
	    lambda[g][b][k][NONE][WT] * z[XN] + lambda[g][b][k][PREP_HIGH][WT] * z[XH] + lambda[g][b][k][PREP_POOR][WT] * z[XP],
	    lambda[g][b][k][NONE][R1] * z[XN] + lambda[g][b][k][PREP_HIGH][R1] * z[XH] + lambda[g][b][k][PREP_POOR][R1] * z[XP],
	    lambda[g][b][k][NONE][C1] * z[XN] + lambda[g][b][k][PREP_HIGH][C1] * z[XH] + lambda[g][b][k][PREP_POOR][C1] * z[XP],
	    lambda[g][b][k][NONE][Q1] * z[XN] + lambda[g][b][k][PREP_HIGH][Q1] * z[XH] + lambda[g][b][k][PREP_POOR][Q1] * z[XP],
	    lambda[g][b][k][NONE][S1] * z[XN] + lambda[g][b][k][PREP_HIGH][S1] * z[XH] + lambda[g][b][k][PREP_POOR][S1] * z[XP]};

	  // update new infection accumulators
	  if (b < BANDS - 1) {
	    dy[ACTIVE_WT] += inci[WT];
	    dy[ACTIVE_R1] += inci[R1];
	    dy[ACTIVE_C1] += inci[C1];
	    dy[ACTIVE_Q1] += inci[Q1];
	    dy[ACTIVE_S1] += inci[S1];

	    dy[ACTIVE_W_WT] += inci[WT] * (g == F);
	    dy[ACTIVE_W_R1] += inci[R1] * (g == F);
	    dy[ACTIVE_W_C1] += inci[C1] * (g == F);
	    dy[ACTIVE_W_Q1] += inci[Q1] * (g == F);
	    dy[ACTIVE_W_S1] += inci[S1] * (g == F);

	    dy[ACTIVE_M_WT] += inci[WT] * (g != F);
	    dy[ACTIVE_M_R1] += inci[R1] * (g != F);
	    dy[ACTIVE_M_C1] += inci[C1] * (g != F);
	    dy[ACTIVE_M_Q1] += inci[Q1] * (g != F);
	    dy[ACTIVE_M_S1] += inci[S1] * (g != F);

	    dy[ACTDSC_WT] += discount * inci[WT];
	    dy[ACTDSC_R1] += discount * inci[R1];
	    dy[ACTDSC_C1] += discount * inci[C1];
	    dy[ACTDSC_Q1] += discount * inci[Q1];
	    dy[ACTDSC_S1] += discount * inci[S1];

	    dy[ACTDSC_W_WT] += discount * inci[WT] * (g == F);
	    dy[ACTDSC_W_R1] += discount * inci[R1] * (g == F);
	    dy[ACTDSC_W_C1] += discount * inci[C1] * (g == F);
	    dy[ACTDSC_W_Q1] += discount * inci[Q1] * (g == F);
	    dy[ACTDSC_W_S1] += discount * inci[S1] * (g == F);

	    dy[ACTDSC_M_WT] += discount * inci[WT] * (g != F);
	    dy[ACTDSC_M_R1] += discount * inci[R1] * (g != F);
	    dy[ACTDSC_M_C1] += discount * inci[C1] * (g != F);
	    dy[ACTDSC_M_Q1] += discount * inci[Q1] * (g != F);
	    dy[ACTDSC_M_S1] += discount * inci[S1] * (g != F);
	  }

	  if (u == ANTE) {
	    dy[COHORT_WT] += inci[WT];
	    dy[COHORT_R1] += inci[R1];
	    dy[COHORT_C1] += inci[C1];
	    dy[COHORT_Q1] += inci[Q1];
	    dy[COHORT_S1] += inci[S1];

	    dy[COHORT_W_WT] += inci[WT] * (g == F);
	    dy[COHORT_W_R1] += inci[R1] * (g == F);
	    dy[COHORT_W_C1] += inci[C1] * (g == F);
	    dy[COHORT_W_Q1] += inci[Q1] * (g == F);
	    dy[COHORT_W_S1] += inci[S1] * (g == F);

	    dy[COHORT_M_WT] += inci[WT] * (g != F);
	    dy[COHORT_M_R1] += inci[R1] * (g != F);
	    dy[COHORT_M_C1] += inci[C1] * (g != F);
	    dy[COHORT_M_Q1] += inci[Q1] * (g != F);
	    dy[COHORT_M_S1] += inci[S1] * (g != F);

	    dy[COHDSC_WT] += discount * inci[WT];
	    dy[COHDSC_R1] += discount * inci[R1];
	    dy[COHDSC_C1] += discount * inci[C1];
	    dy[COHDSC_Q1] += discount * inci[Q1];
	    dy[COHDSC_S1] += discount * inci[S1];

	    dy[COHDSC_W_WT] += discount * inci[WT] * (g == F);
	    dy[COHDSC_W_R1] += discount * inci[R1] * (g == F);
	    dy[COHDSC_W_C1] += discount * inci[C1] * (g == F);
	    dy[COHDSC_W_Q1] += discount * inci[Q1] * (g == F);
	    dy[COHDSC_W_S1] += discount * inci[S1] * (g == F);

	    dy[COHDSC_M_WT] += discount * inci[WT] * (g != F);
	    dy[COHDSC_M_R1] += discount * inci[R1] * (g != F);
	    dy[COHDSC_M_C1] += discount * inci[C1] * (g != F);
	    dy[COHDSC_M_Q1] += discount * inci[Q1] * (g != F);
	    dy[COHDSC_M_S1] += discount * inci[S1] * (g != F);
	  }
	}
      }
    }
  }

  // 2. handle demographic processes
  for (int u(0); u < TRACKS; ++u) {
    for (int k(0); k < LEVELS; ++k) {
      for (int g(0); g < SEXCIRC; ++g) {
	z = offset(g, k, 0, u) + y;
	d = offset(g, k, 0, u) + dy;
	for (int s(0); s < STATES; ++s) {
	  d[s] -= (p.rate_mort[sex[g]][0] + rate_age[0]) * z[s];
	}
	for (int a(1); a < AGES; ++a) {
	  b = band[a];
	  w = z; // previous age compartment magnitude
	  z = offset(g, k, a, u) + y; // current age compartment magnitude
	  d = offset(g, k, a, u) + dy; // current age compartment derivative
	  for (int s(0); s < STATES; ++s) {
	    d[s] += rate_age[a-1] * w[s] - (p.rate_mort[sex[g]][b] + rate_age[a]) * z[s];
	  }
	}
      }

      // People are tracked seperately if they debut before versus
      // after PrEP enrollment closes to allow tracking of lifetime
      // outcomes among people alive during PrEP rollout
      if ((u == 0 && t < p.prep_time_closure) || (u == 1 && t >= p.prep_time_closure)) {
	// sexual debut: susceptible women aged 15
	z = offset(F, k, 0, u) + y;
	d = offset(F, k, 0, u) + dy;
	d[XN] += debut[F][k];

	// sexual debut: uncircumcised, susceptible men aged 15
	z = offset(MU, k, 0, u) + y;
	d = offset(MU, k, 0, u) + dy;
	d[XN] += debut[M][k] * (1.0 - mmc_prop);

	// sexual debut: circumcised, susceptible men aged 15
	z = offset(MC, k, 0, u) + y;
	d = offset(MC, k, 0, u) + dy;
	d[XN] += debut[M][k] * mmc_prop;
      }
    }
  }

  // Male medical circumcision
  rates_mmc(t, y, dy, p, cuptake);
  for (int u(0); u < TRACKS; ++u) {
    for (int a(0); a < AGES-1; ++a) { // No MMC at ages >= 55
      b = band[a];
      for (int k(0); k < LEVELS; ++k) {
	z = offset(MU, k, a, u) + y;

	d = offset(MU, k, a, u) + dy;
	for (int s(0); s < STATES; ++s) d[s] -= cuptake[b] * z[s];

	d = offset(MC, k, a, u) + dy;
	for (int s(0); s < STATES; ++s) d[s] += cuptake[b] * z[s];
      }
    }
  }

  if (p.replace == REPLACE_NONE) {
    // no action required
  } else if (p.replace == REPLACE_ALL) {
    // Men and women in lower sexual activity levels move into higher
    // activity levels at rates sufficient to maintain a constant
    // activity level distribution throughout simulation

    // Simulate exit from sex work
    for (int u(0); u < TRACKS; ++u) {
      for (int a(0); a < AGES-1; ++a) {
	z = offset(F, HIGH, a, u) + y;
	d = offset(F, HIGH, a, u) + dy;
	for (int s(0); s < STATES; ++s) d[s] -= p.prop_risk_exit[s] * p.rate_csw_exit * z[s];
	for (int k(0); k < HIGH; ++k) {
	  d = offset(F, k, a, u) + dy;
	  const double prop(p.prop_risk[F][k] / (p.prop_sex[F] - p.prop_risk[F][HIGH]));
	  for (int s(0); s < STATES; ++s) {
	    d[s] += p.prop_risk_exit[s] * p.rate_csw_exit * prop * z[s];
	  }
	}
      }
    }

    double rate_replace[SEXES][LEVELS];
    rates_replace(t, y, dy, p, rate_replace);

    // update derivatives using replenishment rates
    for (int u(0); u < TRACKS; ++u) {
      for (int g(0); g < SEXCIRC; ++g) {
	for (int a(0); a < AGES-1; ++a) {
	  for (int k(0); k < HIGH; ++k) {
	    z = offset(g, k, a, u) + y;
	    d = offset(g, k, a, u) + dy;
	    for (int s(0); s < STATES; ++s) d[s] -= p.prop_risk_init[s] * rate_replace[sex[g]][k] * z[s];

	    d = offset(g, k+1, a, u) + dy;
	    for (int s(0); s < STATES; ++s) d[s] += p.prop_risk_init[s] * rate_replace[sex[g]][k] * z[s];
	  }
	}
      }
    }

  } else {
    double csw_init[BANDS];

    // Simulate exit from sex work
    for (int u(0); u < TRACKS; ++u) {
      for (int a(0); a < AGES-1; ++a) {
	z = offset(F, HIGH, a, u) + y;
	d = offset(F, HIGH, a, u) + dy;
	for (int s(0); s < STATES; ++s) d[s] -= p.prop_risk_exit[s] * p.rate_csw_exit * z[s];
	for (int k(0); k < HIGH; ++k) {
	  d = offset(F, k, a, u) + dy;
	  const double prop(p.prop_risk[F][k] / (p.prop_sex[F] - p.prop_risk[F][HIGH]));
	  for (int s(0); s < STATES; ++s) {
	    d[s] += p.prop_risk_exit[s] * p.rate_csw_exit * prop * z[s];
	  }
	}
      }
    }

    if (p.replace == REPLACE_STATIC) {
      // Relative rates of sex work initiation are based on the age
      // distribution reported in rees2000sajs
      const double prop_csw_init[] = {1.0, 1.5, 0.9, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0};
      const double rate_csw_base(-log(1.0 - p.prop_risk[F][HIGH] / p.prop_sex[F]));
      for (int bw(0); bw < BANDS; ++bw) csw_init[bw] = rate_csw_base * prop_csw_init[bw];
    } else if (p.replace == REPLACE_CSW) {
      rates_csw_prop(t, y, dy, p, csw_init);
    } else if (p.replace == REPLACE_DEMAND) {
      rates_csw_need(t, y, dy, p, csw_init);
    } else {
      fprintf(stderr, "Warning: replacement scheme %d not implemented\n", p.replace);
      std::fill(csw_init, csw_init + BANDS, 0.0);
    }

    // negative sex work initiation rates are disallowed, since
    // negative rates encode additional exit from sex work. The code
    // is not designed to allow this; exit rates must be scaled by
    // p.prop_risk_exit instead of p.prop_risk_init
    for (b = 0; b < BANDS; ++b) assert(csw_init[b] >= 0.0);

    // update derivatives using replenishment rates
    for (int u(0); u < TRACKS; ++u) {
      for (int k(0); k < HIGH; ++k) {
	for (int a(0); a < AGES-1; ++a) {
	  b = band[a];
	  z = offset(F, k, a, u) + y;
	  d = offset(F, k, a, u) + dy;
	  for (int s(0); s < STATES; ++s) d[s] -= p.prop_risk_init[s] * csw_init[b] * z[s];

	  d = offset(F, HIGH, a, u) + dy;
	  for (int s(0); s < STATES; ++s) d[s] += p.prop_risk_init[s] * csw_init[b] * z[s];
	}
      }
    }

  }

  // BEGIN PREP UPTAKE
  rates_prep(t, y, dy, p, puptake);
  for (int u(0); u < TRACKS; ++u) {
    for (int g(0); g < SEXCIRC; ++g) {
      for (int a(0); a < AGES; ++a) {
	b = band[a];
	for (int k(0); k < LEVELS; ++k) {
	  z = offset(g, k, a, u) + y;
	  d = offset(g, k, a, u) + dy;

	  assert(t < p.prep_time_closure || puptake[sex[g]][b][k] == 0.0);

	  d[XN] -= puptake[sex[g]][b][k] * z[XN];
	  d[XH] += prep_propH * puptake[sex[g]][b][k] * z[XN];
	  d[XP] += prep_propP * puptake[sex[g]][b][k] * z[XN];

	  d[Y1NWT] -= puptake[sex[g]][b][k] * z[Y1NWT];
	  d[Y1HWT] += prep_propH * puptake[sex[g]][b][k] * z[Y1NWT];
	  d[Y1PWT] += prep_propP * puptake[sex[g]][b][k] * z[Y1NWT];

	  d[Y1NR1] -= puptake[sex[g]][b][k] * z[Y1NR1];
	  d[Y1HR1] += prep_propH * puptake[sex[g]][b][k] * z[Y1NR1];
	  d[Y1PR1] += prep_propP * puptake[sex[g]][b][k] * z[Y1NR1];

	  d[Y1NC1] -= puptake[sex[g]][b][k] * z[Y1NC1];
	  d[Y1HC1] += prep_propH * puptake[sex[g]][b][k] * z[Y1NC1];
	  d[Y1PC1] += prep_propP * puptake[sex[g]][b][k] * z[Y1NC1];

	  d[Y1NQ1] -= puptake[sex[g]][b][k] * z[Y1NQ1];
	  d[Y1HQ1] += prep_propH * puptake[sex[g]][b][k] * z[Y1NQ1];
	  d[Y1PQ1] += prep_propP * puptake[sex[g]][b][k] * z[Y1NQ1];

	  d[Y1NS1] -= puptake[sex[g]][b][k] * z[Y1NS1];
	  d[Y1HS1] += prep_propH * puptake[sex[g]][b][k] * z[Y1NS1];
	  d[Y1PS1] += prep_propP * puptake[sex[g]][b][k] * z[Y1NS1];

	  d[Y1Nr1] -= puptake[sex[g]][b][k] * z[Y1Nr1];
	  d[Y1Hr1] += prep_propH * puptake[sex[g]][b][k] * z[Y1Nr1];
	  d[Y1Pr1] += prep_propP * puptake[sex[g]][b][k] * z[Y1Nr1];

	  d[Y1Nc1] -= puptake[sex[g]][b][k] * z[Y1Nc1];
	  d[Y1Hc1] += prep_propH * puptake[sex[g]][b][k] * z[Y1Nc1];
	  d[Y1Pc1] += prep_propP * puptake[sex[g]][b][k] * z[Y1Nc1];

	  d[Y1Nq1] -= puptake[sex[g]][b][k] * z[Y1Nq1];
	  d[Y1Hq1] += prep_propH * puptake[sex[g]][b][k] * z[Y1Nq1];
	  d[Y1Pq1] += prep_propP * puptake[sex[g]][b][k] * z[Y1Nq1];

	  d[Y1Ns1] -= puptake[sex[g]][b][k] * z[Y1Ns1];
	  d[Y1Hs1] += prep_propH * puptake[sex[g]][b][k] * z[Y1Ns1];
	  d[Y1Ps1] += prep_propP * puptake[sex[g]][b][k] * z[Y1Ns1];

	  d[Y1Nq2] -= puptake[sex[g]][b][k] * z[Y1Nq2];
	  d[Y1NQ2] -= puptake[sex[g]][b][k] * z[Y1NQ2];
	  d[Y1HQ2] += prep_propH * puptake[sex[g]][b][k] * (z[Y1NQ2] + z[Y1Nq2]);
	  d[Y1PQ2] += prep_propP * puptake[sex[g]][b][k] * (z[Y1NQ2] + z[Y1Nq2]);
	}
      }
    }
  }
  // END PREP UPTAKE

  // update life-years and person-years of ART and PrEP
  double bg_mort, py_prp, a_init;
  double py_art[STAGES], py_2nd[STAGES], lyears[STAGES];
  double py_arv[STAGES][ARV_STATES];

  for (int g(0); g < SEXCIRC; ++g) {
    for (int a(0); a < AGES; ++a) {
      b = band[a];
      for (int k(0); k < LEVELS; ++k) {
	z = offset(g, k, a, ANTE) + y;

	bg_mort = p.rate_mort[sex[g]][b] * std::accumulate(z, z + STATES, 0.0);
	py_prp = 0.0;
	a_init = 0.0;

	for (int h(0); h < STAGES; ++h) std::fill(py_arv[h], py_arv[h] + ARV_STATES, 0.0);
	for (int s(Y1NWT); s < STATES; ++s) py_arv[attr_stage[s]][attr_drug[s]] += z[s];
	for (int h(0); h < STAGES; ++h) {
	  lyears[h] = std::accumulate(py_arv[h], py_arv[h] + ARV_STATES, 0.0);
	  py_art[h] = py_arv[h][ART1_NEW] + py_arv[h][ART1_OLD] + py_arv[h][ART1_NA];
	  py_2nd[h] = py_arv[h][ART2] + py_arv[h][ART2_NA];
	  py_prp += py_arv[h][PREP_HIGH] + py_arv[h][PREP_POOR];
	  a_init += py_arv[h][NONE] * auptake[h];
	}

	if (b < BANDS - 1) {
	  dy[ACTIVE_LYL_X] += z[XN] + z[XH] + z[XP];
	  dy[ACTIVE_PRP] += z[XH] + z[XP] + py_prp;
	  dy[ACTIVE_BAK_MORT] += bg_mort;
	  dy[ACTIVE_AINIT] += a_init;

	  dy[ACTDSC_LYL_X] += discount * (z[XN] + z[XH] + z[XP]);
	  dy[ACTDSC_PRP] += discount * (z[XH] + z[XP] + py_prp);
	  dy[ACTDSC_BAK_MORT] += discount * bg_mort;
	  dy[ACTDSC_AINIT] += discount * a_init;
	}

	dy[COHORT_LYL_X] += z[XN] + z[XH] + z[XP];
	dy[COHORT_PRP] += z[XH] + z[XP] + py_prp;
	dy[COHORT_BAK_MORT] += bg_mort;
	dy[COHORT_AINIT] += a_init;

	dy[COHDSC_LYL_X] += discount * (z[XN] + z[XH] + z[XP]);
	dy[COHDSC_PRP] += discount * (z[XH] + z[XP] + py_prp);
	dy[COHDSC_BAK_MORT] += discount * bg_mort;
	dy[COHDSC_AINIT] += discount * a_init;

 	for (int h(0); h < STAGES; ++h) {
 	  dy[active_lyl[h]] += (b < BANDS - 1) * lyears[h];
 	  dy[active_art[h]] += (b < BANDS - 1) * py_art[h];
 	  dy[active_2nd[h]] += (b < BANDS - 1) * py_2nd[h];
 	  dy[cohort_lyl[h]] += lyears[h];
 	  dy[cohort_art[h]] += py_art[h];
 	  dy[cohort_2nd[h]] += py_2nd[h];

 	  dy[actdsc_lyl[h]] += discount * (b < BANDS - 1) * lyears[h];
 	  dy[actdsc_art[h]] += discount * (b < BANDS - 1) * py_art[h];
 	  dy[actdsc_2nd[h]] += discount * (b < BANDS - 1) * py_2nd[h];
 	  dy[cohdsc_lyl[h]] += discount * lyears[h];
 	  dy[cohdsc_art[h]] += discount * py_art[h];
 	  dy[cohdsc_2nd[h]] += discount * py_2nd[h];
 	}

	if (g == F) {
	  if (b < BANDS - 1) {
	    dy[ACTIVE_W_LYL_X] += z[XN] + z[XH] + z[XP];
	    dy[ACTIVE_W_PRP] += z[XH] + z[XP] + py_prp;
	    dy[ACTIVE_W_BAK_MORT] += bg_mort;
	    dy[ACTIVE_W_AINIT] += a_init;

	    dy[ACTDSC_W_LYL_X] += discount * (z[XN] + z[XH] + z[XP]);
	    dy[ACTDSC_W_PRP] += discount * (z[XH] + z[XP] + py_prp);
	    dy[ACTDSC_W_BAK_MORT] += discount * bg_mort;
	    dy[ACTDSC_W_AINIT] += discount * a_init;
	  }
	  dy[COHORT_W_LYL_X] += z[XN] + z[XH] + z[XP];
	  dy[COHORT_W_PRP] += z[XH] + z[XP] + py_prp;
	  dy[COHORT_W_BAK_MORT] += bg_mort;
	  dy[COHORT_W_AINIT] += a_init;

	  dy[COHDSC_W_LYL_X] += discount * (z[XN] + z[XH] + z[XP]);
	  dy[COHDSC_W_PRP] += discount * (z[XH] + z[XP] + py_prp);
	  dy[COHDSC_W_BAK_MORT] += discount * bg_mort;
	  dy[COHDSC_W_AINIT] += discount * a_init;
	  for (int h(0); h < STAGES; ++h) {
	    dy[active_w_lyl[h]] += (b < BANDS - 1) * lyears[h];
	    dy[active_w_art[h]] += (b < BANDS - 1) * py_art[h];
	    dy[active_w_2nd[h]] += (b < BANDS - 1) * py_2nd[h];
	    dy[cohort_w_lyl[h]] += lyears[h];
	    dy[cohort_w_art[h]] += py_art[h];
	    dy[cohort_w_2nd[h]] += py_2nd[h];

	    dy[actdsc_w_lyl[h]] += discount * (b < BANDS - 1) * lyears[h];
	    dy[actdsc_w_art[h]] += discount * (b < BANDS - 1) * py_art[h];
	    dy[actdsc_w_2nd[h]] += discount * (b < BANDS - 1) * py_2nd[h];
	    dy[cohdsc_w_lyl[h]] += discount * lyears[h];
	    dy[cohdsc_w_art[h]] += discount * py_art[h];
	    dy[cohdsc_w_2nd[h]] += discount * py_2nd[h];
	  }
	} else {
	  if (b < BANDS - 1) {
	    dy[ACTIVE_M_LYL_X] += z[XN] + z[XH] + z[XP];
	    dy[ACTIVE_M_PRP] += z[XH] + z[XP] + py_prp;
	    dy[ACTIVE_M_BAK_MORT] += bg_mort;
	    dy[ACTIVE_M_AINIT] += a_init;

	    dy[ACTDSC_M_LYL_X] += discount * (z[XN] + z[XH] + z[XP]);
	    dy[ACTDSC_M_PRP] += discount * (z[XH] + z[XP] + py_prp);
	    dy[ACTDSC_M_BAK_MORT] += discount * bg_mort;
	    dy[ACTDSC_M_AINIT] += discount * a_init;
	  }
	  dy[COHORT_M_LYL_X] += z[XN] + z[XH] + z[XP];
	  dy[COHORT_M_PRP] += z[XH] + z[XP] + py_prp;
	  dy[COHORT_M_BAK_MORT] += bg_mort;
	  dy[COHORT_M_AINIT] += a_init;

	  dy[COHDSC_M_LYL_X] += discount * (z[XN] + z[XH] + z[XP]);
	  dy[COHDSC_M_PRP] += discount * (z[XH] + z[XP] + py_prp);
	  dy[COHDSC_M_BAK_MORT] += discount * bg_mort;
	  dy[COHDSC_M_AINIT] += discount * a_init;
	  for (int h(0); h < STAGES; ++h) {
	    dy[active_m_lyl[h]] += (b < BANDS - 1) * lyears[h];
	    dy[active_m_art[h]] += (b < BANDS - 1) * py_art[h];
	    dy[active_m_2nd[h]] += (b < BANDS - 1) * py_2nd[h];
	    dy[cohort_m_lyl[h]] += lyears[h];
	    dy[cohort_m_art[h]] += py_art[h];
	    dy[cohort_m_2nd[h]] += py_2nd[h];

	    dy[actdsc_m_lyl[h]] += discount * (b < BANDS - 1) * lyears[h];
	    dy[actdsc_m_art[h]] += discount * (b < BANDS - 1) * py_art[h];
	    dy[actdsc_m_2nd[h]] += discount * (b < BANDS - 1) * py_2nd[h];
	    dy[cohdsc_m_lyl[h]] += discount * lyears[h];
	    dy[cohdsc_m_art[h]] += discount * py_art[h];
	    dy[cohdsc_m_2nd[h]] += discount * py_2nd[h];
	  }
	}

	if (b < BANDS - 1) {
	  z = offset(g, k, a, POST) + y;

	  bg_mort = p.rate_mort[sex[g]][b] * std::accumulate(z, z + STATES, 0.0);
	  py_prp = 0.0;
	  a_init = 0.0;

	  for (int h(0); h < STAGES; ++h) std::fill(py_arv[h], py_arv[h] + ARV_STATES, 0.0);
	  for (int s(Y1NWT); s < STATES; ++s) py_arv[attr_stage[s]][attr_drug[s]] += z[s];
	  for (int h(0); h < STAGES; ++h) {
	    lyears[h] = std::accumulate(py_arv[h], py_arv[h] + ARV_STATES, 0.0);
	    py_art[h] = py_arv[h][ART1_NEW] + py_arv[h][ART1_OLD] + py_arv[h][ART1_NA];
	    py_2nd[h] = py_arv[h][ART2] + py_arv[h][ART2_NA];
	    py_prp += py_arv[h][PREP_HIGH] + py_arv[h][PREP_POOR];
	    a_init += py_arv[h][NONE] * auptake[h];
	  }

	  dy[ACTIVE_LYL_X] += z[XN] + z[XH] + z[XP];
	  dy[ACTIVE_PRP] += z[XH] + z[XP] + py_prp;
	  dy[ACTIVE_BAK_MORT] += bg_mort;
	  dy[ACTIVE_AINIT] += a_init;

	  dy[ACTDSC_LYL_X] += discount * (z[XN] + z[XH] + z[XP]);
	  dy[ACTDSC_PRP] += discount * (z[XH] + z[XP] + py_prp);
	  dy[ACTDSC_BAK_MORT] += discount * bg_mort;
	  dy[ACTDSC_AINIT] += discount * a_init;
	  for (int h(0); h < STAGES; ++h) {
	    dy[active_lyl[h]] += lyears[h];
	    dy[active_art[h]] += py_art[h];
	    dy[active_2nd[h]] += py_2nd[h];
	    dy[actdsc_lyl[h]] += discount * lyears[h];
	    dy[actdsc_art[h]] += discount * py_art[h];
	    dy[actdsc_2nd[h]] += discount * py_2nd[h];
	  }

	  if (g == F) {
	    dy[ACTIVE_W_LYL_X] += z[XN] + z[XH] + z[XP];
	    dy[ACTIVE_W_PRP] += z[XH] + z[XP] + py_prp;
	    dy[ACTIVE_W_BAK_MORT] += bg_mort;
	    dy[ACTIVE_W_AINIT] += a_init;

	    dy[ACTDSC_W_LYL_X] += discount * (z[XN] + z[XH] + z[XP]);
	    dy[ACTDSC_W_PRP] += discount * (z[XH] + z[XP] + py_prp);
	    dy[ACTDSC_W_BAK_MORT] += discount * bg_mort;
	    dy[ACTDSC_W_AINIT] += discount * a_init;
	    for (int h(0); h < STAGES; ++h) {
	      dy[active_w_lyl[h]] += lyears[h];
	      dy[active_w_art[h]] += py_art[h];
	      dy[active_w_2nd[h]] += py_2nd[h];
	      dy[actdsc_w_lyl[h]] += discount * lyears[h];
	      dy[actdsc_w_art[h]] += discount * py_art[h];
	      dy[actdsc_w_2nd[h]] += discount * py_2nd[h];
	    }
	  } else {
	    dy[ACTIVE_M_LYL_X] += z[XN] + z[XH] + z[XP];
	    dy[ACTIVE_M_PRP] += z[XH] + z[XP] + py_prp;
	    dy[ACTIVE_M_BAK_MORT] += bg_mort;
	    dy[ACTIVE_M_AINIT] += a_init;

	    dy[ACTDSC_M_LYL_X] += discount * (z[XN] + z[XH] + z[XP]);
	    dy[ACTDSC_M_PRP] += discount * (z[XH] + z[XP] + py_prp);
	    dy[ACTDSC_M_BAK_MORT] += discount * bg_mort;
	    dy[ACTDSC_M_AINIT] += discount * a_init;
	    for (int h(0); h < STAGES; ++h) {
	      dy[active_m_lyl[h]] += lyears[h];
	      dy[active_m_art[h]] += py_art[h];
	      dy[active_m_2nd[h]] += py_2nd[h];
	      dy[actdsc_m_lyl[h]] += discount * lyears[h];
	      dy[actdsc_m_art[h]] += discount * py_art[h];
	      dy[actdsc_m_2nd[h]] += discount * py_2nd[h];
	    }
	  }
	}
      }
    }
  }

  // update circumcision accounting
  for (int k(0); k < LEVELS; ++k) {
    dy[ACTIVE_MMC_DEBUT] += debut[M][k] * mmc_prop;
    dy[COHORT_MMC_DEBUT] += debut[M][k] * mmc_prop * (t < p.prep_time_closure);
    dy[ACTDSC_MMC_DEBUT] += discount * debut[M][k] * mmc_prop;
    dy[COHDSC_MMC_DEBUT] += discount * debut[M][k] * mmc_prop * (t < p.prep_time_closure);
    for (int a(0); a < AGES; ++a) {
      b = band[a];
      z = offset(MU, k, a, ANTE) + y;
      dy[ACTIVE_MMC_UPTAKE] += (b < BANDS - 1) * cuptake[b] * std::accumulate(z, z + STATES, 0.0);
      dy[COHORT_MMC_UPTAKE] += cuptake[b] * std::accumulate(z, z + STATES, 0.0);

      dy[ACTDSC_MMC_UPTAKE] += discount * (b < BANDS - 1) * cuptake[b] * std::accumulate(z, z + STATES, 0.0);
      dy[COHDSC_MMC_UPTAKE] += discount * cuptake[b] * std::accumulate(z, z + STATES, 0.0);

      z = offset(MU, k, a, POST) + y;
      dy[ACTIVE_MMC_UPTAKE] += (b < BANDS - 1) * cuptake[b] * std::accumulate(z, z + STATES, 0.0);
      dy[ACTDSC_MMC_UPTAKE] += discount * (b < BANDS - 1) * cuptake[b] * std::accumulate(z, z + STATES, 0.0);
    }
  }

  return GSL_SUCCESS;
}

void force(double const t, double const y[],
	   Parameters const& p, double lambda[][BANDS][LEVELS][PREP_STATES][TRANSMIT]) {

  int g, a, b, k, s;
  double *w, *wm, *wf;
  double const* z;
  double mm, mf, contact;
  double numer, denom;

  // maps sex and circumcision status to sex (F->F, MU->M, MC->M)
  const Sex sex[] = {F, M, M};

  // this is used to cache the infectivity from a partner in a given
  // state that applies to either women and uncircumcised men
  // (beta[0]) or circumcised men (beta[1])
  //  double beta[2][STATES];

  // beta[u][s] cache per-partnership transmission probabilities given
  // the partnership's protection status u and the infected partner's
  // disease state s.
  // u=0 no MMC, no PrEP
  // u=1 no MMC, high PrEP adherence 
  // u=2 no MMC, poor PrEP adherence
  // u=3 MMC, no PrEP
  // u=4 MMC, high PrEP adherence 
  // u=5 MMC, poor PrEP adherence
  double beta[6][STATES];

  // eps is used to avoid divide-by-zero problems when evaluating
  // mixing terms. This is unlikely to come up except during
  // debugging, since mixing term denominators are zero only if a risk
  // group has contact rate 0 or is empty. When either of these are
  // false, eps should be rounded off during floating-point addition.
  const double eps(std::numeric_limits<double>::denorm_min());

  // bands stores the ODE state summed over ages within each age band
  const int nb(SEXES * BANDS * LEVELS * STATES);
  double bands[nb];

  double active[SEXES][BANDS][LEVELS];
  double supply[SEXES][BANDS][LEVELS];
  double supply_act[SEXES][LEVELS];
  double supply_age[SEXES][BANDS];
  double supply_sex[SEXES];

  // suminf[gi][ki][qi][kj][bj][vj] caches the cumulative
  // infectiousness presented to sex gi, level ki persons with PrEP
  // status qi from level kj, age bj partners with virus vj. This is
  // the sum of compartment size times per-partnership transmission
  // probability. Since this latter does not depend on the susceptible
  // partners' age in our model, we do not need to specify that here.
  double suminf[SEXCIRC][LEVELS][PREP_STATES][LEVELS][BANDS][TRANSMIT];

  double condom[LEVELS][LEVELS];
  double condom_factor;


  if (t < p.condom_change_time) {
    condom_factor = 1.0;
  } else if (t <= p.condom_change_time + p.condom_change_span) {
    condom_factor = 1 - (1 - p.condom_change) * (t - p.condom_change_time) / p.condom_change_span;
  } else {
    condom_factor = p.condom_change;
  }

  for (int kf(0); kf < LEVELS; ++kf) {
    for (int km(0); km < LEVELS; ++km) {
      condom[kf][km] = 1 - condom_factor * (1 - p.condom_prop[kf][km]);
    }
  }

  // distribution of coital acts by protection state, stratified by
  // partnership type and the susceptible partner's PrEP adherence
  // level
  // 0 -PrEP, -condom
  // 1 -PrEP, +condom
  // 2 +PrEP, -condom
  // 3 +PrEP, +condom
  double protect[LEVELS][LEVELS][PREP_STATES][4];
  for (int kf(0); kf < LEVELS; ++kf) {
    for (int km(0); km < LEVELS; ++km) {
      protect[kf][km][NONE][0] = 1.0 - condom[kf][km];
      protect[kf][km][NONE][1] = condom[kf][km];
      protect[kf][km][NONE][2] = 0.0;
      protect[kf][km][NONE][3] = 0.0;

      protect[kf][km][PREP_HIGH][0] = (1.0 - p.prep_adh_high) * (1.0 - condom[kf][km]);
      protect[kf][km][PREP_HIGH][1] = (1.0 - p.prep_adh_high) * condom[kf][km];
      protect[kf][km][PREP_HIGH][2] = p.prep_adh_high * (1.0 - condom[kf][km]);
      protect[kf][km][PREP_HIGH][3] = p.prep_adh_high * condom[kf][km];

      protect[kf][km][PREP_POOR][0] = (1.0 - p.prep_adh_poor) * (1.0 - condom[kf][km]);
      protect[kf][km][PREP_POOR][1] = (1.0 - p.prep_adh_poor) * condom[kf][km];
      protect[kf][km][PREP_POOR][2] = p.prep_adh_poor * (1.0 - condom[kf][km]);
      protect[kf][km][PREP_POOR][3] = p.prep_adh_poor * condom[kf][km];
    }
  }

  std::fill(bands, bands + nb, 0.0);

  for (g = 0; g < SEXCIRC; ++g) {
    for (k = 0; k < LEVELS; ++k) {
      for (b = 0; b < BANDS; ++b) {
	for (int q(0); q < PREP_STATES; ++q) {
	  std::fill(lambda[g][b][k][q], lambda[g][b][k][q] + TRANSMIT, 0.0);
	}
      }
    }
  }

  // collapse age groups into age bands. Mixing is random by
  // circumcision status, so circumcised and uncircumcised men are
  // aggregated here
  for (int u(0); u < TRACKS; ++u) {
    for (g = 0; g < SEXCIRC; ++g) {
      for (k = 0; k < LEVELS; ++k) {
	for (a = 0; a < AGES; ++a) {
	  b = band[a];
	  z = offset(g, k, a, u) + y;
	  w = STATES * (BANDS * (sex[g] * LEVELS + k) + b) + bands;
	  for (s = 0; s < STATES; ++s) w[s] += z[s];
	}
      }
    }
  }

  for (g = 0; g < SEXES; ++g) {
    for (k = 0; k < LEVELS; ++k) {
      for (b = 0; b < BANDS; ++b) {
	w = STATES * (BANDS * (g * LEVELS + k) + b) + bands;
	active[g][b][k] = std::accumulate(w, w + STATES, 0.0);
	supply[g][b][k] = p.partner[g][b][k] * std::inner_product(w, w + STATES, p.prop_risk_decr, 0.0);
      }
    }
  }

  for (g = 0; g < SEXES; ++g) {
    supply_sex[g] = 0.0;
    for (b = 0; b < BANDS; ++b) {
      supply_age[g][b] = 0.0;
      for (k = 0; k < LEVELS; ++k) supply_age[g][b] += supply[g][b][k];
    }
    for (k = 0; k < LEVELS; ++k) {
      supply_act[g][k] = 0.0;
      for (b = 0; b < BANDS; ++b) supply_act[g][k] += supply[g][b][k];
      supply_sex[g] += supply_act[g][k];
    }
  }

  for (int kf(0); kf < LEVELS; ++kf) {
    for (int km(0); km < LEVELS; ++km) {
      // calculate per-partnership transmission
      // probabilities. Exploits the identity a = 1 - exp(ln(1 - a))
      // to speed up calculation
      for (s = Y1NWT; s < STATES; ++s) {
	beta[0][s] = 1.0 - exp(p.acts[kf][km] * (  protect[kf][km][NONE     ][0] * p.infect[0][0][0][s]
						 + protect[kf][km][NONE     ][1] * p.infect[1][0][0][s]));
	beta[1][s] = 1.0 - exp(p.acts[kf][km] * (  protect[kf][km][PREP_HIGH][0] * p.infect[0][0][0][s]
						 + protect[kf][km][PREP_HIGH][1] * p.infect[1][0][0][s]
                                                 + protect[kf][km][PREP_HIGH][2] * p.infect[0][0][1][s]
						 + protect[kf][km][PREP_HIGH][3] * p.infect[1][0][1][s]));
	beta[2][s] = 1.0 - exp(p.acts[kf][km] * (  protect[kf][km][PREP_POOR][0] * p.infect[0][0][0][s]
						 + protect[kf][km][PREP_POOR][1] * p.infect[1][0][0][s]
                                                 + protect[kf][km][PREP_POOR][2] * p.infect[0][0][1][s]
						 + protect[kf][km][PREP_POOR][3] * p.infect[1][0][1][s]));

	beta[3][s] = 1.0 - exp(p.acts[kf][km] * (  protect[kf][km][NONE     ][0] * p.infect[0][1][0][s]
						 + protect[kf][km][NONE     ][1] * p.infect[1][1][0][s]));
	beta[4][s] = 1.0 - exp(p.acts[kf][km] * (  protect[kf][km][PREP_HIGH][0] * p.infect[0][1][0][s]
						 + protect[kf][km][PREP_HIGH][1] * p.infect[1][1][0][s]
                                                 + protect[kf][km][PREP_HIGH][2] * p.infect[0][1][1][s]
						 + protect[kf][km][PREP_HIGH][3] * p.infect[1][1][1][s]));
	beta[5][s] = 1.0 - exp(p.acts[kf][km] * (  protect[kf][km][PREP_POOR][0] * p.infect[0][1][0][s]
						 + protect[kf][km][PREP_POOR][1] * p.infect[1][1][0][s]
                                                 + protect[kf][km][PREP_POOR][2] * p.infect[0][1][1][s]
						 + protect[kf][km][PREP_POOR][3] * p.infect[1][1][1][s]));
      }

      for (int b(0); b < BANDS; ++b) {
	for (int q(0); q < PREP_STATES; ++q) {
	  std::fill(suminf[F ][kf][q][km][b], suminf[F ][kf][q][km][b] + TRANSMIT, 0.0);
	  std::fill(suminf[MU][km][q][kf][b], suminf[MU][km][q][kf][b] + TRANSMIT, 0.0);
	  std::fill(suminf[MC][km][q][kf][b], suminf[MC][km][q][kf][b] + TRANSMIT, 0.0);
	}

	wf = STATES * (BANDS * (F * LEVELS + kf) + b) + bands;
	wm = STATES * (BANDS * (M * LEVELS + km) + b) + bands;

	for (int s = Y1NWT; s < STATES; ++s) {
	  int v(attr_transmit[s]);
 	  suminf[F ][kf][NONE][km][b][v] += beta[0][s] * p.prop_risk_decr[s] * wm[s];
 	  suminf[MU][km][NONE][kf][b][v] += beta[0][s] * p.prop_risk_decr[s] * wf[s];
 	  suminf[MC][km][NONE][kf][b][v] += beta[3][s] * p.prop_risk_decr[s] * wf[s];

 	  suminf[F ][kf][PREP_HIGH][km][b][v] += beta[1][s] * p.prop_risk_decr[s] * wm[s];
 	  suminf[MU][km][PREP_HIGH][kf][b][v] += beta[1][s] * p.prop_risk_decr[s] * wf[s];
 	  suminf[MC][km][PREP_HIGH][kf][b][v] += beta[4][s] * p.prop_risk_decr[s] * wf[s];

 	  suminf[F ][kf][PREP_POOR][km][b][v] += beta[2][s] * p.prop_risk_decr[s] * wm[s];
 	  suminf[MU][km][PREP_POOR][kf][b][v] += beta[2][s] * p.prop_risk_decr[s] * wf[s];
 	  suminf[MC][km][PREP_POOR][kf][b][v] += beta[5][s] * p.prop_risk_decr[s] * wf[s];
	}
      }
    }
  }

  // kf, km: female and male partners' activity levels
  // bf, bm: female and male partners' age bands
  for (int kf(0); kf < LEVELS; ++kf) {
    for (int km(0); km < LEVELS; ++km) {
      for (int bf(0); bf < BANDS - 1; ++bf) {   // exclude ages 55+
	for (int bm(0); bm < BANDS - 1; ++bm) { // exclude ages 55+

	  // Calculate mixing terms
	  double wk, wb, wd;

	  // Men
	  wb = (1 - p.assort_age) * supply_age[F][bf] / supply_sex[F] + p.assort_age * (bm == bf);
	  wk = (1 - p.assort_act) * supply[F][bf][kf] / (eps + supply_age[F][bf]) + p.assort_act * (km == kf);
	  if (bm == bf && bm > 0) {
	    mm = (1 - p.assort_dif) * wb * wk;
	  } else if (bm == bf + 1) {
	    wd = (1 - p.assort_age) * supply_age[F][bf+1] / supply_sex[F] + p.assort_age;
	    mm = (wb + p.assort_dif * wd) * wk;
	  } else {
	    mm = wb * wk;
	  }

	  // Women
	  wb = (1 - p.assort_age) * supply_age[M][bm] / supply_sex[M] + p.assort_age * (bm == bf);
	  wk = (1 - p.assort_act) * supply[M][bm][km] / (eps + supply_age[M][bm]) + p.assort_act * (km == kf);
	  if (bm == bf && bf + 1 < BANDS - 1) {
	    mf = (1 - p.assort_dif) * wb * wk;
	  } else if (bm == bf + 1) {
	    wd = (1 - p.assort_age) * supply_age[M][bm-1] / supply_sex[M] + p.assort_age;
	    mf = (wb + p.assort_dif * wd) * wk;
	  } else {
	    mf = wb * wk;
	  }

	  // this assumes that men and women compensate for imbalance
	  // equally, in which case the balancing term simplifies to
	  // the square root of the supply ratio, and the contact rate
	  // simplifies to the value below
	  numer = mf * mm;
	  denom = eps + supply[F][bf][kf] * supply[M][bm][km];
	  contact = std::sqrt(numer / denom) * p.partner[F][bf][kf] * p.partner[M][bm][km];

	  for (int v(0); v < TRANSMIT; ++v) {
	    for (int q(0); q < PREP_STATES; ++q) {
	      // force acting on women
	      lambda[ F][bf][kf][q][v] += contact * suminf[ F][kf][q][km][bm][v];

	      // force acting on uncircumcised men
	      lambda[MU][bm][km][q][v] += contact * suminf[MU][km][q][kf][bf][v];

	      // force acting on circumcised men
	      lambda[MC][bm][km][q][v] += contact * suminf[MC][km][q][kf][bf][v];
	    }
	  }
	}
      }
    }
  }
}

double prep_cover_ft(double const t, Parameters const& p) {
  if (t < p.prep_time_rollout) {
    return 0.0;
  } else {
    return std::min(1.0, (t - p.prep_time_rollout) / p.prep_span_rollout);
  }
}

double prep_cover_dt(double const t, Parameters const& p) {
  const double x((t - p.prep_time_rollout) / p.prep_span_rollout);
  return ((x > 0) && (x <= 1.0)) / p.prep_span_rollout;
}

void rates_prep(double const t,
		double const y[],
		double const dy[],
		Parameters const& p,
		double uptake[][BANDS][LEVELS]) {

  for (int g(0); g < SEXES; ++g) {
    for (int b(0); b < BANDS; ++b) {
      std::fill(uptake[g][b], uptake[g][b] + LEVELS, 0.0);
    }
  }

  if (t < p.prep_time_rollout || t >= p.prep_time_closure) return;

  // maps sex and circumcision status to sex (F->F, MU->M, MC->M)
  const Sex sex[] = {F, M, M};

  const double prop_ft(prep_cover_ft(t,p));
  const double prop_dt(prep_cover_dt(t,p));
  double goal_ft, goal_dt;

  int g, a, b, k;
  double nX[SEXES][BANDS][LEVELS]; // number susceptible and not on PrEP
  double nP[SEXES][BANDS][LEVELS]; // number susceptible and on PrEP
  double dX[SEXES][BANDS][LEVELS]; // time-derivative of nX
  double dP[SEXES][BANDS][LEVELS]; // time-derivative of nY
  double nN, dN;
  double nNall, dNall, nPall, dPall;
  double prop, diff, scale;
  double nXp, nXe, nXu; // number prioritized / eligible / ineligble for PrEP
  double nPp, nPe, nPu;
  double const *z;
  double const *d;

  // number covered when coverage target in prioritized group met
  double priority_max(0.0);

  for (g = 0; g < SEXES; ++g) {
    for (b = 0; b < BANDS; ++b) {
      for (k = 0; k < LEVELS; ++k) {
	std::fill(nX[g][b], nX[g][b] + LEVELS, 0.0);
	std::fill(nP[g][b], nP[g][b] + LEVELS, 0.0);
	std::fill(dX[g][b], dX[g][b] + LEVELS, 0.0);
	std::fill(dP[g][b], dP[g][b] + LEVELS, 0.0);
      }
    }
  }

  for (int u(0); u < TRACKS; ++u) {
    for (g = 0; g < SEXCIRC; ++g) {
      for (a = 0; a < AGES - 1; ++a) {
	b = band[a];
	for (k = 0; k < LEVELS; ++k) {
	  z = offset(g, k, a, u) + y;
	  d = offset(g, k, a, u) + dy;
	  nX[sex[g]][b][k] += z[XN];
	  nP[sex[g]][b][k] += z[XH] + z[XP];
	  dX[sex[g]][b][k] += d[XN];
	  dP[sex[g]][b][k] += d[XH] + d[XP];
	}
      }
    }
  }

  nXp = nXe = nXu = 0.0;
  nPp = nPe = nPu = 0.0;  
  nNall = dNall = nPall = dPall = 0.0;
  for (g = 0; g < SEXES; ++g) {
    for (b = 0; b < BANDS - 1; ++b) {
      for (k = 0; k < LEVELS; ++k) {
	nXp += (p.prep_init[g][b][k] >  0.0) * nX[g][b][k];
	nXe += (p.prep_init[g][b][k] == 0.0) * nX[g][b][k];
	nXu += (p.prep_init[g][b][k] <  0.0) * nX[g][b][k];

	nPp += (p.prep_init[g][b][k] >  0.0) * nP[g][b][k];
	nPe += (p.prep_init[g][b][k] == 0.0) * nP[g][b][k];
	nPu += (p.prep_init[g][b][k] <  0.0) * nP[g][b][k];

	dNall += dX[g][b][k] + dP[g][b][k];
	dPall += dP[g][b][k];

	nN = nX[g][b][k] + nP[g][b][k];
	priority_max += (p.prep_init[g][b][k] > 0) * p.prep_init[g][b][k] * nN;
      }
    }
  }
  nPall = nPp + nPe + nPu;
  nNall = nXp + nXe + nXu + nPall;

  // scale reduces the level of coverage in prioritized groups to
  // satisfy the population-level coverage constraint, if necessary
  scale = std::min(1.0, (p.prep_target * nNall - nPe - nPu)  / priority_max);

  // Calculate uptake rates in prioritized populations
  for (g = 0; g < SEXES; ++g) {
    for (b = 0; b < BANDS - 1; ++b) {
      for (k = 0; k < LEVELS; ++k) {
	if (p.prep_init[g][b][k] > 0.0) {
	  goal_ft = scale * p.prep_init[g][b][k] * prop_ft;
	  goal_dt = scale * p.prep_init[g][b][k] * prop_dt;
	  nN = nX[g][b][k] + nP[g][b][k];
	  dN = dX[g][b][k] + dP[g][b][k];
	  prop = goal_ft * nN - nP[g][b][k];
	  diff = goal_dt * nN + goal_ft * dN - dP[g][b][k];
	  uptake[g][b][k] = std::max(0.0, (prop + diff) / nX[g][b][k]);

	  // update change in number on PrEP
	  dPall += uptake[g][b][k] * nX[g][b][k];
	}
      }
    }
  }

  // Calculate uptake rates in unprioritized, eligible populations
  goal_ft = prop_ft * p.prep_target;
  goal_dt = prop_dt * p.prep_target;
  prop = goal_ft * nNall - nPall;
  diff = goal_dt * nNall + goal_ft * dNall - dPall;
  for (g = 0; g < SEXES; ++g) {
    for (b = 0; b < BANDS - 1; ++b) {
      for (k = 0; k < LEVELS; ++k) {
	if (p.prep_init[g][b][k] == 0.0) {	
	  uptake[g][b][k] = std::max(0.0, (prop + diff) / nXe);
	}
      }
    }
  }

}

void rates_art(double const t, Parameters const& p, double uptake[]) {
  std::fill(uptake, uptake + STAGES, 0.0);
  double tbgn, tmid, tlim;
  double pbgn, plim, prop;
  for (size_t h(L2); h < STAGES; ++h) {
    tbgn = p.art_time_start[h];
    tlim = p.art_time_limit[h];
    tmid = 0.5 * (tbgn + tlim);
    pbgn = p.art_prop_start[h];
    plim = p.art_prop_limit[h];
    if (t < tbgn) {
      prop = 0.0;
    } else if (t < tmid) {
      prop = pbgn * (t - tbgn) / (tmid - tbgn);
    } else if (t < tlim) {
      prop = pbgn + (plim - pbgn) * (t - tmid) / (tlim - tmid);
    } else {
      prop = plim;
    }
    uptake[h] = fabs(log(1.0 - prop));
  }
}

// returns (by reference) the target circumcision coverage ft and the
// change in target coverage dt at time t
void mmc_cover(double const t, Parameters const& p, double& ft, double& dt) {
  const double t1(p.mmc_time_scaleup);
  const double t2(2012);
  const double t3(p.mmc_time_scaleup + p.mmc_span_scaleup);

  const double p1(p.mmc_prop_init);
  const double p2(0.232);
  const double p3(p.mmc_prop_late);

  const double dt1(t2 - t1), dt2(t3 - t2);
  const double dp1(p2 - p1), dp2(p3 - p2);

  if (t < t1) {
    ft = p1;
    dt = 0.0;
  } else if (t < t2) {
    ft = p1 + dp1 * (t - t1) / dt1;
    dt = dp1 / dt1;
  } else if (t < t3) {
    ft = p2 + dp2 * (t - t2) / dt2;
    dt = dp2 / dt2;
  } else {
    ft = p3;
    dt = 0.0;
  }
}

void rates_mmc(double const t,
	       double const y[],
	       double const dy[],
	       Parameters const& p,
	       double uptake[]) {
  double goal_ft, goal_dt;
  mmc_cover(t, p, goal_ft, goal_dt);

  double nnum, nmmc;
  double dnum, dmmc;
  double const* z;
  double const* d;

  nnum = dnum = nmmc = dmmc = 0.0;
  for (int u(0); u < TRACKS; ++u) {
    for (int k(0); k < LEVELS; ++k) {
      for (int a(0); a < AGES - 1; ++a) {
	z = offset(MU, k, a, u) + y;
	d = offset(MU, k, a, u) + dy;
	nnum += std::accumulate(z, z + STATES, 0.0);
	dnum += std::accumulate(d, d + STATES, 0.0);

	z = offset(MC, k, a, u) + y;
	d = offset(MC, k, a, u) + dy;
	nmmc += std::accumulate(z, z + STATES, 0.0);
	dmmc += std::accumulate(d, d + STATES, 0.0);
      }
    }
  }
  nnum += nmmc;
  dnum += dmmc;

  // no uptake at ages 55+
  uptake[BANDS-1] = 0.0;
  for (int b(0); b < BANDS - 1; ++b) {
    uptake[b] = ((goal_ft * nnum - nmmc) + (goal_dt * nnum + goal_ft * dnum - dmmc)) / (nnum - nmmc);
    if (uptake[b] < 0.0) uptake[b] = 0.0;
  }
}

void rates_csw_prop(double const t,
		    double const y[],
		    double const dy[],
		    Parameters const& p,
		    double csw_init[]) {
  double nnum[LEVELS], dnum[LEVELS], mnum[LEVELS];
  double nsum(0.0), dsum(0.0), numer, denom, target;
  double const* z;
  double const* d;
  for (int k(0); k < LEVELS; ++k) {
    nnum[k] = dnum[k] = mnum[k] = 0.0;
    for (int u(0); u < TRACKS; ++u) {
      for (int a(0); a < AGES - 1; ++a) {
	z = offset(F, k, a, u) + y;
	d = offset(F, k, a, u) + dy;
	nnum[k] += std::accumulate(z, z + STATES, 0.0);
	dnum[k] += std::accumulate(d, d + STATES, 0.0);
	mnum[k] += std::inner_product(z, z + STATES, p.prop_risk_init, 0.0);
      }
    }
    nsum += nnum[k];
    dsum += dnum[k];
  }
  target = p.prop_risk[F][HIGH] / p.prop_sex[F];
  numer = target * nsum - nnum[HIGH] + target * dsum - dnum[HIGH];
  //  numer = target * dsum - dnum[HIGH];
  denom = std::accumulate(mnum, mnum + HIGH, 0.0);
  std::fill(csw_init, csw_init + BANDS, max(numer / denom, 0.0));
}

void rates_csw_need(double const t,
		    double const y[],
		    double const dy[],
		    Parameters const& p,
		    double csw_init[]) {

  int g, k, a, b;
  int bm, bw, km, kw(HIGH);
  double wk, wb, wd;
  double const *z;
  double const *d;

  // maps sex and circumcision status to sex (F->F, MU->M, MC->M)
  const Sex sex[] = {F, M, M};

  double active[SEXES][BANDS][LEVELS];
  double supply[SEXES][BANDS][LEVELS];
  double supply_act[SEXES][LEVELS];
  double supply_age[SEXES][BANDS];
  double supply_sex[SEXES];

  // time-derivative of number sexually active
  double ddtact[SEXES][BANDS][LEVELS];

  // Mixing coefficients for highly active women (commercial sex workers)
  // mix_w[bw][bm][km]: preference of age bw highly active women for
  //   age bm men in activity level km
  // mix_m[bw][bm][km]: preference of age bm men in activity level km
  //    for age bw highly active women
  double mix_w[BANDS][BANDS][LEVELS];
  double mix_m[BANDS][BANDS][LEVELS];

  // The deficit is the number of sex workers required to accomodate
  // unmet demand in each age band
  double deficit[BANDS];

  std::fill(deficit, deficit + BANDS, 0.0);

  const double eps(std::numeric_limits<double>::denorm_min());

  for (g = 0; g < SEXES; ++g) {
    for (b = 0; b < BANDS; ++b) {
      std::fill(active[g][b], active[g][b] + LEVELS, 0.0);
      std::fill(ddtact[g][b], ddtact[g][b] + LEVELS, 0.0);
    }
  }

  for (int u(0); u < TRACKS; ++u) {
    for (a = 0; a < AGES; ++a) {
      b = band[a];
      for (g = 0; g < SEXCIRC; ++g) {
	for (k = 0; k < LEVELS; ++k) {
	  z = offset(g, k, a, u) + y;
	  d = offset(g, k, a, u) + dy;
	  active[sex[g]][b][k] += std::accumulate(z, z + STATES, 0.0);
	  ddtact[sex[g]][b][k] += std::accumulate(d, d + STATES, 0.0);
	}
      }
    }
  }

  for (g = 0; g < SEXES; ++g) {
    for (k = 0; k < LEVELS; ++k) {
      for (b = 0; b < BANDS; ++b) {
	supply[g][b][k] = active[g][b][k] * p.partner[g][b][k];
      }
    }
  }

  for (g = 0; g < SEXES; ++g) {
    supply_sex[g] = 0.0;
    for (b = 0; b < BANDS; ++b) {
      supply_age[g][b] = 0.0;
      for (k = 0; k < LEVELS; ++k) supply_age[g][b] += supply[g][b][k];
    }
    for (k = 0; k < LEVELS; ++k) {
      supply_act[g][k] = 0.0;
      for (b = 0; b < BANDS; ++b) supply_act[g][k] += supply[g][b][k];
      supply_sex[g] += supply_act[g][k];
    }
  }

  for (bw = 0; bw < BANDS - 1; ++bw) {
    for (bm = 0; bm < BANDS - 1; ++bm) {
      for (km = 0; km < LEVELS; ++km) {
	// men -> high activity women
	wb = (1 - p.assort_age) * supply_age[F][bw] / supply_sex[F] + p.assort_age * (bm == bw);
	wk = (1 - p.assort_act) * supply[F][bw][kw] / (eps + supply_age[F][bw]) + p.assort_act * (km == kw);
	if (bm == bw && bm > 0) {
	  mix_m[bw][bm][km] = (1 - p.assort_dif) * wb * wk;
	} else if (bm == bw + 1) {
	  wd = (1 - p.assort_age) * supply_age[F][bw+1] / supply_sex[F] + p.assort_age;
	  mix_m[bw][bm][km] = (wb + p.assort_dif * wd) * wk;
	} else {
	  mix_m[bw][bm][km] = wb * wk;
	}

	// high activity women -> men
	wb = (1 - p.assort_age) * supply_age[M][bm] / supply_sex[M] + p.assort_age * (bm == bw);
	wk = (1 - p.assort_act) * supply[M][bm][km] / (eps + supply_age[M][bm]) + p.assort_act * (km == kw);
	if (bm == bw && bw + 1 < BANDS - 1) {
	  mix_w[bw][bm][km] = (1 - p.assort_dif) * wb * wk;
	} else if (bm == bw + 1) {
	  wd = (1 - p.assort_age) * supply_age[M][bm-1] / supply_sex[M] + p.assort_age;
	  mix_w[bw][bm][km] = (wb + p.assort_dif * wd) * wk;
	} else {
	  mix_w[bw][bm][km] = wb * wk;
	}
      }
    }
  }

  // Pool of women in lower activity levels who may become sex
  // workers, weighted by propensity to increase sexual activity
  double pool[BANDS];
  std::fill(pool, pool + BANDS, 0.0);
  for (int u(0); u < TRACKS; ++u) {
    for (a = 0; a < AGES - 1; ++a) {
      b = band[a];
      for (k = 0; k < HIGH; ++k) {
	z = offset(F, k, a, u) + y;
	pool[b] += std::inner_product(z, z + STATES, p.prop_risk_init, 0.0);
      }
    }
  }

  double numer, denom;
  std::fill(csw_init, csw_init + BANDS, 0.0);
  for (bw = 0; bw < BANDS - 1; ++bw) {
    numer = 0.0;
    denom = eps;
    for (bm = 0; bm < BANDS - 1; ++bm) {
      for (km = 0; km < LEVELS; ++km) {
	numer += mix_w[bw][bm][km] * mix_m[bw][bm][km] * supply[M][bm][km];
	denom += mix_w[bw][bm][km] * mix_w[bw][bm][km];
      }
    }
    denom *= p.partner[F][bw][HIGH];
    csw_init[bw] = (numer / denom - active[F][bw][HIGH]) / pool[bw];

    // The upper limit of 1 was introduced because some
    // parameterizations considered during calibration may lead to
    // extreme imbalances that induce high sex work initiation rates
    // and destabilize numerical integration. One was chosen because
    // then it takes about one year on average to meet demand
    csw_init[bw] = std::min(std::max(0.0, csw_init[bw]), 1.0);
  }
}

void rates_replace(double const t,
		   double const y[],
		   double const dy[],
		   Parameters const& p,
		   double rate[][LEVELS]) {

  // nnum - number of individuals in each sex and activity level 
  // dnum - time-derivative of nnum
  // mnum - number of individuals, scaled by propensity to increase risk behavior
  double nnum[SEXES][LEVELS], dnum[SEXES][LEVELS], mnum[SEXES][LEVELS];
  double dsum[SEXES];
  double const* z;
  double const* d;

  // maps sex and circumcision status to sex
  // F->F, MU->M, MC->M
  const Sex sex[] = {F, M, M};

  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      nnum[g][k] = dnum[g][k] = mnum[g][k] = 0.0;
    }
  }

  // Sum population sizes and their derivatives
  for (int u(0); u < TRACKS; ++u) {
    for (int g(0); g < SEXCIRC; ++g) {
      for (int k(0); k < LEVELS; ++k) {
	for (int a(0); a < AGES - 1; ++a) {
	  z = offset(g, k, a, u) + y;
	  d = offset(g, k, a, u) + dy;
	  nnum[sex[g]][k] += std::accumulate(z, z + STATES, 0.0);
	  dnum[sex[g]][k] += std::accumulate(d, d + STATES, 0.0);
	  mnum[sex[g]][k] += std::inner_product(z, z + STATES, p.prop_risk_init, 0.0);
	}
      }
    }
  }

  for (int g(0); g < SEXES; ++g) {
    dsum[g] = std::accumulate(dnum[g], dnum[g] + LEVELS, 0.0);
    rate[g][HIGH] = 0.0;
    for (int k(HIGH); k > 0; --k) {
      const double target(p.prop_risk[g][k] / p.prop_sex[g]);
      rate[g][k-1] = (dsum[g] * target - dnum[g][k] + rate[g][k] * mnum[g][k]) / mnum[g][k-1];
      if (rate[g][k-1] < 0.0) rate[g][k-1] = 0.0;
    }
  }
}

void initialize(double y[], Parameters const& p) {
  std::fill(y, y + COMPARTMENTS, 0.0);
  // Simulation starts with no HIV. The seed cases must be manually
  // entered into the population at the epidemic start date.

  // Age band "width" ignores the 55+ age band
  const int span((AGES - 1) / (BANDS - 1));
  double num;
  double *z;

  // maps sex and circumcision status to sex
  // F->F, MU->M, MC->M
  const Sex sex[] = {F, M, M};

  for (int g(0); g < SEXCIRC; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      // We ignore the initial population aged 55+, since we assume
      // negligible transmission in this group. They are included in
      // the model to capture cohort outcomes from 2015, at which
      // point the majority of people aged 55+ at 1978 will have died
      for (int b(0); b < BANDS - 1; ++b) {
  	num = p.init_size * p.init_risk[sex[g]][b][k] / static_cast<double>(span);
      	for (int a(0); a < span; ++a) {
  	  z = offset(g, k, b * span + a, ANTE) + y;
  	  if (g == F) {
  	    z[XN] = num;
  	  } else if (g == MU) {
  	    z[XN] = num * (1.0 - p.mmc_prop_init);
	  } else {
  	    z[XN] = num * p.mmc_prop_init;
  	  }
  	}
      }
    }
  }
}
