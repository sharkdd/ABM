#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <numeric>
#include <vector>
#include <gsl/gsl_math.h>
#include <ModelTypes.H>
#include <Person.H>

double infinity = std::numeric_limits<double>::infinity();
Parameters params;

Parameters::Parameters() {
  // This only automatically assigns parameters that have default
  // values that are not application-specific. 

  // We allocate storage for all possible values of enum values that
  // stratify the population. This sometimes means that we have unused
  // storage for absurd parameters (like progression while
  // susceptible). This constructor initializes this extra storage to
  // safe values.

  // Susceptible individuals don't experience disease progression
  wait_stage[SUSCEPTIBLE] = infinity;

  // Initialize the times to emergence or reversion to infinity so
  // that the programmer only has to initialize parameters for the
  // combinations of virus and ARV state that are allowed
  for (unsigned int vi(0); vi < VARIANTS; ++vi) {
    for (unsigned int ti(0); ti < ART_STATES; ++ti) {
      for (unsigned int qi(0); qi < PREP_STATES; ++qi) {
	wait_emerge[vi][ti][qi] = infinity;
	wait_revert[vi][ti][qi] = infinity;
      }
    }
  }

  // The base transmission probability assumes wild-type infection
  prop_transmit[WT ] = 1.0;
  prop_transmit[WR1] = 1.0;
  prop_transmit[WC1] = 1.0;
  prop_transmit[WQ1] = 1.0;
  prop_transmit[WR2] = 1.0;
  prop_transmit[WC2] = 1.0;
  prop_transmit[WQ2] = 1.0;

  // The base progression rates assume wild-type infection
  prop_progress[WT ] = 1.0;
  prop_progress[WR1] = 1.0;
  prop_progress[WC1] = 1.0;
  prop_progress[WQ1] = 1.0;
  prop_progress[WR2] = 1.0;
  prop_progress[WQ2] = 1.0;
  prop_progress[WC2] = 1.0;

  // establish safe default ART availability assumptions so that users
  // only need to specify relevant values.
  for (unsigned int si(0); si < STAGES; ++si) {
    art.year_available[si] = infinity;
  }

  // default to no ART uptake
  std::fill(art.rate_uptake_max, art.rate_uptake_max + STAGES, 0.0);

  // no ART is the baseline for transmission reduction on ART
  art.prop_transmit[ART_NONE] = 1.0;

  // baseline for ART-related mortality is AIDS
  art.prop_mortality[AIDS] = 1.0;

  // several ART-related parameters are only appropriate when actively
  // receiving or adherent to ART
  art.rate_dropout[ART_NONE] = 0.0;
  art.rate_nonadhere[ART_NONE] = 0.0;
  art.rate_nonadhere[ART_FAIL] = 0.0;
  art.prop_mortality[SUSCEPTIBLE] = 0.0;
  art.prop_progress[ART_NONE] = 1.0;

  // default to no PrEP uptake
  prep.prop_debut = Parameters::zero_prep;
  prep.rate_uptake = Parameters::zero_prep;
  prep.rate_uptake_max = 0.0;
}

int Parameters::read_epidemic(std::string const& filename) {
  std::ifstream in(filename.c_str());
  if (in.fail()) return -1;

  std::istream_iterator<double> start(in), end;
  std::vector<double> input(start, end);
  in.close();

  if (input.size() != 119) return -2;

#warning "Hard-coded constant parameters"
  const double base_year(1978);
  span_age_band = 5;
  std::fill(damp_concurrency[ MALE ], damp_concurrency[ MALE ] + LEVELS, 0.0); // default: no concurrency
  std::fill(damp_concurrency[FEMALE], damp_concurrency[FEMALE] + LEVELS, 0.0); // default: no concurrency

  // Partnership dissolution rate. This will eventually be stratified
  // by partnership type (i.e., partner sexual activity
  // levels). Currently homogeneous for ease of debugging.
  for (int ki(0); ki < LEVELS; ++ki) {
    for (int kj(0); kj < LEVELS; ++kj) {
      rate_breakup[ki][kj] = 12.0;
    }
  }

  // The simulator is based on year zero (year zero = base_year). The
  // input in year_init is the year that HIV is introduced into the
  // population.
  year_init = round(input[0]) - base_year;

  size_population = input[1];

  rate_grow_init = input[2];
  rate_grow_late = input[3];
  rate_grow_decr = input[4];

  prop_female = 1.0 - input[5];

  std::vector<double> prop_level(LEVELS);
  std::vector<double> prop_aband(AGES);

  prop_level[LEAST] = input[6];
  prop_level[LOW  ] = 1 - (input[6] + input[7] + input[8]);
  prop_level[MID  ] = input[7];
  prop_level[HIGH ] = input[8];
  prop_activity(FEMALE, prop_level);

  prop_level[LEAST] = input[9];
  prop_level[LOW  ] = 1 - (input[9] + input[10] + input[11]);
  prop_level[MID  ] = input[10];
  prop_level[HIGH ] = input[11];
  prop_activity(MALE, prop_level);

  for (int b(0); b < AGES; ++b) {
    prop_aband[b] = input[12 + b];
  }
  prop_age(prop_aband);

  rate_csw_exit = 1.0 / input[21];

  for (int h(0); h < STAGES; ++h) {
    prop_risk_exit[h][ART_NONE ] = 1.00;
    prop_risk_exit[h][ART_EARLY] = 2.50;
    prop_risk_exit[h][ART_LATE ] = 2.50;
    prop_risk_exit[h][ART_FAIL ] = 2.50;

    prop_risk_init[h][ART_NONE ] = 1.00;
    prop_risk_init[h][ART_EARLY] = 0.40;
    prop_risk_init[h][ART_LATE ] = 0.40;
    prop_risk_init[h][ART_FAIL ] = 0.40;

    prop_seek_decr[h][ART_NONE ] = 1.00;
    prop_seek_decr[h][ART_EARLY] = 1.00;
    prop_seek_decr[h][ART_LATE ] = 1.00;
    prop_seek_decr[h][ART_FAIL ] = 1.00;
  }
  prop_risk_exit[CHRONIC_LATE][ART_NONE] = 2.50;
  prop_risk_exit[AIDS        ][ART_NONE] = 4.00;

  prop_risk_init[CHRONIC_LATE][ART_NONE] = 0.40;
  prop_risk_init[AIDS        ][ART_NONE] = 0.25;

  prop_seek_decr[AIDS][ART_NONE] = 0.65;

#warning "Tracking of individuals aged 55+ not supported"
  // #warning "Hard-coded survival at ages 55+"
  for (int b(0); b < AGES; ++b) wait_death[FEMALE][b] = 1.0 / input[22 + b];
  for (int b(0); b < AGES; ++b) wait_death[ MALE ][b] = 1.0 / input[30 + b];
  //   rate_mort[F][AGES-1] = 0.04694; // reciprocal of WHO remaining life expectancy at age 55
  //   rate_mort[M][AGES-1] = 0.05828; // reciprocal of WHO remaining life expectancy at age 55

#warning "Partner change rates for individuals aged 55+ not initialized"
  for (int k(0); k < LEVELS; ++k) {
    for (int b(0); b < AGES; ++b) {
      rate_partner[FEMALE][k][b] = input[38] * input[40+b] * input[56+k];
      rate_partner[ MALE ][k][b] = input[39] * input[48+b] * input[60+k];
    }
  }

  for (int kf(0); kf < LEVELS; ++kf) {
    for (int km(0); km < LEVELS; ++km) {
      num_acts[kf][km] = input[64 + kf * LEVELS + km];
    }
  }

  for (int kf(0); kf < LEVELS; ++kf) {
    for (int km(0); km < LEVELS; ++km) {
      prop_condom[kf][km] = input[80 + kf * LEVELS + km];
    }
  }

  prop_condom_efficacy = 1.0 - input[96];
  prop_condom_nonuse = input[97];
  time_condom_change = input[98] - base_year;
  span_condom_change = input[99];

  // The parameters of MMC scale-up are set in read_intervention. The
  // values below assume MMC stabilizes at 2012 HSRC levels for KZN.
  mmc.prop_base = input[100];
  mmc.efficacy = 1.0 - input[101];

  mmc.target_max = 0.232;
  mmc.year_init = 2011 - base_year;
  mmc.year_done = 2012 - base_year;
  mmc.change_max = (mmc.target_max - mmc.prop_base) / (mmc.year_done - mmc.year_init);

  assort_act = input[102];
  assort_age = input[103];
  assort_dif = input[104];

  // NOTE: input[108] not used (reserved for CD4>500 vs. 500-350)
  prob_transmit[ACUTE_WINDOW ] = input[105] * input[106];
  prob_transmit[ACUTE_DETECT ] = input[105] * input[107];
  prob_transmit[CHRONIC_EARLY] = input[105] * input[109];
  prob_transmit[CHRONIC_LATE ] = input[105] * input[110];
  prob_transmit[AIDS         ] = input[105] * input[111];

  // Infectivity of non-WT variants is set by read_intervention.
  std::fill(prop_transmit, prop_transmit + VARIANTS, 0.0);
  prop_transmit[WT] = 1.0;

  // EDIT POINT 2015-02-02
  wait_stage[ACUTE_WINDOW ] = input[112];
  wait_stage[ACUTE_DETECT ] = input[113];
  wait_stage[CHRONIC_EARLY] = input[114] + input[115];
  wait_stage[CHRONIC_LATE ] = input[116];
  wait_stage[AIDS         ] = input[117];

#warning "Sex worker replacement mode ignored"
  // Input 118 specifies the sex worker replacement algorithm used.
  // This input is currently ignored.

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

  const double base_year(1978);

  // +-+ PrEP parameters +-----------------------------------------------------+
  prep.year_available = input[0] - base_year;
  prep.span_scale = input[1];
  prep.year_closure = input[2] - base_year;
  prep.target_overall = input[3];
  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int a(0); a < AGES; ++a) {
	prep.target[g][a][k] = input[4 + (g * LEVELS + k) * AGES + a];
      }
    }
  }

  prep.wait_dose = 1.0 / input[68]; // dosing frequency
  prep.prop_test = input[70]; // proportion of doses coupled with HIV testing
  prep.wait_exit = input[71]; // duration of PrEP use
  prep.prop_high = input[72]; // proportion with high PrEP concentrations

  std::fill(prep.efficacy, prep.efficacy + VARIANTS, 1.0 - input[73]);
  prep.efficacy[C1] = 1.0 - input[76];
  prep.efficacy[Q1] = 1.0 - input[77];

  // "adherence" is a misnomer for TMC278LA, here used to encode that
  // individuals with PREP_HIGH receive full efficacy, PREP_POOR
  // receive no efficacy.
  prep.adherence[PREP_NONE] = std::numeric_limits<double>::quiet_NaN();
  prep.adherence[PREP_HIGH] = 1.0;
  prep.adherence[PREP_POOR] = 0.0;

  // waiting times until PrEP resistance emergence. The absolute value
  // is taken below because otherwise if PrEP does not cause
  // resistance to emerge (inputs 74 or 75 = 0) then the wait time
  // evaluates to -inf instead of +inf.
  const double wait_prep_emerge[] = {
    infinity,
    fabs(input[78] / -log(1.0 - input[74] * 0.99)), // PREP_HIGH
    fabs(input[78] / -log(1.0 - input[75] * 0.99))  // PREP_POOR
  };

  for (int qi(PREP_HIGH); qi < PREP_STATES; ++qi) {
    wait_emerge[WT ][ART_NONE][qi] = wait_prep_emerge[qi];
    wait_emerge[R1 ][ART_NONE][qi] = wait_prep_emerge[qi];
    wait_emerge[WR1][ART_NONE][qi] = wait_prep_emerge[qi];
    //    wait_emerge[S1 ][ART_NONE][qi] = wait_prep_emerge[qi];
    //    wait_emerge[WS1][ART_NONE][qi] = wait_prep_emerge[qi];
    wait_emerge[WQ1][ART_NONE][qi] = wait_prep_emerge[qi] * 2.0;
    wait_emerge[WC1][ART_NONE][qi] = wait_prep_emerge[qi] * 2.0;

    wait_emerge[R2 ][ART_NONE][qi] = std::numeric_limits<double>::quiet_NaN();
    wait_emerge[C2 ][ART_NONE][qi] = std::numeric_limits<double>::quiet_NaN();
    wait_emerge[WR2][ART_NONE][qi] = std::numeric_limits<double>::quiet_NaN();
    wait_emerge[WC2][ART_NONE][qi] = std::numeric_limits<double>::quiet_NaN();

    // Ensure that the model outputs obviously incorrect results if
    // individuals are ever taking ART and PrEP simultaneously
    for (int ti(ART_EARLY); ti < ART_STATES; ++ti) {
      wait_emerge[WT ][ti][qi] = std::numeric_limits<double>::quiet_NaN();
      wait_emerge[R1 ][ti][qi] = std::numeric_limits<double>::quiet_NaN();
      wait_emerge[WR1][ti][qi] = std::numeric_limits<double>::quiet_NaN();
      //      wait_emerge[S1 ][ti][qi] = std::numeric_limits<double>::quiet_NaN();
      //      wait_emerge[WS1][ti][qi] = std::numeric_limits<double>::quiet_NaN();
      wait_emerge[WQ1][ti][qi] = std::numeric_limits<double>::quiet_NaN();
      wait_emerge[WC1][ti][qi] = std::numeric_limits<double>::quiet_NaN();
    }
  }

  // +-+ ART parameters +------------------------------------------------------+

  // delay ART implementation indefinitely in ART-ineligible stages
  std::fill(art.year_available, art.year_available + STAGES, std::numeric_limits<double>::infinity());
  std::fill(art.year_scale, art.year_scale + STAGES, std::numeric_limits<double>::infinity());

  art.year_available[AIDS] = 2004 - base_year;
  art.year_scale[AIDS] = 2012 - base_year;
  art.prop_scale[AIDS] = 0.5555;
  art.prop_level[AIDS] = 0.5563;

  art.year_available[CHRONIC_LATE] = 2010 - base_year;
  art.year_scale[CHRONIC_LATE] = 2017 - base_year;
  art.prop_scale[CHRONIC_LATE] = 0.5426;
  art.prop_level[CHRONIC_LATE] = 0.3870;

  // Calculate the maximum uptake rate. This exploits exact knowledge
  // of the family of ART uptake functions used: maximum uptake occurs
  // either midway through scale-up or once scale-up completes.
  for (int h(ACUTE_WINDOW); h < STAGES; ++h) {
    art.rate_uptake_max[h] = 0.0;
    double tlim(art.year_scale[h]);
    double tmid(0.5 * (art.year_scale[h] + art.year_available[h]));
    art.rate_uptake_max[h] = std::max(art.uptake(tmid, Stage(h)), art.uptake(tlim, Stage(h)));
  }

  art.rate_advance = 1.0;

  art.rate_dropout[ART_EARLY] = input[79];
  art.rate_dropout[ART_LATE ] = input[80];
  art.rate_dropout[ART_FAIL ] = input[81];

  avfailE[WT] = input[82];
  avfailL[WT] = input[83];
  avfailE[R1] = avfailE[WR1] = avfailE[C1] = avfailE[WC1] = input[84];
  avfailL[R1] = avfailL[WR1] = avfailL[C1] = avfailL[WC1] = input[85];
  //  avfailE[S1] = avfailE[WS1] = input[86];
  //  avfailL[S1] = avfailL[WS1] = input[87];
  avfailE[Q1] = avfailE[WQ1] = input[88];
  avfailL[Q1] = avfailL[WQ1] = input[89];
  avfailE[Q2] = avfailE[WQ2] = input[90];
  avfailL[Q2] = avfailL[WQ2] = input[91];

  art.rate_nonadhere[ART_EARLY] = input[92] * -log(1.0 - avfailE[WT]);
  art.rate_nonadhere[ART_LATE ] = input[93] * -log(1.0 - avfailL[WT]);

  for (int v(0); v < VARIANTS; ++v) {
    if (v == R2 || v == C2) { // no emergence, resistance already present
      wait_emerge[v][ART_EARLY][PREP_NONE] = infinity;
      wait_emerge[v][ART_LATE ][PREP_NONE] = infinity;
    } else if (v == WR2 || v == WC2) { // immediate re-emergence
      wait_emerge[v][ART_EARLY][PREP_NONE] = 0.0;
      wait_emerge[v][ART_LATE ][PREP_NONE] = 0.0;
    } else {
      wait_emerge[v][ART_EARLY][PREP_NONE] = 1.0 / (-log(1.0 - avfailE[v]) - art.rate_nonadhere[ART_EARLY]);
      wait_emerge[v][ART_LATE ][PREP_NONE] = 1.0 / (-log(1.0 - avfailL[v]) - art.rate_nonadhere[ART_LATE ]);
    }
    for (int qi(PREP_HIGH); qi < PREP_STATES; ++qi) {
      wait_emerge[v][ART_EARLY][qi] = std::numeric_limits<double>::quiet_NaN();
      wait_emerge[v][ART_LATE ][qi] = std::numeric_limits<double>::quiet_NaN();
      wait_emerge[v][ART_FAIL ][qi] = std::numeric_limits<double>::quiet_NaN();
    }
  }

  art.prop_cross_resist = input[94];

  // It is an error to evaluate ART mortality for an individuals not
  // on ART. There is no treatment-associated mortality associated
  // with non-adherence (non-adherent individuals experience disease
  // progression as though not on ART)
  art.rate_mortality[ART_NONE ] = std::numeric_limits<double>::quiet_NaN();
  art.rate_mortality[ART_EARLY] = input[95];
  art.rate_mortality[ART_LATE ] = input[96];
  art.rate_mortality[ART_FAIL ] = 0.0;

  // set so that the ART mortality rate is undefined for individuals
  // who should never be on ART.
  std::fill(art.prop_mortality, art.prop_mortality + STAGES, std::numeric_limits<double>::quiet_NaN());
  art.prop_mortality[CHRONIC_LATE] = input[97];
  art.prop_mortality[AIDS        ] = 1.0;

#warning "Stratify art efficacy by variant then simplify rate_transmit"
  art.prop_transmit[ART_EARLY] = 1.0 - input[98];
  art.prop_transmit[ART_LATE ] = 1.0 - input[98];
  art.prop_transmit[ART_FAIL ] = 1.0;

#warning "Input 99 (ART efficacy against DR virus) not used"
  art.prop_progress[ART_NONE ] = 1.0;
  art.prop_progress[ART_EARLY] = 0.0;
  art.prop_progress[ART_LATE ] = 0.0;
  art.prop_progress[ART_FAIL ] = 1.0;

#warning "Second-line ART not modeled yet"
  //   sec_init_dr = input[100];
  //   sec_init_na = input[101];
  //   for (int v(0); v < VARIANTS; ++v) {
  //     emerge[ART2][v] = -input[104] * log(1.0 - input[102]);
  //   }
  //   emerge[ART2][ S1] = emerge[ART2][WS1] = input[105] * emerge[ART2][WT];

  //   sec_fail_dr = -log(1.0 - input[102]) - emerge[ART2][WT];
  //   sec_fail_na = -log(1.0 - input[103]) - emerge[ART2][WT];

  //   std::fill(sec_mort, sec_mort + STAGES, std::numeric_limits<double>::quiet_NaN());
  //   sec_mort[CHRONIC_LATE] = input[106] * input[97];
  //   sec_mort[AIDS        ] = input[106];

  //   sec_exit = input[107];

  // +-+ Drug resistance fitness parameters +----------------------------------+
  wait_revert[C1][ART_NONE][PREP_NONE] = input[108];
  wait_revert[C2][ART_NONE][PREP_NONE] = input[109];
  wait_revert[R1][ART_NONE][PREP_NONE] = input[110];
  wait_revert[R2][ART_NONE][PREP_NONE] = input[111];
  //  wait_revert[S1][ART_NONE][PREP_NONE] = input[112];
  //  wait_revert[S2][ART_NONE][PREP_NONE] = input[113];
  wait_revert[Q1][ART_NONE][PREP_NONE] = input[114];
  wait_revert[Q2][ART_NONE][PREP_NONE] = input[115];

  wait_revert[R1][ART_NONE][PREP_HIGH] = wait_revert[R1][ART_NONE][PREP_NONE];
  wait_revert[R1][ART_NONE][PREP_POOR] = wait_revert[R1][ART_NONE][PREP_NONE];

  //  wait_revert[S1][ART_NONE][PREP_HIGH] = wait_revert[S1][ART_NONE][PREP_NONE];
  //  wait_revert[S1][ART_NONE][PREP_POOR] = wait_revert[S1][ART_NONE][PREP_NONE];

#warning "Reversion on second-line ART not implemented"
  //   revert[ART2][R1] = revert[NONE][R1];
  //   revert[ART2][C1] = revert[NONE][C1];
  //   revert[ART2][Q1] = revert[NONE][Q1];
  //   revert[ART2][Q2] = revert[NONE][Q2];
  //   revert[ART2][R2] = revert[NONE][R2];
  //   revert[ART2][C2] = revert[NONE][C2];

  for (int v(0); v < VARIANTS; ++v) {
    wait_revert[v][ART_FAIL][PREP_NONE] = wait_revert[v][ART_NONE][PREP_NONE];
  }

  prop_progress[WT] = 1.0;
  prop_progress[C1] = input[116];
  prop_progress[C2] = input[117];
  prop_progress[R1] = input[118];
  prop_progress[R2] = input[119];
  //  prop_progress[S1] = input[120];
  //  prop_progress[S2] = input[121];
  prop_progress[Q1] = input[122];
  prop_progress[Q2] = input[123];
  prop_progress[WQ1] = prop_progress[WQ2] = 1.0;
  prop_progress[WC1] = prop_progress[WC2] = 1.0;
  prop_progress[WR1] = prop_progress[WR2] = 1.0;
  //  prop_progress[WS1] = prop_progress[WS2] = 1.0;

  std::fill(prop_transmit, prop_transmit + VARIANTS, 1.0);
  prop_transmit[C1] = input[124];
  prop_transmit[C2] = input[125];
  prop_transmit[R1] = input[126];
  prop_transmit[R2] = input[127];
  //  prop_transmit[S1] = input[128];
  //  prop_transmit[S2] = input[129];
  prop_transmit[Q1] = input[130];
  prop_transmit[Q2] = input[131];

  mmc.target_max = input[132];
  mmc.year_init = input[133] - base_year;
  mmc.year_done = input[134] - base_year;
  mmc.change_max = (input[132] - mmc.prop_base) / (mmc.year_done - mmc.year_init);

  return 0;
}

void Parameters::prop_activity(Sex sex, std::vector<double> const& proportions) {
  alias_setup(proportions,
	      m_activity_table[sex],
	      m_activity_alias[sex],
	      LEVELS);
}

double Parameters::prop_activity(Sex sex, Activity activity) const {
  const unsigned int act(activity); // cast to clean up warnings.
  double sum(m_activity_table[sex][act]);
  for (int ai(0); ai < LEVELS; ++ai)
    if (m_activity_alias[sex][ai] == act)
      sum += 1.0 - m_activity_table[sex][ai];
  return sum / static_cast<double>(LEVELS);
}

Activity Parameters::sample_activity_level(gsl_rng *rng, Sex sex) const {
  const Activity A(static_cast<Activity>(gsl_rng_uniform_int(rng, LEVELS)));
  const double U(gsl_rng_uniform(rng));
  return static_cast<Activity>((U <= m_activity_table[sex][A]) ? A : m_activity_alias[sex][A]);
}

void Parameters::prop_age(std::vector<double> const& proportions) {
  alias_setup(proportions, m_age_table, m_age_alias, AGES);
}

double Parameters::prop_age(unsigned int band) const {
  double sum(m_age_table[band]);
  for (int ai(0); ai < AGES; ++ai)
    if (m_age_alias[ai] == band)
      sum += 1.0 - m_age_table[ai];
  return sum / static_cast<double>(AGES);
}

std::pair<double, unsigned int> Parameters::sample_age(gsl_rng* rng) const {
  const unsigned int A(gsl_rng_uniform_int(rng, AGES));
  const double U(gsl_rng_uniform(rng));
  const double W(gsl_rng_uniform(rng));
  const unsigned int band((U <= m_age_table[A]) ? A : m_age_alias[A]);
  const double age((W + band) * span_age_band);
  return std::pair<double,unsigned int>(age, band);
}

// This recalculates the rate of transmission from scratch. This does
// involve several memory accesses and math operations, including
// taking powers, so it may be appropriate to cache as much of the
// calculation as possible if much time is spent here (particularly as
// the model expands)
//
// The transmission rate depends on time because of changes in condom use
double Parameters::rate_transmit(Person const* donor, Person const* recipient, double const time) const {
  const Activity kr(recipient->activity());
  const Activity kd(donor->activity());

  const double a(num_acts[kr][kd]);

  // calculate proportion of acts in which condoms are used
  const double cu(prop_condom[kr][kd]);

  double eta;
  eta = std::max(prop_condom_nonuse,
		 std::min(1.0, 1 - (1 - prop_condom_nonuse) * (time - time_condom_change) / span_condom_change));

  const double an(a * (1 - cu) * eta); // number of unprotected sex acts
  const double ac(a - an);             // number of protected sex acts

  // virus multiplier
  const double vm(prop_transmit[donor->virus()]);

  // ART multiplier, applies to per-act transmission probability. No
  // benefit from ART if the virus is ART-resistant
  const double tm((donor->virus() == R2 || donor->virus() == C2) ? 1.0 : art.prop_transmit[donor->treat()]);

  // PrEP multiplier, applies to transmission rate. A person is PrEP
  // protected when he or she is enrolled in PrEP and responds to PrEP
  const bool drug(recipient->prep() != PREP_NONE);
  const double pm(1.0 - drug * prep.efficacy[donor->virus()]);

  // circumcision multiplier. requires females circumcised() == false
  const double cm(1.0 - recipient->circumcised() * (1 - params.mmc.efficacy));

  // Per-act transmission probability
  const double pt(cm * vm * tm * prob_transmit[donor->stage()]);

  // Per-partnership transmission probability
  const double prob(1.0 - pow(1.0 - pt * prop_condom_efficacy, ac) * pow(1.0 - pt, an));

  // The Per-partnership transmission rate (pm * prob / (1 - pm *
  // prob)) is multiplied by the partnership dissolution rate to get
  // the transmission rate per unit time
  return params.rate_breakup[kr][kd] * pm * prob / (1 - pm * prob);
}

double Parameters::wait_progress(Person const* p) const {
  const double base(wait_stage[p->stage()]);

  // ART multiplier; ART has no effect on progression if the virus has
  // acquired resistance to ART
  const double tm((p->virus() == R2 || p->virus() == C2) ? 1.0 : 1.0 / art.prop_progress[p->treat()]);

  // virus multiplier
  const double vm(1.0 / prop_progress[p->virus()]);

  return vm * tm * base;
}


void Parameters::alias_setup(std::vector<double> const& weights,
			     double* table,
			     unsigned int* alias,
			     size_t n) {
  assert(weights.size() == n);
  const double sum(std::accumulate(weights.begin(), weights.begin() + n, 0.0));

  std::list<int> smaller, greater;
  for (size_t k(0); k < n; ++k) {
    table[k] = n * weights[k] / sum;
    if (table[k] < 1.0) {
      smaller.push_back(k);
    } else {
      greater.push_back(k);
    }
  }

  int s, g;
  while (!smaller.empty() && !greater.empty()) {
    s = smaller.front();
    g = greater.front();
    alias[s] = g;
    table[g] -= (1.0 - table[s]);
    if (table[g] < 1.0) {
      smaller.push_back(g);
      greater.pop_front();
    }
    smaller.pop_front();
  }

}

double Parameters::params_art::uptake(double const t, Stage const h) const {
  double tbgn, tmid, tlim;
  double pbgn, plim, prop;
  tbgn = year_available[h];
  tlim = year_scale[h];
  tmid = 0.5 * (tbgn + tlim);
  pbgn = prop_scale[h];
  plim = prop_level[h];
  if (t < tbgn) {
    prop = 0.0;
  } else if (t < tmid) {
    prop = pbgn * (t - tbgn) / (tmid - tbgn);
  } else if (t < tlim) {
    prop = pbgn + (plim - pbgn) * (t - tmid) / (tlim - tmid);
  } else {
    prop = plim;
  }
  return fabs(log(1.0 - prop));
}

std::pair<double,double> Parameters::params_mmc::coverage(double const t) const {
  // The hard-coded 0.232 value at 2012 is observed coverage in KZN at
  // 2012 from the 2012 HSRC survey
  const double t1(year_init), t2(2012-1978), t3(year_done);
  const double p1(prop_base), p2(0.232), p3(target_max);
  const double dt1(t2 - t1), dt2(t3 - t2);
  const double dp1(p2 - p1), dp2(p3 - p2);
  double target, change;
  if (t < t1) {
    target = p1;
    change = 0.0;
  } else if (t < t2) {
    target = p1 + dp1 * (t - t1) / dt1;
    change = dp1 / dt1;
  } else if (t < t3) {
    target = p2 + dp2 * (t - t2) / dt2;
    change = dp2 / dt2;
  } else {
    target = p3;
    change = 0.0;
  }
  return std::pair<double,double>(target, change);
}
