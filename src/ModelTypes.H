#ifndef MODELTYPES_H
#define MODELTYPES_H

#include <fcntl.h>
#include <cmath>
#include <limits>
#include <vector>
#include <tr1/unordered_set>
#include <gsl/gsl_rng.h>

// ModelTypes is used to define data types, constants and functions
// that are global in some sense to the entire individual-based model.

namespace std {
  namespace tr1 {
    // Defines a hash functor for std::pairs of hashable types
    template<typename X, typename Y> 
    struct hash< std::pair<X,Y> > {
      size_t operator()(std::pair<X,Y> const& z) const throw() {
	return hash<X>()(z.first) ^ hash<Y>()(z.second);
      }
    };

  }
}

// Forward declaration
struct Person;

// Events
//
// DEBUT       - entry into the sexually active population
// AGE         - transition between age bands
// PARTNER     - partnership formation
// BREAKUP     - partnership dissolution
// DEATH_NAT   - death unrelated to HIV infection or age-out of sexual activity
// DEATH_HIV   - HIV-related death
// INFECT      - transmission of infection
// PROGRESS    - progress to the next stage of HIV infection
// ART_INIT    - treatment initiation
// ART_CHANGE  - change in treatment status (other than INIT)
// PREP_INIT   - chemoprophylaxis initiation
// PREP_CHANGE - change in chemoprophylaxis status (other than INIT)
// RESISTANCE  - resistance emergence or reversion
// RISK_CHANGE - change sexual activity level
// MMC         - man newly circumcised
//
// We do not distinguish between ARV (ART or PrEP) change and
// termination because ARV termination can be determined using the
// individual's state, e.g., ART_NONE after an ART_CHANGE event
// implies the individual is no longer receiving ART.
enum EventType {DEBUT=0, AGE=1, PARTNER=2, BREAKUP=3, DEATH_NAT=4, DEATH_HIV=5,
		INFECT=6, PROGRESS=7,
		ART_INIT=8, ART_CHANGE=9,
		PREP_INIT=10, PREP_CHANGE=11,
		RESISTANCE=12, RISK_CHANGE=13,
		MMC=14,
		EVENTS=15};

// Person attributes
enum Sex {MALE=0, FEMALE=1, SEXES=2};
enum Activity {LEAST=0, LOW=1, MID=2, HIGH=3, LEVELS=4};

enum Stage {SUSCEPTIBLE=0,
	    ACUTE_WINDOW=1, ACUTE_DETECT=2,
	    CHRONIC_EARLY=3, CHRONIC_LATE=4,
	    AIDS=5,
	    STAGES=6};

// ART_NONE  - not on ART
// ART_EARLY - early (first-year) ART
// ART_LATE  - late (subsequent years of) ART
// ART_FAIL  - failed ART due to non-adherence
enum Treat {ART_NONE=0, ART_EARLY=1, ART_LATE=2, ART_FAIL=3, ART_STATES=4};

// PREP_NONE - not on PrEP
// PREP_HIGH - on PrEP with high drug concentrations on average
// PREP_POOR - on PrEP with poor drug concentrations on average
enum PrEP {PREP_NONE=0, PREP_HIGH=1, PREP_POOR=2, PREP_STATES=3};

// Drug sensitive virus
//  WT - wild-type
// Transmitted drug resistant variants:
//  R1, WR1 - majority, minority ART resistance (no PrEP cross-resistance)
//  C1, WC1 - majority, minority ART resistance (PrEP cross-resistance)
//  Q1, WQ1 - majority, minority PrEP resistance
// Acquired drug resistant variants
//  R2, WR2 - majority, minority ART resistance (no PrEP cross-resistance)
//  C2, WC2 - majority, minority ART resistance (PrEP cross-resistance)
//  Q2, WQ2 - majority, minority PrEP resistance
enum Virus {WT=0,
	     R1=1,  C1=2,  Q1=3,  R2=4,   C2=5,   Q2=6,
	    WR1=7, WC1=8, WQ1=9, WR2=10, WC2=11, WQ2=12,
	    VARIANTS=13};

// Age bands are not explicitly enumerated, since they depend on the
// span_age_band model parameter. For readability and for consistency
// with other state variables we define the number of age bands
// represented in the model.
enum AgeBand {AGES=8};

// Model parameter structure
struct Parameters {
public:
  Parameters();

  // Initialize epidemic parameters from a configuration file
  int read_epidemic(std::string const& filename);

  // Initialize epidemic parameters from a configuration file
  int read_intervention(std::string const& filename);

  // Set the proportion of model entrants in each activity
  // level. Proportions can be unnormalized but must be non-negative
  void prop_activity(Sex sex, std::vector<double> const& proportions);

  // Get the proportion of model entrants in each activity level
  double prop_activity(Sex sex, Activity activity) const;

  // Sample an activity level at random proportionally to the fraction
  // of the given sex who enter the model with that activity level
  Activity sample_activity_level(gsl_rng* rng, Sex sex) const;

  // Set the proportion of the initial population in each age band.
  // Proportions can be unnormalized but must be non-negative.
  void prop_age(std::vector<double> const& proportions);

  // Get the proportion of the initial population in the given age band.
  double prop_age(unsigned int band) const;

  // Sample an age and age band for an individual in the initial
  // population. This is returned as an std::pair with first
  // containing the age and second the band. These are selected so
  // that the age band probability is proportional to the weights
  // specified using prop_age. The actual age is distributed uniformly
  // at random within that band; e.g., a person in age band [15,20)
  // has P(age <= x) = (x-15)/5
  std::pair<double, unsigned int> sample_age(gsl_rng* rng) const;

  // Get the rate of transmission within a partnership.
  double rate_transmit(Person const* donor, Person const* recipient, double const time) const;

  // Get the waiting time until the next progression time for an
  // HIV-infected host
  double wait_progress(Person const* host) const;

  // initial population size
  unsigned int size_population;

  // parameters of the population growth rate
  double rate_grow_init;
  double rate_grow_late;
  double rate_grow_decr;

  // average time until non-HIV mortality by age. Inverse of
  // age-specific mortality rates
  double wait_death[SEXES][AGES];

  // age-related parameters and functions. Ages are continuous, but
  // several parameters and modules care about discrete age bands. Age
  // bands have uniform span (e.g., 5-year age groups). Ages reflect
  // the time *since sexual debut*, so a person with "age" zero is
  // assumed to have been living for 15 years
  unsigned int span_age_band;

  // proportion of population entrants who are female
  double prop_female;

  // Rate that women engaged in commercial sex work transition to
  // lower sexual activity levels
  double rate_csw_exit;

  // Relative rate of sex work ceasation by HIV stage and ART status
  double prop_risk_exit[STAGES][ART_STATES];

  // Relative propensity to initiate sex work by HIV stage and ART
  // status
  double prop_risk_init[STAGES][ART_STATES];

  // relative partner seek rate by HIV stage and ART status
  double prop_seek_decr[STAGES][ART_STATES];

  // degree of assortative mixing in the population
  double assort_act; // activity level-based assortativity
  double assort_age; // age-based assortativity (same age)
  double assort_dif; // age-based assortativity (older men with younger women)

  // base partnership seek parameter. This is the seek rate that
  // determines the average time from the end of one partnership to
  // the beginning of the next.
  double rate_partner[SEXES][LEVELS][AGES];

  // Partnership dissolution rate, stratified by partnership
  // type. This should be symmetric.
  double rate_breakup[LEVELS][LEVELS];

  // concurrency damping factors. Each person has their own damping
  // factor. damp_concurrency is used to assign the individual-level
  // parameter when a new person enters the model
  double damp_concurrency[SEXES][LEVELS];

  // Per-act transmission probability by stage of infection (assuming
  // untreated infection with wild-type virus)
  double prob_transmit[STAGES];

  // Number of acts per partnership per year. This is used to
  // calculate transmission risks. This should be symmetric.
  double num_acts[LEVELS][LEVELS];

  // Proportion of sex acts protected by condom use. This should be
  // symmetric.
  double prop_condom[LEVELS][LEVELS];

  // Relative HIV transmission risk during condom-protected sex acts
  double prop_condom_efficacy;

  // Condom-related behavioral change parameters
  double time_condom_change; // start of behavioral change (years after simulation begins)
  double span_condom_change; // time span of behavioral change (in years)
  double prop_condom_nonuse; // relative non-use after change

  // Mean duration of infection stages (waiting time in each stage)
  double wait_stage[STAGES];

  // year the HIV epidemic begins
  double year_init;

  // simulation termination time in years since the simulation began
  double year_done;

  // Average time to resistance emergence with different variants and
  // ARV states. These values are initialized to infinity by the
  // constructor, the programmer is responsible for setting the values
  // of reachable virus and ARV states where resistance emerges.
  double wait_emerge[VARIANTS][ART_STATES][PREP_STATES];

  // Average time to resistance reversion with different variants and
  // ARV states. These values are initialized to infinity by the
  // constructor, the programmer is responsible for setting the values
  // of reachable virus and ARV states where resistance reverts.
  double wait_revert[VARIANTS][ART_STATES][PREP_STATES];

  // Relative infectiousness of each viral variant
  double prop_transmit[VARIANTS];

  // Relative progression rate for each viral variant
  double prop_progress[VARIANTS];

  // Parameters specific to treatment
  struct params_art {

    // Calculates ART uptake rates in HIV disease stage h at time t
    double uptake(double const t, Stage const h) const;

    // Upper bound on the uptake rate in each stage.. Simulation
    // results may be incorrect if uptake(t,h) > rate_uptake_max at
    // any time or stage of infection. Simulation may be more
    // efficient when this bound is as tight as possible.
    double rate_uptake_max[STAGES];

    // First year of ART availability by stage of infection. For
    // simulation results to be correct, rate_uptake(t,s) must be zero
    // before ART becomes available to people in stage s. Simulation
    // may be more efficient when year_available[s] is the earliest
    // time that rate_uptake(t,s) > 0. By default,
    // year_available[s]=infinity for all stages s.
    double year_available[STAGES];

    // year scale-up ends; uptake rates stabilize after scale-up ends
    double year_scale[STAGES];

    // proportion initiating ART annually at time of peak scale-up
    double prop_scale[STAGES];

    // proportion initiating ART annually after scale-up ends
    double prop_level[STAGES];

    // rate of transition from early to late ART
    double rate_advance;

    // drop out rate
    double rate_dropout[ART_STATES];

    // non-adherence rate
    double rate_nonadhere[ART_STATES];

    // ART-related mortality rate with symptomatic infection
    double rate_mortality[ART_STATES];
    // Relative ART-related mortality by stage of infection
    double prop_mortality[STAGES];

    // Relative per-act probability of transmission according to the
    // donor's ART status
    double prop_transmit[ART_STATES];

    // Relative rate of progression by ART state
    double prop_progress[ART_STATES];

    // Proportion who acquire ART resistance who are cross-resistant
    // to PrEP
    double prop_cross_resist;

  } art;

  struct params_prep {
    // Probability of PrEP initiation at sexual debut. This must be a
    // non-negative function with value between 0 and 1.
    double (*prop_debut)(double, Sex, Activity);

    // PrEP uptake rate. This must be a non-negative bounded function
    // of time. Two members of params_prep, rate_uptake_max and
    // year_available, are used to flexibly support efficient
    // simulation of PrEP uptake
    //
    // PrEP uptake rates may be stratified by sex, sexual activity
    // level and age band
    double (*rate_uptake)(double, Sex, Activity, unsigned int);

    // maximum rate of PrEP uptake. For correct simulation results,
    // rate_uptake(t) <= rate_uptake_max at all times t. Simulation
    // may be most efficient when this bound is as tight as possible
    double rate_uptake_max;

    // time that PrEP scale-up begins. For correct simulation results,
    // rate_uptake(t) must be zero before PrEP becomes
    // available. Simulation may be most efficient when year_available
    // is the earliest time that rate_uptake(t) > 0,
    //
    //   year_available = inf{t : rate_uptake(t) > 0}
    double year_available;

    // duration of scale-up; coverage stabilizes once scale-up ends 
    double span_scale;

    // time PrEP becomes unavailable
    double year_closure;

    // maximum population-level PrEP coverage (actual coverage may be
    // lower depending on eligibility restrictions).
    double target_overall;

    // coverage targets by risk group. Eligibility is encoded as follows:
    //
    // target[g][a][k] > 0: prioritized for maximum coverage of target[g][a][k]
    // target[g][a][k] = 0: unprioritized but eligible
    // target[g][a][k] < 0: not eligible
    //
    // Prioritized individuals may initiate PrEP. Eligible individuals
    // may not initiate PrEP, but individuals in this population may
    // stay on PrEP (e.g., if women aged 15-19 are prioritized and
    // women 20-24 are eligible, then women who start at age 15-19 may
    // remain on PrEP after age 20). Ineligible individuals may not
    // initiate PrEP, and individuals on PrEP in this group stop PrEP
    // immediately at their next dosing visit
    double target[SEXES][AGES][LEVELS];

    // average time from PrEP initiation until stopping PrEP. Factors
    // in dropout and PrEP completion
    double wait_exit;

    // average time between PrEP dosing visits
    double wait_dose;

    // probability of HIV testing at a given dosing visit
    double prop_test;

    // PrEP efficacy, the relative rate of HIV acquisition when fully
    // adherent to PrEP, stratified by the donor's virus
    double efficacy[VARIANTS];

    // PrEP adherence, determines drug concentrations in each PrEP
    // state (for RPV, PrEP adherence is a misnomer; drug
    // concentrations are determined by PrEP reliability instead).
    double adherence[PREP_STATES];

    // Proportion of PrEP doses that produce high PrEP levels
    double prop_high;

  } prep;

  struct params_mmc {
    // calculate target MMC coverage and its derivative at time t
    std::pair<double, double> coverage(double const t) const;

    // Relative HIV acquisition risk among circumcised males compared to
    // uncircumcised males
    double efficacy;

    // Baseline MMC prevalence at sexual debut
    double prop_base;

    // Maximum target coverage. For correct simulation results, the
    // cover_target(t) <= target_max at all times t. Simulation may be
    // most efficient when this bound is as tight as possible.
    double target_max;

    // Maximum change in target coverage. For correct simulation
    // results, the cover_change(t) <= change_max at all times
    // t. Simulation may be most efficient when this bound is as tight
    // as possible.
    double change_max;

    // Time that MMC scale-up among sexually active men begins
    double year_init;

    // Time that MMC scale-up among sexually active men ends; coverage
    // stable after this time
    double year_done;

  } mmc;

private:
  // helper function that sets up data structures used by the alias
  // method. n is the number of elements in 'table' and 'alias'. It
  // is an error if proportions.size() != n.
  static void alias_setup(std::vector<double> const& proportions,
			  double* table, unsigned int* alias, size_t n);

  // m_activity_table and m_activity alias are the probability and
  // alias tables used to sample activity levels proportionally using
  // the alias method
  double m_activity_table[SEXES][LEVELS];
  unsigned int m_activity_alias[SEXES][LEVELS];

  // m_age_table and m_age_alias are the probability and alias tables
  // used to sample the age distribution at the start of simulation
  // proportionally using the alias method
  double m_age_table[AGES];
  unsigned int m_age_alias[AGES];

  // default ART uptake, PrEP debut and PrEP uptake rate functions
  static inline double zero_art(double, Stage) {return 0.0;}
  static inline double zero_prep(double, Sex, Activity) {return 0.0;}
  static inline double zero_prep(double, Sex, Activity, unsigned int) {return 0.0;}
};

// global model parameters
extern Parameters params;

// global alias
extern double infinity;

#endif // MODELTYPES_H
