#include <getopt.h>
#include <sys/time.h>
#include <cassert>
#include <limits>
#include <numeric>
#include <cmath>
#include <fstream>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <Model.H>

using namespace std;

static struct options_t {
  options_t() : epidemic(""), arvs(""), year(2030), art(0), sec(0), prep(0) {}

  void print_help() {
    cout << "USAGE: simulate [options]\n\n"
	 << "Options:\n"
	 << "-h --help            show this help message and exit\n"
	 << "-e --epidemic=FILE   initialize epidemic parameters from FILE\n"
	 << "-a --arvs=FILE       initialize intervention parameters from FILE\n"
	 << "Interventions:\n"
	 << "--art                enable treatment with ART\n"
	 << "--sec                enable second-line ART\n"
	 << "--prep               enable PrEP (includes window use)\n"
	 << "Simulation parameters:\n"
	 << "-y Y --year=Y        simulate until year Y (default: 2030)\n";
  }

  void parse_args(int argc, char **argv) {
    int c(0), index(0);
    option longopts[] = {
      {"help",     no_argument,       0, 'h'},
      {"epidemic", required_argument, 0, 'e'},
      {"arvs",     required_argument, 0, 'a'},
      {"art",      no_argument,       &art, 1},
      {"sec",      no_argument,       &sec, 1},
      {"prep",     no_argument,       &prep, 1},
      {"year",     required_argument, 0, 'y'},
      {0, 0, 0, 0}};
    while(c >= 0) {
      c = getopt_long(argc, argv, "ha:e:y:", longopts, &index);
      switch(c) {
      case 'h': print_help(); exit(0); break;
      case 'a': arvs = optarg; break;
      case 'e': epidemic = optarg; break;
      case 'y': year = strtoul(optarg, NULL, 0); break;
      default: break;
      }
    }
  }

  string epidemic; // filename for epidemic inputs
  string arvs;     // filename for intervention inputs

  unsigned int year; // year simulation ends

  int art;  // 1 if ART available, 0 otherwise
  int sec;  // 1 if second-line ART available, 0 otherwise
  int prep; // 1 if PrEP available, 0 otherwise

} options;

void header(ostream &out);
void report(ostream &out, double const t, double const y[], Parameters const& p);

int main(int argc, char** argv) {
  options.parse_args(argc, argv);

  cout << '%';
  for (int arg(0); arg < argc; ++arg) cout << ' ' << argv[arg];
  cout << '\n';

  char hostname[256];
  gethostname(hostname, 256);
  cout << "% Host: " << hostname << '\n';

  Parameters parameters;
  int status;

  status = parameters.read_epidemic(options.epidemic);
  if (status < 0) {
    cerr << "Error: could not read epidemic inputs (code=" << status << ")\n";
    return -1;
  } else {
    cout << "% Epidemic inputs: " << options.epidemic << '\n';
  }

  status = parameters.read_intervention(options.arvs);
  if (status < 0) {
    cout << "% Intervention inputs: " << options.arvs << '\n';
    cerr << "Error: could not read intervention inputs (code=" << status << ")\n";
    return -1;
  } else {
    cout << "% Intervention inputs: " << options.arvs << '\n';
  }

  const unsigned int init_year(1978);
  const unsigned int n(COMPARTMENTS);
  const unsigned int years(options.year - init_year);
  double t(init_year), y[n];

  std::fill(y, y + n, 0.0);

  // initialize parameters for outcome discounting
#warning "Cost-effectiveness parameters hard-coded"
  parameters.ce_discount = 0.03;
  parameters.ce_tinit = parameters.prep_time_rollout;

  if (options.art) {
    cout << "% With ART\n";
  } else {
    // disable ART by delaying rollout indefinitely
    std::fill(parameters.art_time_start, parameters.art_time_start + STAGES, std::numeric_limits<double>::infinity());
    cout << "% Without ART\n";
  }

  if (options.sec) {
    cout << "% With 2nd-line ART\n";
  } else {
    // disable second-line ART by setting the switch rate to zero
    parameters.sec_init_dr = 0.0;
    parameters.sec_init_na = 0.0;
    cout << "% Without 2nd-line ART\n";
  }

  if (options.prep) {
    cout << "% With PrEP\n";
  } else {
    // disable PrEP by delaying rollout indefinitely
    parameters.prep_time_rollout = std::numeric_limits<double>::infinity();
    cout << "% Without PrEP\n";
  }

  header(cout);
  initialize(y, parameters);

  // Age band "width" ignores the 55+ age band
  const int span((AGES - 1) / (BANDS - 1));
  double num;
  double *z;

  // maps sex and circumcision status to sex
  // F->F, MU->M, MC->M
  const Sex sex[] = {F, M, M};

  std::fill(y, y + COMPARTMENTS, 0.0);
  for (int g(0); g < SEXCIRC; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS - 1; ++b) {
  	num = parameters.init_size * parameters.init_risk[sex[g]][b][k] / static_cast<double>(span);
      	for (int a(0); a < span; ++a) {
  	  z = offset(g, k, b * span + a, ANTE) + y;
  	  if (g == F) {
  	    z[Y1NWT] = num;
  	  } else if (g == MU) {
  	    z[Y1NWT] = num * (1.0 - parameters.mmc_prop_init);
	  } else {
  	    z[Y1NWT] = num * parameters.mmc_prop_init;
  	  }
  	}
      }
    }
  }

  report(cout, t, y, parameters);

  timeval st0, st1;
  gettimeofday(&st0, NULL);

  timeval ut0, ut1;

  const double hstart(0.05); // suggested initial step size
  gsl_odeiv2_system sys = {model, NULL, n, &parameters};
  gsl_odeiv2_driver* drv(gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, hstart, 1e-4, 1e-4));
  for (unsigned int step(0); step < years; ++step) {
    gettimeofday(&ut0, NULL);
    status = gsl_odeiv2_driver_apply(drv, &t, t + 1, y);
    gettimeofday(&ut1, NULL);
    fprintf(stderr, "%0.0f %0.4f\n", t, 1e-6 * (ut1.tv_usec - ut0.tv_usec) + (ut1.tv_sec - ut0.tv_sec));

    if (status != GSL_SUCCESS) {
      cerr << "Warning: GSL is complaining at step " << step << ", s=" << status << '\n';
    }

    report(cout, t, y, parameters);
  }

  gsl_odeiv2_driver_free(drv);

  gettimeofday(&st1, NULL);

  printf("%% Simulation time: %0.4f\n", 1e-6 * (st1.tv_usec - st0.tv_usec) + (st1.tv_sec - st0.tv_sec));

  return 0;
}

void header(ostream &out) {
  const char sex[] = {'w','m'};

  // outcomes among 15-54
  out << "year X P1 P2 L1 L2 L3 A1 inci";

  // infected+susceptible by sex and age
  for (int g(0); g < SEXES; ++g) {
    for (int b(0); b < BANDS; ++b) out << " n" << sex[g] << 'b' << b;
  }

  // infected by sex and age
  for (int g(0); g < SEXES; ++g) {
    for (int b(0); b < BANDS; ++b) out << " y" << sex[g] << 'b' << b;
  }

  // incident cases by sex and age
  for (int g(0); g < SEXES; ++g) {
    for (int b(0); b < BANDS; ++b) out << " i" << sex[g] << 'b' << b;
  }

  // infected+susceptible by sex and activity level (age 15-54)
  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) out << " n" << sex[g] << 'k' << k;
  }

  // infected by sex and activity level (age 15-54)
  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) out << " y" << sex[g] << 'k' << k;
  }

  // incident cases by sex and activity level (age 15-54)
  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) out << " i" << sex[g] << 'k' << k;
  }

  // These outcomes are reported among persons age 15-54
  out << " prp.X"
      << " prp.P1 prp.P2 prp.L1 prp.L2 prp.L3 prp.A1"
      << " art.P1 art.P2 art.L1 art.L2 art.L3 art.A1" // both lines
      << " sec.P1 sec.P2 sec.L1 sec.L2 sec.L3 sec.A1" // second line
      << " mmc.X mmc.Y"
      << " DR R1 R2 C1 C2 Q1 Q2 S1 S2 WT r1 r2 c1 c2 q1 q2 s1 s2";

  // cumulative outcomes since simulation start date
  // men and women
  out << " active.B.lyl.X"
      << " active.B.lyl.P1 active.B.lyl.P2 active.B.lyl.L1 active.B.lyl.L2 active.B.lyl.L3 active.B.lyl.A1" // life-years lived
      << " active.B.inf.WT active.B.inf.R1 active.B.inf.C1 active.B.inf.Q1 active.B.inf.S1" // infections by transmitted variant
      << " active.B.mrt.D active.B.mrt.B" // disease (D) and background (B) mortality
      << " active.B.1st.P1 active.B.1st.P2 active.B.1st.L1 active.B.1st.L2 active.B.1st.L3 active.B.1st.A1" // person-years of 1st-line ART
      << " active.B.2nd.P1 active.B.2nd.P2 active.B.2nd.L1 active.B.2nd.L2 active.B.2nd.L3 active.B.2nd.A1" // person-years of 2nd-line ART
      << " active.B.ainit" // ART initiation
      << " active.B.prp"; // person-years of PrEP

  // women
  out << " active.W.lyl.X"
      << " active.W.lyl.P1 active.W.lyl.P2 active.W.lyl.L1 active.W.lyl.L2 active.W.lyl.L3 active.W.lyl.A1" // life-years lived
      << " active.W.inf.WT active.W.inf.R1 active.W.inf.C1 active.W.inf.Q1 active.W.inf.S1" // infections by transmitted variant
      << " active.W.mrt.D active.W.mrt.B" // disease (D) and background (B) mortality
      << " active.W.1st.P1 active.W.1st.P2 active.W.1st.L1 active.W.1st.L2 active.W.1st.L3 active.W.1st.A1" // person-years of 1st-line ART
      << " active.W.2nd.P1 active.W.2nd.P2 active.W.2nd.L1 active.W.2nd.L2 active.W.2nd.L3 active.W.2nd.A1" // person-years of 2nd-line ART
      << " active.W.ainit" // ART initiation
      << " active.W.prp"; // person-years of PrEP

  // men
  out << " active.M.lyl.X"
      << " active.M.lyl.P1 active.M.lyl.P2 active.M.lyl.L1 active.M.lyl.L2 active.M.lyl.L3 active.M.lyl.A1" // life-years lived
      << " active.M.inf.WT active.M.inf.R1 active.M.inf.C1 active.M.inf.Q1 active.M.inf.S1" // infections by transmitted variant
      << " active.M.mrt.D active.M.mrt.B" // disease (D) and background (B) mortality
      << " active.M.1st.P1 active.M.1st.P2 active.M.1st.L1 active.M.1st.L2 active.M.1st.L3 active.M.1st.A1" // person-years of 1st-line ART
      << " active.M.2nd.P1 active.M.2nd.P2 active.M.2nd.L1 active.M.2nd.L2 active.M.2nd.L3 active.M.2nd.A1" // person-years of 2nd-line ART
      << " active.M.ainit" // ART initiation
      << " active.M.prp" // person-years of PrEP
      << " active.M.mmc.D active.M.mmc.U"; // circumcisions at debut (D) or uptake after debut (U)

  // cumulative outcomes prior to PrEP closure
  // men and women
  out << " cohort.B.lyl.X"
      << " cohort.B.lyl.P1 cohort.B.lyl.P2 cohort.B.lyl.L1 cohort.B.lyl.L2 cohort.B.lyl.L3 cohort.B.lyl.A1" // life-years lived
      << " cohort.B.inf.WT cohort.B.inf.R1 cohort.B.inf.C1 cohort.B.inf.Q1 cohort.B.inf.S1" // infections by transmitted variant
      << " cohort.B.mrt.D cohort.B.mrt.B" // disease (D) and background (B) mortality
      << " cohort.B.1st.P1 cohort.B.1st.P2 cohort.B.1st.L1 cohort.B.1st.L2 cohort.B.1st.L3 cohort.B.1st.A1" // person-years of 1st-line ART
      << " cohort.B.2nd.P1 cohort.B.2nd.P2 cohort.B.2nd.L1 cohort.B.2nd.L2 cohort.B.2nd.L3 cohort.B.2nd.A1" // person-years of 2nd-line ART
      << " cohort.B.ainit" // ART initiation
      << " cohort.B.prp"; // person-years of PrEP

  // cumulative outcomes prior to PrEP closure
  // women
  out << " cohort.W.lyl.X"
      << " cohort.W.lyl.P1 cohort.W.lyl.P2 cohort.W.lyl.L1 cohort.W.lyl.L2 cohort.W.lyl.L3 cohort.W.lyl.A1" // life-years lived
      << " cohort.W.inf.WT cohort.W.inf.R1 cohort.W.inf.C1 cohort.W.inf.Q1 cohort.W.inf.S1" // infections by transmitted variant
      << " cohort.W.mrt.D cohort.W.mrt.B" // disease (D) and background (B) mortality
      << " cohort.W.1st.P1 cohort.W.1st.P2 cohort.W.1st.L1 cohort.W.1st.L2 cohort.W.1st.L3 cohort.W.1st.A1" // person-years of 1st-line ART
      << " cohort.W.2nd.P1 cohort.W.2nd.P2 cohort.W.2nd.L1 cohort.W.2nd.L2 cohort.W.2nd.L3 cohort.W.2nd.A1" // person-years of 2nd-line ART
      << " cohort.W.ainit" // ART initiation
      << " cohort.W.prp"; // person-years of PrEP

  // cumulative outcomes prior to PrEP closure
  // men
  out << " cohort.M.lyl.X"
      << " cohort.M.lyl.P1 cohort.M.lyl.P2 cohort.M.lyl.L1 cohort.M.lyl.L2 cohort.M.lyl.L3 cohort.M.lyl.A1" // life-years lived
      << " cohort.M.inf.WT cohort.M.inf.R1 cohort.M.inf.C1 cohort.M.inf.Q1 cohort.M.inf.S1" // infections by transmitted variant
      << " cohort.M.mrt.D cohort.M.mrt.B" // disease (D) and background (B) mortality
      << " cohort.M.1st.P1 cohort.M.1st.P2 cohort.M.1st.L1 cohort.M.1st.L2 cohort.M.1st.L3 cohort.M.1st.A1" // person-years of 1st-line ART
      << " cohort.M.2nd.P1 cohort.M.2nd.P2 cohort.M.2nd.L1 cohort.M.2nd.L2 cohort.M.2nd.L3 cohort.M.2nd.A1" // person-years of 2nd-line ART
      << " cohort.M.ainit" // ART initiation
      << " cohort.M.prp" // person-years of PrEP
      << " cohort.M.mmc.D cohort.M.mmc.U"; // circumcisions at debut (D) or uptake after debut (U)

  // discounted cumulative outcomes
  out << " actdsc.B.lyl.X"
      << " actdsc.B.lyl.P1 actdsc.B.lyl.P2 actdsc.B.lyl.L1 actdsc.B.lyl.L2 actdsc.B.lyl.L3 actdsc.B.lyl.A1" // life-years lived
      << " actdsc.B.inf.WT actdsc.B.inf.R1 actdsc.B.inf.C1 actdsc.B.inf.Q1 actdsc.B.inf.S1" // infections by transmitted variant
      << " actdsc.B.mrt.D actdsc.B.mrt.B" // disease (D) and background (B) mortality
      << " actdsc.B.1st.P1 actdsc.B.1st.P2 actdsc.B.1st.L1 actdsc.B.1st.L2 actdsc.B.1st.L3 actdsc.B.1st.A1" // person-years of 1st-line ART
      << " actdsc.B.2nd.P1 actdsc.B.2nd.P2 actdsc.B.2nd.L1 actdsc.B.2nd.L2 actdsc.B.2nd.L3 actdsc.B.2nd.A1" // person-years of 2nd-line ART
      << " actdsc.B.ainit" // ART initiation
      << " actdsc.B.prp"; // person-years of PrEP

  out << " actdsc.W.lyl.X"
      << " actdsc.W.lyl.P1 actdsc.W.lyl.P2 actdsc.W.lyl.L1 actdsc.W.lyl.L2 actdsc.W.lyl.L3 actdsc.W.lyl.A1" // life-years lived
      << " actdsc.W.inf.WT actdsc.W.inf.R1 actdsc.W.inf.C1 actdsc.W.inf.Q1 actdsc.W.inf.S1" // infections by transmitted variant
      << " actdsc.W.mrt.D actdsc.W.mrt.B" // disease (D) and background (B) mortality
      << " actdsc.W.1st.P1 actdsc.W.1st.P2 actdsc.W.1st.L1 actdsc.W.1st.L2 actdsc.W.1st.L3 actdsc.W.1st.A1" // person-years of 1st-line ART
      << " actdsc.W.2nd.P1 actdsc.W.2nd.P2 actdsc.W.2nd.L1 actdsc.W.2nd.L2 actdsc.W.2nd.L3 actdsc.W.2nd.A1" // person-years of 2nd-line ART
      << " actdsc.W.ainit" // ART initiation
      << " actdsc.W.prp"; // person-years of PrEP

  out << " actdsc.M.lyl.X"
      << " actdsc.M.lyl.P1 actdsc.M.lyl.P2 actdsc.M.lyl.L1 actdsc.M.lyl.L2 actdsc.M.lyl.L3 actdsc.M.lyl.A1" // life-years lived
      << " actdsc.M.inf.WT actdsc.M.inf.R1 actdsc.M.inf.C1 actdsc.M.inf.Q1 actdsc.M.inf.S1" // infections by transmitted variant
      << " actdsc.M.mrt.D actdsc.M.mrt.B" // disease (D) and background (B) mortality
      << " actdsc.M.1st.P1 actdsc.M.1st.P2 actdsc.M.1st.L1 actdsc.M.1st.L2 actdsc.M.1st.L3 actdsc.M.1st.A1" // person-years of 1st-line ART
      << " actdsc.M.2nd.P1 actdsc.M.2nd.P2 actdsc.M.2nd.L1 actdsc.M.2nd.L2 actdsc.M.2nd.L3 actdsc.M.2nd.A1" // person-years of 2nd-line ART
      << " actdsc.M.ainit" // ART initiation
      << " actdsc.M.prp" // person-years of PrEP
      << " actdsc.M.mmc.D actdsc.M.mmc.U"; // circumcisions at debut (D) or uptake after debut (U)

  out << " cohdsc.B.lyl.X"
      << " cohdsc.B.lyl.P1 cohdsc.B.lyl.P2 cohdsc.B.lyl.L1 cohdsc.B.lyl.L2 cohdsc.B.lyl.L3 cohdsc.B.lyl.A1" // life-years lived
      << " cohdsc.B.inf.WT cohdsc.B.inf.R1 cohdsc.B.inf.C1 cohdsc.B.inf.Q1 cohdsc.B.inf.S1" // infections by transmitted variant
      << " cohdsc.B.mrt.D cohdsc.B.mrt.B" // disease (D) and background (B) mortality
      << " cohdsc.B.1st.P1 cohdsc.B.1st.P2 cohdsc.B.1st.L1 cohdsc.B.1st.L2 cohdsc.B.1st.L3 cohdsc.B.1st.A1" // person-years of 1st-line ART
      << " cohdsc.B.2nd.P1 cohdsc.B.2nd.P2 cohdsc.B.2nd.L1 cohdsc.B.2nd.L2 cohdsc.B.2nd.L3 cohdsc.B.2nd.A1" // person-years of 2nd-line ART
      << " cohdsc.B.ainit" // ART initiation
      << " cohdsc.B.prp"; // person-years of PrEP

  out << " cohdsc.W.lyl.X"
      << " cohdsc.W.lyl.P1 cohdsc.W.lyl.P2 cohdsc.W.lyl.L1 cohdsc.W.lyl.L2 cohdsc.W.lyl.L3 cohdsc.W.lyl.A1" // life-years lived
      << " cohdsc.W.inf.WT cohdsc.W.inf.R1 cohdsc.W.inf.C1 cohdsc.W.inf.Q1 cohdsc.W.inf.S1" // infections by transmitted variant
      << " cohdsc.W.mrt.D cohdsc.W.mrt.B" // disease (D) and background (B) mortality
      << " cohdsc.W.1st.P1 cohdsc.W.1st.P2 cohdsc.W.1st.L1 cohdsc.W.1st.L2 cohdsc.W.1st.L3 cohdsc.W.1st.A1" // person-years of 1st-line ART
      << " cohdsc.W.2nd.P1 cohdsc.W.2nd.P2 cohdsc.W.2nd.L1 cohdsc.W.2nd.L2 cohdsc.W.2nd.L3 cohdsc.W.2nd.A1" // person-years of 2nd-line ART
      << " cohdsc.W.ainit" // ART initiation
      << " cohdsc.W.prp"; // person-years of PrEP

  out << " cohdsc.M.lyl.X"
      << " cohdsc.M.lyl.P1 cohdsc.M.lyl.P2 cohdsc.M.lyl.L1 cohdsc.M.lyl.L2 cohdsc.M.lyl.L3 cohdsc.M.lyl.A1" // life-years lived
      << " cohdsc.M.inf.WT cohdsc.M.inf.R1 cohdsc.M.inf.C1 cohdsc.M.inf.Q1 cohdsc.M.inf.S1" // infections by transmitted variant
      << " cohdsc.M.mrt.D cohdsc.M.mrt.B" // disease (D) and background (B) mortality
      << " cohdsc.M.1st.P1 cohdsc.M.1st.P2 cohdsc.M.1st.L1 cohdsc.M.1st.L2 cohdsc.M.1st.L3 cohdsc.M.1st.A1" // person-years of 1st-line ART
      << " cohdsc.M.2nd.P1 cohdsc.M.2nd.P2 cohdsc.M.2nd.L1 cohdsc.M.2nd.L2 cohdsc.M.2nd.L3 cohdsc.M.2nd.A1" // person-years of 2nd-line ART
      << " cohdsc.M.ainit" // ART initiation
      << " cohdsc.M.prp" // person-years of PrEP
      << " cohdsc.M.mmc.D cohdsc.M.mmc.U"; // circumcisions at debut (D) or uptake after debut (U)

  // infected+susceptible by sex, activity level and age  
  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS; ++b) out << " n" << sex[g] << 'k' << k << 'b' << b;
    }
  }

  // infected by sex, activity level and age  
  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS; ++b) out << " y" << sex[g] << 'k' << k << 'b' << b;
    }
  }

  // susceptible on PrEP by sex, activity level and age  
  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS; ++b) out << " p" << sex[g] << 'k' << k << 'b' << b;
    }
  }

  out << "\n";
}

void report(ostream &out, double const t, double const y[], Parameters const& p) {
  int b;

  double const* z;
  double f[SEXCIRC][BANDS][LEVELS][PREP_STATES][TRANSMIT];
  double num[SEXES][BANDS][LEVELS]; // infected+susceptible
  double inf[SEXES][BANDS][LEVELS]; // infected
  double inc[SEXES][BANDS][LEVELS]; // incident cases
  double prp[SEXES][BANDS][LEVELS]; // susceptible on PrEP
  double inf_stage[STAGES+1]; // age 15-54 infected by stage
  double prp_stage[STAGES+1]; // age 15-54 on PrEP by stage
  double art_stage[STAGES]; // age 15-54 on 1st-line ART by stage
  double sec_stage[STAGES]; // age 15-54 on 2nd-line ART by stage
  double mmc_sus(0.0), mmc_inf(0.0); // age 15-54 circumcised men
  double virus[VARIANTS]; // age 15-54 resistance prevalence
  double inci(0.0); // age 15-54 incident cases
  double sum;

  // maps sex and circumcision status to sex
  // F->F, MU->M, MC->M
  const Sex sex[] = {F, M, M};

  std::fill(inf_stage, inf_stage + STAGES + 1, 0.0);
  std::fill(prp_stage, prp_stage + STAGES + 1, 0.0);
  std::fill(art_stage, art_stage + STAGES, 0.0);
  std::fill(sec_stage, sec_stage + STAGES, 0.0);
  std::fill(virus, virus + VARIANTS, 0.0);

  for (b = 0; b < BANDS; ++b) {
    std::fill(num[F][b], num[F][b] + LEVELS, 0.0);
    std::fill(num[M][b], num[M][b] + LEVELS, 0.0);
    std::fill(inf[F][b], inf[F][b] + LEVELS, 0.0);
    std::fill(inf[M][b], inf[M][b] + LEVELS, 0.0);
    std::fill(inc[F][b], inc[F][b] + LEVELS, 0.0);
    std::fill(inc[M][b], inc[M][b] + LEVELS, 0.0);
    std::fill(prp[F][b], prp[F][b] + LEVELS, 0.0);
    std::fill(prp[M][b], prp[M][b] + LEVELS, 0.0);
  }

  force(t, y, p, f);

  // main accounting loop
  for (int u(0); u < TRACKS; ++u) {
    for (int a(0); a < AGES; ++a) {
      b = band[a];
      for (int g(0); g < SEXCIRC; ++g) {
	for (int k(0); k < LEVELS; ++k) {
	  z = offset(g, k, a, u) + y;

	  num[sex[g]][b][k] += z[XH] + z[XP] + z[XN];
	  prp[sex[g]][b][k] += z[XH] + z[XP];

	  sum = 0.0;
	  for (int v(0); v < TRANSMIT; ++v) {
	    sum += f[g][b][k][NONE][v] * z[XN] + f[g][b][k][PREP_HIGH][v] * z[XH] + f[g][b][k][PREP_POOR][v] * z[XP];
	  }
	  inc[sex[g]][b][k] += sum;

	  if (a < AGES - 1) {
	    inci += sum;
	    inf_stage[0] += z[XH] + z[XP] + z[XN];
	    prp_stage[0] += z[XH] + z[XP];
	    mmc_sus += (g == MC) * (z[XH] + z[XP] + z[XN]);
	  }

	  for (int s(Y1NWT); s < STATES; ++s) {
	    num[sex[g]][b][k] += z[s];
	    inf[sex[g]][b][k] += z[s];
	    if (a < AGES - 1) {
	      inf_stage[attr_stage[s]+1] += z[s];
	      virus[attr_virus[s]] += z[s];
	      mmc_inf += (g == MC) * z[s];
	      if (attr_drug[s] == PREP_HIGH || attr_drug[s] == PREP_POOR) {
		prp_stage[attr_stage[s]+1] += z[s];
	      } else if (attr_drug[s] == ART1_NEW || attr_drug[s] == ART1_OLD || attr_drug[s] == ART1_NA) {
		art_stage[attr_stage[s]] += z[s];
	      } else if (attr_drug[s] == ART2 || attr_drug[s] == ART2_NA) {
		art_stage[attr_stage[s]] += z[s];
		sec_stage[attr_stage[s]] += z[s];
	      }
	    }
	  }

	}
      }
    }
  }

  // print outcomes
  out << std::fixed;
  out.precision(4);
  out << t;
  for (int h(0); h < STAGES + 1; ++h) out << ' ' << inf_stage[h];
  out << ' ' << inci;

  for (int g(0); g < SEXES; ++g) {
    for (int b(0); b < BANDS; ++b) {
      sum = std::accumulate(num[g][b], num[g][b] + LEVELS, 0.0);
      out << ' ' << sum;
    }
  }

  for (int g(0); g < SEXES; ++g) {
    for (int b(0); b < BANDS; ++b) {
      sum = std::accumulate(inf[g][b], inf[g][b] + LEVELS, 0.0);
      out << ' ' << sum;
    }
  }

  for (int g(0); g < SEXES; ++g) {
    for (int b(0); b < BANDS; ++b) {
      sum = std::accumulate(inc[g][b], inc[g][b] + LEVELS, 0.0);
      out << ' ' << sum;
    }
  }

  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      sum = 0.0;
      for (int b(0); b < BANDS - 1; ++b) sum += num[g][b][k];
      out << ' ' << sum;
    }
  }

  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      sum = 0.0;
      for (int b(0); b < BANDS - 1; ++b) sum += inf[g][b][k];
      out << ' ' << sum;
    }
  }

  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      sum = 0.0;
      for (int b(0); b < BANDS - 1; ++b) sum += inc[g][b][k];
      out << ' ' << sum;
    }
  }

  for (int h(0); h < STAGES + 1; ++h) out << ' ' << prp_stage[h];
  for (int h(0); h < STAGES; ++h) out << ' ' << art_stage[h];
  for (int h(0); h < STAGES; ++h) out << ' ' << sec_stage[h];
  out << ' ' << mmc_sus << ' ' << mmc_inf
      << ' ' << virus[R1] + virus[R2] + virus[C1] + virus[C2] + virus[Q1] + virus[Q2] + virus[S1] + virus[S2]
      << ' ' << virus[R1] << ' ' << virus[R2]
      << ' ' << virus[C1] << ' ' << virus[C2]
      << ' ' << virus[Q1] << ' ' << virus[Q2]
      << ' ' << virus[S1] << ' ' << virus[S2]
      << ' ' << virus[WT]
      << ' ' << virus[WR1] << ' ' << virus[WR2]
      << ' ' << virus[WC1] << ' ' << virus[WC2]
      << ' ' << virus[WQ1] << ' ' << virus[WQ2]
      << ' ' << virus[WS1] << ' ' << virus[WS2];

  // undiscounted cumulative outcomes: since simulation start date
  out << ' ' << y[ACTIVE_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[active_lyl[h]];
  out << ' ' << y[ACTIVE_WT] << ' ' << y[ACTIVE_R1] << ' ' << y[ACTIVE_C1] << ' ' << y[ACTIVE_Q1] << ' ' << y[ACTIVE_S1]
      << ' ' << y[ACTIVE_HIV_MORT] << ' ' << y[ACTIVE_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[active_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[active_2nd[h]];
  out << ' ' << y[ACTIVE_AINIT]
      << ' ' << y[ACTIVE_PRP];

  out << ' ' << y[ACTIVE_W_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[active_w_lyl[h]];
  out << ' ' << y[ACTIVE_W_WT] << ' ' << y[ACTIVE_W_R1] << ' ' << y[ACTIVE_W_C1] << ' ' << y[ACTIVE_W_Q1] << ' ' << y[ACTIVE_W_S1]
      << ' ' << y[ACTIVE_W_HIV_MORT] << ' ' << y[ACTIVE_W_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[active_w_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[active_w_2nd[h]];
  out << ' ' << y[ACTIVE_W_AINIT]
      << ' ' << y[ACTIVE_W_PRP];

  out << ' ' << y[ACTIVE_M_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[active_m_lyl[h]];
  out << ' ' << y[ACTIVE_M_WT] << ' ' << y[ACTIVE_M_R1] << ' ' << y[ACTIVE_M_C1] << ' ' << y[ACTIVE_M_Q1] << ' ' << y[ACTIVE_M_S1]
      << ' ' << y[ACTIVE_M_HIV_MORT] << ' ' << y[ACTIVE_M_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[active_m_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[active_m_2nd[h]];
  out << ' ' << y[ACTIVE_M_AINIT]
      << ' ' << y[ACTIVE_M_PRP]
      << ' ' << y[ACTIVE_MMC_DEBUT] << ' ' << y[ACTIVE_MMC_UPTAKE];

  // undiscounted cumulative outcomes: in cohort alive prior to PrEP closure
  out << ' ' << y[COHORT_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohort_lyl[h]];
  out << ' ' << y[COHORT_WT] << ' ' << y[COHORT_R1] << ' ' << y[COHORT_C1] << ' ' << y[COHORT_Q1] << ' ' << y[COHORT_S1]
      << ' ' << y[COHORT_HIV_MORT] << ' ' << y[COHORT_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohort_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohort_2nd[h]];
  out << ' ' << y[COHORT_AINIT]
      << ' ' << y[COHORT_PRP];

  out << ' ' << y[COHORT_W_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohort_w_lyl[h]];
  out << ' ' << y[COHORT_W_WT] << ' ' << y[COHORT_W_R1] << ' ' << y[COHORT_W_C1] << ' ' << y[COHORT_W_Q1] << ' ' << y[COHORT_W_S1]
      << ' ' << y[COHORT_W_HIV_MORT] << ' ' << y[COHORT_W_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohort_w_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohort_w_2nd[h]];
  out << ' ' << y[COHORT_W_AINIT]
      << ' ' << y[COHORT_W_PRP];

  out << ' ' << y[COHORT_M_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohort_m_lyl[h]];
  out << ' ' << y[COHORT_M_WT] << ' ' << y[COHORT_M_R1] << ' ' << y[COHORT_M_C1] << ' ' << y[COHORT_M_Q1] << ' ' << y[COHORT_M_S1]
      << ' ' << y[COHORT_M_HIV_MORT] << ' ' << y[COHORT_M_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohort_m_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohort_m_2nd[h]];
  out << ' ' << y[COHORT_M_AINIT]
      << ' ' << y[COHORT_M_PRP]
      << ' ' << y[COHORT_MMC_DEBUT] << ' ' << y[COHORT_MMC_UPTAKE];

  // discounted cumulative outcomes: since simulation start date
  out << ' ' << y[ACTDSC_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[actdsc_lyl[h]];
  out << ' ' << y[ACTDSC_WT] << ' ' << y[ACTDSC_R1] << ' ' << y[ACTDSC_C1] << ' ' << y[ACTDSC_Q1] << ' ' << y[ACTDSC_S1]
      << ' ' << y[ACTDSC_HIV_MORT] << ' ' << y[ACTDSC_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[actdsc_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[actdsc_2nd[h]];
  out << ' ' << y[ACTDSC_AINIT]
      << ' ' << y[ACTDSC_PRP];

  out << ' ' << y[ACTDSC_W_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[actdsc_w_lyl[h]];
  out << ' ' << y[ACTDSC_W_WT] << ' ' << y[ACTDSC_W_R1] << ' ' << y[ACTDSC_W_C1] << ' ' << y[ACTDSC_W_Q1] << ' ' << y[ACTDSC_W_S1]
      << ' ' << y[ACTDSC_W_HIV_MORT] << ' ' << y[ACTDSC_W_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[actdsc_w_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[actdsc_w_2nd[h]];
  out << ' ' << y[ACTDSC_W_AINIT]
      << ' ' << y[ACTDSC_W_PRP];

  out << ' ' << y[ACTDSC_M_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[actdsc_m_lyl[h]];
  out << ' ' << y[ACTDSC_M_WT] << ' ' << y[ACTDSC_M_R1] << ' ' << y[ACTDSC_M_C1] << ' ' << y[ACTDSC_M_Q1] << ' ' << y[ACTDSC_M_S1]
      << ' ' << y[ACTDSC_M_HIV_MORT] << ' ' << y[ACTDSC_M_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[actdsc_m_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[actdsc_m_2nd[h]];
  out << ' ' << y[ACTDSC_M_AINIT]
      << ' ' << y[ACTDSC_M_PRP]
      << ' ' << y[ACTDSC_MMC_DEBUT] << ' ' << y[ACTDSC_MMC_UPTAKE];

  // discounted cumulative outcomes: in cohort alive prior to PrEP closure
  out << ' ' << y[COHDSC_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohdsc_lyl[h]];
  out << ' ' << y[COHDSC_WT] << ' ' << y[COHDSC_R1] << ' ' << y[COHDSC_C1] << ' ' << y[COHDSC_Q1] << ' ' << y[COHDSC_S1]
      << ' ' << y[COHDSC_HIV_MORT] << ' ' << y[COHDSC_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohdsc_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohdsc_2nd[h]];
  out << ' ' << y[COHDSC_AINIT]
      << ' ' << y[COHDSC_PRP];

  out << ' ' << y[COHDSC_W_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohdsc_w_lyl[h]];
  out << ' ' << y[COHDSC_W_WT] << ' ' << y[COHDSC_W_R1] << ' ' << y[COHDSC_W_C1] << ' ' << y[COHDSC_W_Q1] << ' ' << y[COHDSC_W_S1]
      << ' ' << y[COHDSC_W_HIV_MORT] << ' ' << y[COHDSC_W_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohdsc_w_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohdsc_w_2nd[h]];
  out << ' ' << y[COHDSC_W_AINIT]
      << ' ' << y[COHDSC_W_PRP];

  out << ' ' << y[COHDSC_M_LYL_X];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohdsc_m_lyl[h]];
  out << ' ' << y[COHDSC_M_WT] << ' ' << y[COHDSC_M_R1] << ' ' << y[COHDSC_M_C1] << ' ' << y[COHDSC_M_Q1] << ' ' << y[COHDSC_M_S1]
      << ' ' << y[COHDSC_M_HIV_MORT] << ' ' << y[COHDSC_M_BAK_MORT];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohdsc_m_art[h]];
  for (int h(0); h < STAGES; ++h) out << ' ' << y[cohdsc_m_2nd[h]];
  out << ' ' << y[COHDSC_M_AINIT]
      << ' ' << y[COHDSC_M_PRP]
      << ' ' << y[COHDSC_MMC_DEBUT] << ' ' << y[COHDSC_MMC_UPTAKE];

  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS; ++b) out << ' ' << num[g][b][k];
    }
  }

  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS; ++b) out << ' ' << inf[g][b][k];
    }
  }

  for (int g(0); g < SEXES; ++g) {
    for (int k(0); k < LEVELS; ++k) {
      for (int b(0); b < BANDS; ++b) out << ' ' << prp[g][b][k];
    }
  }

  out << endl;
}
