#include <getopt.h>
#include <sys/time.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <Transmission.H>
#include "adevs.h"

using namespace std;

static struct options_t {
  options_t() : epidemic(""), seed(0) {}

  void print_help() {
    cout << "USAGE: simulate [options]\n\n"
	 << "Options:\n"
	 << "-h --help            show this help message and exit\n"
	 << "-e --epidemic=FILE   initialize naive epidemic parameters from FILE\n"
	 << "Simulation parameters:\n"
	 << "-c C --collect=C     memory management policy parameter\n"
	 << "-s SEED --seed=SEED  seed the random number generator using SEED,\n"
	 << "                     otherwise a seed is obtained from the operating system.\n"
         << "                     Ignored during deterministic simulation.\n";
  }

  void parse_args(int argc, char **argv) {
    int c(0), index(0);
    option longopts[] = {
      {"help",     no_argument,       0, 'h'},
      {"epidemic", required_argument, 0, 'e'},
      {"collect",  required_argument, 0, 'c'},
      {"seed",     required_argument, 0, 's'},
      {0, 0, 0, 0}};
    while(c >= 0) {
      c = getopt_long(argc, argv, "hc:f:s:", longopts, &index);
      switch(c) {
      case 'c': collect = strtoul(optarg, NULL, 0); break;
      case 'h': print_help(); exit(0); break;
      case 'e': epidemic = optarg; break;
      case 's': seed = strtoul(optarg, NULL, 0); break;
      default: break;
      }
    }
    if (seed == 0) {
      std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary);
      urandom.read((char*)&seed, sizeof(seed));
      urandom.close();
    }
  }

  string epidemic; // epidemic inputs

  unsigned long int seed;

  unsigned int collect;

} options;

// Top level model. Handles communication between modules
class Model : public adevs::Network<Message> {
public:
  typedef adevs::Devs<Message> Component;

  Model() : adevs::Network<Message>() {}
  virtual ~Model() {}

  // Must be called prior to simulation
  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);

    transmission.initialize(people, gsl_rng_get(rng));
    transmission.setParent(this);

    gsl_rng_free(rng);
  }

  // adevs::Network interface method
  virtual void getComponents(adevs::Set<Component*>& c) {
    c.insert(&transmission);
  }

  // adevs::Network interface method
  virtual void route(Message const& msg, Component* sender, adevs::Bag< adevs::Event<Message> >& output) {}

private:
  // Modules
  Transmission transmission;
};

unsigned long int get_seed() {
  unsigned long int seed(0);
  const int fd(open("/dev/random", O_RDONLY));
  if (fd > 0) {
    read(fd, &seed, sizeof(seed));
    close(fd);
  }
  return seed;
}

void simulate(std::vector<Person*> &P, unsigned long int seed) {
  timeval t0, t1, t2;
  gettimeofday(&t0, NULL);
  Model M;
  M.initialize(P, seed);
  adevs::Simulator<Message> simulator(&M);
  gettimeofday(&t1, NULL);

  // Distribution testing
  while (simulator.nextEventTime() < DBL_MAX) simulator.execNextEvent();

  // sp: seropositive partner, sn: seronegative partner
  cout << "sp.sex sp.level sp.stage sn.sex sn.level sn.mmc time rate\n";
  Person::History::const_iterator hi;
  for (size_t k(0); k < P.size(); ++k) {
    for (hi = P[k]->history.begin(); hi != P[k]->history.end(); ++hi) {
      if ((*hi)->type() == INFECT) {
	Person* donor(dynamic_cast<InfectRecord*>(*hi)->donor);
	cout << donor->sex() << ' '
	     << donor->activity() << ' '
	     << donor->stage() << ' '
	     << P[k]->sex() << ' '
	     << P[k]->activity() << ' '
	     << P[k]->circumcised() << ' '
	     << (*hi)->time << ' '
	     << params.rate_transmit(donor, P[k], 0.0) << '\n';
      }
    }
  }

  gettimeofday(&t2, NULL);
  printf("%% Simulation time: %0.4f %0.4f %0.4f\n",
 	 1e-6 * (t1.tv_usec - t0.tv_usec) + (t1.tv_sec - t0.tv_sec),
 	 1e-6 * (t2.tv_usec - t1.tv_usec) + (t2.tv_sec - t1.tv_sec),	 
 	 1e-6 * (t2.tv_usec - t0.tv_usec) + (t2.tv_sec - t0.tv_sec));
}

int main(int argc, char **argv) {
  options.parse_args(argc, argv);

  cout << '%';
  for (int arg(0); arg < argc; ++arg) cout << ' ' << argv[arg];
  cout << '\n';

  char hostname[256];
  gethostname(hostname, 256);
  cout << "% Host: " << hostname << '\n'
       << "% Seed: " << options.seed << '\n'
       << "% Memory mode: " << options.collect << '\n';

  if (options.epidemic.size() == 0) {
    cerr << "Error: no epidemic inputs specified, exiting." << endl;
    return -1;
  } else if (params.read_epidemic(options.epidemic) < 0) {
    cerr << "Error: failed to read epidemic inputs, exiting." << endl;
    return -1;
  } else {
    cout << "% Epidemic inputs: " << options.epidemic << '\n';
  }

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, options.seed);

  for (int ki(0); ki < LEVELS; ++ki) {
    for (int kj(0); kj < LEVELS; ++kj) {
      params.rate_breakup[ki][kj] = 1.0;
    }
  }

  const size_t N(params.size_population);

  // Initialize the model population
  std::vector<Person*> P;
  Person* donor;
  Person* recip;
  for (int h(ACUTE_WINDOW); h < STAGES; ++h) {
    for (int j(LEAST); j < LEVELS; ++j) {
      for (int k(LEAST); k < LEVELS; ++k) {

	for (size_t i(0); i < N; ++i) { // Female recipients
	  donor = new Person(0.0);
	  recip = new Person(0.0);

	  donor->sex(MALE);
	  donor->activity(Activity(j));
	  donor->age_band(0);
	  donor->circumcised(gsl_rng_uniform(rng) < 0.5);
	  donor->stage(Stage(h));
	  donor->treat(ART_NONE);
	  donor->prep(PREP_NONE);
	  donor->virus(WT);

	  recip->sex(FEMALE);
	  recip->activity(Activity(k));
	  recip->age_band(0);
	  recip->circumcised(false);
	  recip->stage(SUSCEPTIBLE);
	  recip->treat(ART_NONE);
	  recip->prep(PREP_NONE);
	  recip->virus(VARIANTS);

	  donor->partners.insert(recip);
	  recip->partners.insert(donor);

	  P.push_back(donor);
	  P.push_back(recip);
	}

	for (size_t i(0); i < N; ++i) { // Male recipients
	  // Uncircumcised
	  donor = new Person(0.0);
	  recip = new Person(0.0);

	  donor->sex(FEMALE);
	  donor->activity(Activity(j));
	  donor->age_band(0);
	  donor->circumcised(false);
	  donor->stage(Stage(h));
	  donor->treat(ART_NONE);
	  donor->prep(PREP_NONE);
	  donor->virus(WT);

	  recip->sex(MALE);
	  recip->activity(Activity(k));
	  recip->age_band(0);
	  recip->circumcised(false);
	  recip->stage(SUSCEPTIBLE);
	  recip->treat(ART_NONE);
	  recip->prep(PREP_NONE);
	  recip->virus(VARIANTS);

	  donor->partners.insert(recip);
	  recip->partners.insert(donor);

	  P.push_back(donor);
	  P.push_back(recip);

	  // Circumcised
	  donor = new Person(0.0);
	  recip = new Person(0.0);

	  donor->sex(FEMALE);
	  donor->activity(Activity(j));
	  donor->age_band(0);
	  donor->circumcised(false);
	  donor->stage(Stage(h));
	  donor->treat(ART_NONE);
	  donor->prep(PREP_NONE);
	  donor->virus(WT);

	  recip->sex(MALE);
	  recip->activity(Activity(k));
	  recip->age_band(0);
	  recip->circumcised(true);
	  recip->stage(SUSCEPTIBLE);
	  recip->treat(ART_NONE);
	  recip->prep(PREP_NONE);
	  recip->virus(VARIANTS);

	  donor->partners.insert(recip);
	  recip->partners.insert(donor);

	  P.push_back(donor);
	  P.push_back(recip);

	}

      }
    }
  }

  simulate(P, gsl_rng_get(rng));
  for (size_t n(0); n < P.size(); ++n) delete P[n];

  gsl_rng_free(rng);

  return 0;
}
