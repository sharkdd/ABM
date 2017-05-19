#include <getopt.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <Bookkeeper.H>
#include <Demographer.H>
#include <PshipClassic.H>
#include <Message.H>
#include <Person.H>
#include "adevs.h"

#warning "Using experimental partnership module"

using namespace std;

static struct options_t {
  options_t() : epidemic(""), years(60), seed(0) {}

  void print_help() {
    cout << "USAGE: simulate [options]\n\n"
	 << "Options:\n"
	 << "-h --help            show this help message and exit\n"
	 << "-e --epidemic=FILE   initialize naive epidemic parameters from FILE\n"
	 << "Simulation parameters:\n"
	 << "-c C --collect=C     memory management policy parameter\n"
	 << "-y Y --years=Y       simulate Y years (default: 60)\n"
	 << "-s SEED --seed=SEED  seed the random number generator using SEED,\n"
	 << "                     otherwise a seed is obtained from the operating system.\n"
         << "                     Ignored during deterministic simulation.\n";
  }

  void parse_args(int argc, char **argv) {
    int c(0), index(0);
    option longopts[] = {
      {"help",     no_argument,       0, 'h'},
      {"epidemic", required_argument, 0, 'e'},
      {"years",    required_argument, 0, 'y'},
      {"collect",  required_argument, 0, 'c'},
      {"seed",     required_argument, 0, 's'},
      {0, 0, 0, 0}};
    while(c >= 0) {
      c = getopt_long(argc, argv, "hc:e:y:s:", longopts, &index);
      switch(c) {
      case 'c': collect = strtoul(optarg, NULL, 0); break;
      case 'h': print_help(); exit(0); break;
      case 'e': epidemic = optarg; break;
      case 'y': years = strtoul(optarg, NULL, 0); break;
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

  unsigned int years;
  unsigned long int seed;

  unsigned int collect;

} options;

class Model : public adevs::Network<Message> {
public:
  typedef adevs::Devs<Message> Component;

  Model(unsigned int const collect = 0) : adevs::Network<Message>() {
    bookkeeper.collect(collect);
  }
  virtual ~Model() {}

  // Must be called prior to simulation
  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);

    bookkeeper.report_partnerships(true);
    bookkeeper.initialize(people);
    bookkeeper.setParent(this);

    demographer.initialize(people, gsl_rng_get(rng));
    demographer.setParent(this);

    partnership.initialize(people, gsl_rng_get(rng));
    partnership.setParent(this);

    gsl_rng_free(rng);
  }

  virtual void getComponents(adevs::Set<Component*>& c) {
    c.insert(&bookkeeper);
    c.insert(&demographer);
    c.insert(&partnership);
  }

  virtual void route(Message const& msg, Component* sender, adevs::Bag< adevs::Event<Message> >& output) {
    // fprintf(stderr, "%s[%d]:%s() %u %p\n", __FILE__, __LINE__, __FUNCTION__, msg.event_type, msg.person1);
    switch (msg.event_type) {
    case DEBUT:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&partnership, msg));
      break;
    case AGE:
      output.insert(adevs::Event<Message>(&partnership, msg));
      break;
    case PARTNER:
      break;
    case BREAKUP:
      break;
    case DEATH_NAT:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&partnership, msg));
      break;
    default: break;
    }
  }

private:
  // Modules
  Bookkeeper bookkeeper;
  Demographer demographer;
  Partnership partnership;
};

long get_vm_peak() {
  long value(0);
#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
  string token;
  stringstream stat;
  stat << "/proc/" << getpid() << "/status";
  ifstream in(stat.str().c_str());
  while (in.good()) {
    in >> token;
    if (token == "VmPeak:") {
      in >> value;
      break;
    }
  }
  in.close();
#elif defined(__APPLE__)
  // This does not actually get peak memory usage, instead it gets
  // peak physical memory usage. Not clear how to get peak memory
  // usage under OSX
  rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  value = r_usage.ru_maxrss >> 10;
#endif
  return value;
}

void simulate() {
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, options.seed);

  timeval t0, t1, t2;
  gettimeofday(&t0, NULL);

  vector<Person*> P;
  double prop;
  size_t k, m;
  std::pair<double, unsigned int> age;
  for (size_t si(0); si < SEXES; ++si) {
    for (size_t ai(0); ai < LEVELS; ++ai) {
      prop = params.prop_activity(Sex(si), Activity(ai)) * params.size_population / 2.0;
      m = static_cast<size_t>(round(prop));
      for (size_t i(0); i < m; ++i) {
  	age = params.sample_age(rng);
  	k = P.size();
  	P.push_back(new Person(-age.first));
  	P[k]->sex(Sex(si));
  	P[k]->activity(Activity(ai));
  	P[k]->age_band(age.second);
  	P[k]->seek = params.rate_partner[si][ai][P[k]->age_band()];
  	P[k]->damp = params.damp_concurrency[si][ai];
  	P[k]->circumcised((si == MALE) ? (gsl_rng_uniform(rng) < params.mmc.prop_base) : false);
  	P[k]->history.insert(new DebutRecord(-age.first));

	P[k]->stage(SUSCEPTIBLE);
	P[k]->treat(ART_NONE);
	P[k]->prep(PREP_NONE);
	P[k]->virus(VARIANTS);
      }
    }
  }

  Model M(options.collect);
  M.initialize(P, gsl_rng_get(rng));
  adevs::Simulator<Message> simulator(&M);

  gettimeofday(&t1, NULL);

  while (simulator.nextEventTime() <= options.years) simulator.execNextEvent();

  gettimeofday(&t2, NULL);

  printf("%% Simulation time: %0.4f %0.4f %0.4f %ld\n",
 	 1e-6 * (t1.tv_usec - t0.tv_usec) + (t1.tv_sec - t0.tv_sec),
 	 1e-6 * (t2.tv_usec - t1.tv_usec) + (t2.tv_sec - t1.tv_sec),	 
 	 1e-6 * (t2.tv_usec - t0.tv_usec) + (t2.tv_sec - t0.tv_sec),
	 get_vm_peak());

  gsl_rng_free(rng);
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

  params.year_done = options.years;

  if (options.epidemic.size() == 0) {
    cerr << "Error: no epidemic inputs specified, exiting." << endl;
    return -1;
  } else if (params.read_epidemic(options.epidemic) < 0) {
    cerr << "Error: failed to read epidemic inputs, exiting." << endl;
    return -1;
  } else {
    cout << "% Epidemic inputs: " << options.epidemic << '\n';
  }

  cout << flush;

  params.art.year_available[CHRONIC_LATE] = infinity;
  params.art.year_available[AIDS] = infinity;
  params.prep.year_available = infinity;

  std::fill(params.damp_concurrency[ MALE ], params.damp_concurrency[ MALE ] + LEVELS, 0.0);
  std::fill(params.damp_concurrency[FEMALE], params.damp_concurrency[FEMALE] + LEVELS, 0.0);
  for (int ki(0); ki < LEVELS; ++ki) {
    for (int kj(0); kj < LEVELS; ++kj) {
      params.rate_breakup[ki][kj] = 12.0;
    }
  }

  simulate();

  return 0;
}
