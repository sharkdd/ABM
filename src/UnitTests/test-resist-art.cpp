#include <getopt.h>
#include <sys/time.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <Treatment.H>
#include <Resistance.H>
#include <Progression.H>
#include "adevs.h"

using namespace std;

static struct options_t {
  options_t() : epidemic(""), analysis(""), years(100), seed(0), art(0), prep(0) {}

  void print_help() {
    cout << "USAGE: simulate [options]\n\n"
	 << "Options:\n"
	 << "-h --help            show this help message and exit\n"
	 << "-e --epidemic=FILE   initialize naive epidemic parameters from FILE\n"
	 << "-a --analysis=FILE   initialize ARV-related parameters from FILE (ignored)\n"
	 << "Interventions: simulate an ARV-naive epidemic by default\n"
	 << "--art                include ART rollout\n"
	 << "--prep               include PrEP rollout\n"
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
      {"analysis", required_argument, 0, 'a'},
      {"years",    required_argument, 0, 'y'},
      {"collect",  required_argument, 0, 'c'},
      {"art",      no_argument,       &art, 1},
      {"prep",     no_argument,       &prep, 1},
      {"seed",     required_argument, 0, 's'},
      {0, 0, 0, 0}};
    while(c >= 0) {
      c = getopt_long(argc, argv, "hc:e:a:y:s:", longopts, &index);
      switch(c) {
      case 'c': collect = strtoul(optarg, NULL, 0); break;
      case 'h': print_help(); exit(0); break;
      case 'e': epidemic = optarg; break;
      case 'a': analysis = optarg; break;
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
  string analysis; // intervention inputs

  unsigned int years;
  unsigned long int seed;

  int art; // true if ART rollout is enabled, false otherwise
  int prep; // true if PrEP rollout is enabled, false otherwise

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

    treatment.initialize(people, gsl_rng_get(rng));
    treatment.setParent(this);

    progression.initialize(people, gsl_rng_get(rng));
    progression.setParent(this);

    resistance.initialize(people, gsl_rng_get(rng));
    resistance.setParent(this);

    gsl_rng_free(rng);
  }

  // adevs::Network interface method
  virtual void getComponents(adevs::Set<Component*>& c) {
    c.insert(&treatment);
    c.insert(&progression);
    c.insert(&resistance);
  }

  // adevs::Network interface method
  virtual void route(Message const& msg, Component* sender, adevs::Bag< adevs::Event<Message> >& output) {
    switch (msg.event_type) {
    case PROGRESS:
      output.insert(adevs::Event<Message>(&treatment, msg));
      break;
    case ART_INIT: case ART_CHANGE:
      output.insert(adevs::Event<Message>(&progression, msg));
      output.insert(adevs::Event<Message>(&resistance, msg));
      break;
    case RESISTANCE:
      output.insert(adevs::Event<Message>(&progression, msg));
      output.insert(adevs::Event<Message>(&treatment, msg));
      break;
    case DEATH_HIV:
      if (sender == &progression) {
	output.insert(adevs::Event<Message>(&treatment, msg));
      } else {
	output.insert(adevs::Event<Message>(&progression, msg));
      }
      output.insert(adevs::Event<Message>(&resistance, msg));
      break;
    default: break;
    }
  }

private:
  // Modules
  Treatment treatment;
  Resistance resistance;
  Progression progression;
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

  printf("person sex t.curr t.next e.next h.curr u.curr v.curr h.next u.next v.next\n");

  Stage scurr, snext;
  Virus vcurr, vnext;
  Treat ucurr, unext;
  for (size_t i(0); i < P.size(); ++i) {
    Person::History reversed;    
    Person::History::iterator hcurr, hnext;
    for (hcurr = P[i]->history.begin(); hcurr != P[i]->history.end(); ++hcurr) reversed.insert(*hcurr);

    snext = scurr = CHRONIC_LATE;
    vnext = vcurr = VARIANTS;
    unext = ucurr = ART_EARLY;

    assert(reversed.size() > 1);

    hnext = hcurr = reversed.begin();
    for (++hnext; hnext != reversed.end(); ++hnext, ++hcurr) {

      if ((*hnext)->type() == INFECT) {
	vnext = dynamic_cast<InfectRecord*>(*hnext)->virus;
      } else if (((*hnext)->type() == ART_INIT) || ((*hnext)->type() == ART_CHANGE)) {
	unext = dynamic_cast<TreatmentRecord*>(*hnext)->treat;
      } else if ((*hnext)->type() == PROGRESS) {
	snext = static_cast<Stage>(scurr + 1);
      } else if ((*hnext)->type() == RESISTANCE) {
	vnext = dynamic_cast<ResistanceRecord*>(*hnext)->virus;
      }

      // Print the times of consecutive events and the person's HIV
      // state after the current (curr) and next events.
      printf("%p %d %10.6f %10.6f % 3d %d %d % 3d %d %d % 3d\n",
	     P[i], P[i]->sex(), // person identity
	     (*hcurr)->time, (*hnext)->time, // consecutive event times
	     (*hnext)->type(), // event type
	     scurr, ucurr, vcurr,  // HIV status after current event
	     snext, unext, vnext); // HIV status after next event

      // update current event
      scurr = snext;
      vcurr = vnext;
      ucurr = unext;
    }
  }

  //   Stage stage;
  //   Virus virus;
  //   Treat treat;
  //   for (size_t i(0); i < P.size(); ++i) {
  //     Person::History reversed;
  //     Person::History::iterator hi;
  //     // put history in chronological order
  //     for (hi = P[i]->history.begin(); hi != P[i]->history.end(); ++hi) reversed.insert(*hi);
  //     stage = AIDS;
  //     treat = ART_EARLY;
  //     virus = VARIANTS; // should never see this in a history
  //     for (hi = reversed.begin(); hi != reversed.end(); ++hi) {    
  //       if ((*hi)->type() != DEBUT) {
  // 	// Ignore DEBUT events since they're not relevant to what
  // 	// we're trying to test and just clutter up the output
  // 	if ((*hi)->type() == INFECT) {
  // 	  virus = dynamic_cast<InfectRecord*>(*hi)->virus;
  // 	} else if (((*hi)->type() == ART_INIT) || ((*hi)->type() == ART_CHANGE)) {
  // 	  treat = dynamic_cast<TreatmentRecord*>(*hi)->treat;
  // 	} else if ((*hi)->type() == PROGRESS) {
  // 	  stage = static_cast<Stage>(stage + 1);
  // 	} else if ((*hi)->type() == RESISTANCE) {
  // 	  virus = dynamic_cast<ResistanceRecord*>(*hi)->virus;
  // 	}
  // 	printf("%p %d %10.6f %d %d %d %d\n", P[i], P[i]->sex(), (*hi)->time, (*hi)->type(), stage, treat, virus);
  //       }
  //     }
  //   }

  gettimeofday(&t2, NULL);
  rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);

  printf("%% Simulation time: %0.4f %0.4f %0.4f %ld\n",
 	 1e-6 * (t1.tv_usec - t0.tv_usec) + (t1.tv_sec - t0.tv_sec),
 	 1e-6 * (t2.tv_usec - t1.tv_usec) + (t2.tv_sec - t1.tv_sec),	 
 	 1e-6 * (t2.tv_usec - t0.tv_usec) + (t2.tv_sec - t0.tv_sec),
	 r_usage.ru_maxrss);
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

  if (options.analysis.size() == 0) {
    cerr << "Error: no intervention inputs specified, exiting." << endl;
    return -1;
  } else if (params.read_intervention(options.analysis) < 0) {
    cerr << "Error: failed to read intervention inputs, exiting." << endl;
    return -1;
  } else {
    cout << "% Intervention inputs: " << options.analysis << '\n';
  }

  cout << flush;

  if (argc < 3) {cerr << "USAGE " << argv[0] << " <N> <M>\n"; return -1;}

  const unsigned long int seed(get_seed());
  if (seed) {
    cout << "% Seed: " << seed << '\n';
  } else {
    cerr << "Fatal error: Failed to get a random seed\n";
    return -1;
  }

  const size_t N(params.size_population);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, seed);

  // assume ART scale-up has already completed
  const double prop_uptake(0.55);
  params.art.year_available[ACUTE_WINDOW ] = infinity;
  params.art.year_available[ACUTE_DETECT ] = infinity;
  params.art.year_available[CHRONIC_EARLY] = infinity;
  params.art.year_available[CHRONIC_LATE ] = -10.0;
  params.art.year_available[AIDS         ] = -10.0;
  params.art.year_scale[ACUTE_WINDOW ] = infinity;
  params.art.year_scale[ACUTE_DETECT ] = infinity;
  params.art.year_scale[CHRONIC_EARLY] = infinity;
  params.art.year_scale[CHRONIC_LATE ] = 0.0;
  params.art.year_scale[AIDS         ] = 0.0;
  params.art.prop_scale[CHRONIC_LATE] = prop_uptake;
  params.art.prop_scale[AIDS        ] = prop_uptake;
  params.art.prop_level[CHRONIC_LATE] = prop_uptake;
  params.art.prop_level[AIDS        ] = prop_uptake;
  params.art.rate_uptake_max[CHRONIC_LATE] = -log(1.0 - prop_uptake);
  params.art.rate_uptake_max[AIDS        ] = -log(1.0 - prop_uptake);

  // Build a population of N individuals with AIDS who have newly
  // initiated ART. These individuals are assumed to have WT or
  // transmitted ART resistance at initiation.
  std::vector<Person*> P(N);
  for (size_t i(0); i < N; ++i) {
    P[i] = Person::build_person(rng, 0.0);
    P[i]->stage(CHRONIC_LATE);
    P[i]->treat(ART_EARLY);
    do {
      P[i]->virus(static_cast<Virus>(gsl_rng_uniform_int(rng, VARIANTS)));
    } while ((P[i]->virus() != WT)
	     && (P[i]->virus() != R1) && (P[i]->virus() != WR1)
	     && (P[i]->virus() != C1) && (P[i]->virus() != WC1));
    P[i]->history.insert(new InfectRecord(0.0, NULL, P[i]->virus()));
  }
  simulate(P, gsl_rng_get(rng));
  for (size_t n(0); n < P.size(); ++n) delete P[n];

  gsl_rng_free(rng);

  return 0;
}
