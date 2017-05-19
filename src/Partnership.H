#ifndef PARTNERSHIP_H
#define PARTNERSHIP_H

#include <fstream>
#include <list>
#include <vector>
#include <tr1/unordered_map>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <DynamicQueue.H>
#include <QuantizedSampler.H>
#include <ModelTypes.H>
#include <Person.H>
#include <Message.H>
#include "adevs.h"

using adevs::Atomic;
using adevs::Bag;

// Implements a module simulating partnership formation and
// dissolution. Implemented using the GERMS formalism to simulate a
// blend of random mixing and preferentialy mixing by age and/or
// sexual activity level.
//
// Listens for DEBUT, DEATH_NAT, DEATH_HIV, AGE, PROGRESS, ART_INIT,
// ART_CHANGE and RISK_CHANGE events
//
//   DEBUT: register the new individual to form partnerships.
//
//   DEATH_NAT, DEATH_HIV: end active partnerships of the deceased
//   individual and unregister him or her for new partnership events.
//
//   AGE, RISK_CHANGE: these events cause the affected individual's
//   preferred partner populations to change and generally change his
//   or her partner seeking rate.
//
//   PROGRESS, ART_INIT, ART_CHANGE: individuals may change their
//   sexual behavior due to HIV status awareness or more advanced
//   disease
//
// Emits PARTNER and BREAKUP events. To reduce simulation overheads,
// these events are only emitted when serodiscordant partnerships
// begin or end.
class Partnership : public Atomic<Message> {
public:
  Partnership();
  virtual ~Partnership();

  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    m_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_rng, seed);

    typename Container::const_iterator pi;
    for (pi = people.begin(); pi != people.end(); ++pi) insert(*pi);

    m_time[PARTNER_INT] = tnext();
    m_next = PARTNER_INT;
  }

  // adevs::Atomic<T> interface
  virtual double ta() {return m_time[m_next] - m_tcurr;}
  virtual void delta_int();
  virtual void delta_ext(double dt, Bag<Message> const& msgs);
  virtual void delta_conf(Bag<Message> const& msgs);
  virtual void output_func(Bag<Message>& msgs);
  virtual void gc_output(Bag<Message>& garbage) {}

  // Debugging interface
  double rate();

private:
  // Internal events of the partnership module. Other than PARTNER_INT
  // and BREAKUP_INT, these are scheduled immediately when
  // corresponding external events occur
  enum Event {PARTNER_INT=0,  // partnership forms
	      BREAKUP_INT=1,  // partnership dissolves
	      DEBUT_INT=2,    // new sexual debut
	      DEATH_INT=3,    // individual has died
	      AGE_INT=4,      // individual transitions to a higher age band
	      RISK_INT=5,     // individual changes sexual activity level
	      SEEK_INT=6,     // individual changes partner seek rate
	      EVENTS_INT=7};

  // Local-scope named itegral constants. The definition of m_bin
  // explains the expression for BINS
  enum Constant {
    // index of random mixing bin
    BIN_RAND = 0,                     

    // index of first bin for same-age mixing
    BIN_SAGE = 1,

    // index of first bin for disparate-age mixing
    BIN_DAGE = 1 + AGES,

    // index of first bin for same-risk (e.g., same sexual activity level) mixing
    BIN_ACT = 1 + AGES + (AGES - 1),

    // index of first bin for same-age plus same-risk mixing
    BIN_SAGE_ACT = 1 + AGES + (AGES - 1) + LEVELS,

    // index of first bin for disparate-age plus same-risk mixing
    BIN_DAGE_ACT = 1 + AGES + (AGES - 1) + LEVELS + LEVELS * AGES,

    // total number of bins
    BINS = 1 + AGES + (AGES - 1) + LEVELS + LEVELS * AGES + LEVELS * (AGES - 1)
  };

  // The Match structure is used to store pending partnerships and to
  // schedule dissolution of existing partnerships. Match::first is
  // male, Match::second is female
  typedef std::pair<Person*, Person*> Match;
  typedef DynamicQueue<Match, double, std::greater<double> > Queue;
  typedef __gnu_cxx::slist<Person*> EventCache;

  // Bins to support assortative mixing
  struct Bin {
    // Default constructor
    Bin();

    // Default destructor
    ~Bin();

    // insert a person into the bin
    void insert(Person* const person, double const seek_rate);

    // remove a person from the bin
    void remove(Person* const person);

    // update a person's seek rate
    void update(Person* const person, double const seek_rate);

    // Check if a person mixes in the bin
    bool member(Person* const person);

    // return the rate of partnership formation in the bin
    double rate();

    // return a potential match from the bin. Pass in the random
    // number generator used for sampling.
    Match sample(gsl_rng* rng) const;

    // stores concurrency-damped partnership seek rates for each person
    DynamicSampler<Person*> seek[SEXES];
  private:
    // caches the rate of partnership formation in the bin
    double m_rate;

    // This flag indicates whether the cached partnership formation
    // rate is valid
    bool m_valid;
  };

  // Helper functions for inserting and deleting people from the
  // population. These are responsible for adjusting bookkeeping data
  // regarding formation rates, but do not modify the partnership
  // network. The event handlers (<event>_int, <event>_ext) are
  // responsible for maintaining the network.
  void insert(Person* const person);
  void remove(Person* const person);

  // Helper function to update the seek rates assocated with a person
  // in each bin they mix in.
  void update(Person* const person);

  // convenience function to get a person's partnership seek rate
  static double seek(Person const* const person);

  // convenience function that determines the set of indices for bins
  // a given person forms partnerships in, and the partitioning
  // weights used to allocate his or her seek rate across bins. index
  // must have space for six integers, weight must have space for six
  // doubles.
  //
  // On return, n and index will be set such that the given person p
  // mixes in bins Bin[index[k]] for k=0..n-1, and the seek rate in
  // Bin[index[k]] = weight[k] * seek(p)
  static void bins(Person const* const person, int *index, double *weight, int &n);

  // Handlers for internal events.
  void age_int();
  void partner_int();
  void breakup_int();
  void debut_int();
  void death_int();
  void risk_int();
  void seek_int();

  // Handlers that produce output upon internal events.
  void partner_out(Bag<Message>& msgs);
  void breakup_out(Bag<Message>& msgs);
  void death_out(Bag<Message>& msgs);

  // Schedule internal events to update the module state in response
  // to external events
  void schedule_ext(Person* p, Event e);

  // tnext() returns the time of the next partnership formation event
  double tnext();

  // mnext() returns a Match object storing the male and female
  // involved in the next partnership
  Match mnext() const;

  gsl_rng* m_rng;

  double m_tcurr; // Partnership clock
  Event m_next;   // next (internal) event type
  Match m_mnext;  // next partnership

  // 1. One bin for age-agnostic, activity level-agnostic mixing
  // 2. One bin for same-age-assortative, level-agnostic mixing for
  //    each age band
  // 3. One bin for disparate-age-assortative, level-agnostic mixing
  //    for each age band except the youngest men and oldest women
  // 4. One bin for level-assortative, age-agnostic mixing for each
  //    sexual activity level
  // 5. One bin for same-age-assortative, level-assortative mixing for
  //    each pair of ages and activity levels
  // 6. One bin for disparate-age-assortative, level-assortative
  //    mixing for each pair of ages and activity levels, excluding
  //    the youngest men and oldest women
  Bin m_bin[BINS];

  // caches the rates of partnership formation within each bin
  // between calls to tnext()
  double m_rate[BINS];

  // A priority queue maintaining partnership dissolution times.
  Queue m_breakup;

  // maintains times of scheduled internal transitions
  double m_time[EVENTS_INT+1];

  // The caches below store individuals affected by external events
  // until they are processed by internal event handlers.
  EventCache m_cache[EVENTS_INT];
};

#endif // PARTNERSHIP_H
