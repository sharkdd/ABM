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
// dissolution. Simulates a blend of random mixing and preferentialy
// mixing by age and/or sexual activity level. In contrast to
// Partnership.*, this module exactly implements Garnett and
// Anderson-style partnership formation rates (garnett1993ptrsb).
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

  // The Match structure is used to store pending partnerships and to
  // schedule dissolution of existing partnerships. Match::first is
  // male, Match::second is female
  typedef std::pair<Person*, Person*> Match;
  typedef DynamicQueue<Match, double, std::greater<double> > Queue;
  typedef __gnu_cxx::slist<Person*> EventCache;

  // Data types used to track last-known information about
  // individuals. This is stored within the module because we do not
  // know an individual's sexual activity level prior to RISK_CHANGE
  // events in particular; the age band is also tracked to avoid
  // problems if a person simultaneously ages and changes sexual
  // activity levels.
  struct Record {
    Record(unsigned int const ai = AGES, unsigned int const ki = LEVELS)
      : band(ai), activity(ki) {}
    unsigned int band     : 4;
    unsigned int activity : 3;
  };
  typedef std::tr1::unordered_map<Person*,Record> Registry;

  // Helper functions for inserting and deleting people from the
  // population. These are responsible for adjusting bookkeeping data
  // regarding formation rates, but do not modify the partnership
  // network. The event handlers (<event>_int, <event>_ext) are
  // responsible for maintaining the network.
  void insert(Person* const person);
  void remove(Person* const person);

  // Helper function to update the seek rates assocated with a person
  void update(Person* const person);

  // Helper function used to reassign a person to a different risk
  // group after aging and/or changing sexual activity levels
  void reassign(Person* const person);

  // convenience function to get a person's partnership seek rate
  static double seek(Person const* const person);

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

  // Helper function. Uses current partnership supply (m_supply below)
  // to calculate mixing terms. These are stored in m_mix.
  void calculate_mixing();
  void calculate_mixing_new();

  gsl_rng* m_rng;

  double m_tcurr; // Partnership clock
  Event m_next;   // next (internal) event type
  Match m_mnext;  // next partnership

  // Storage for individuals and their seek rates
  DynamicSampler<Person*> m_seek[SEXES][AGES][LEVELS];

  // m_mix[ma][mk][wa][wk] stores the product of mixing coefficients
  // men in age band ma and activity level mk and women in age band wa
  // and activity level wk
  double m_mix[AGES][LEVELS][AGES][LEVELS];

  // Store the rates of partnership formation between men and women in
  // each risk group. Uses the same indexing scheme as m_mix above
  double m_rate[AGES][LEVELS][AGES][LEVELS];

  // m_supply[gi][ai][ki] tracks the total supply of partnerships from
  // sex gi, age ai and activity level ki individuals
  double m_supply[SEXES][AGES][LEVELS];
  double m_supply_sex[SEXES]; // sum across age and level
  double m_supply_age[SEXES][AGES]; // sum across levels

  // A priority queue maintaining partnership dissolution times.
  Queue m_breakup;

  // maintains times of scheduled internal transitions
  double m_time[EVENTS_INT+1];

  // The caches below store individuals affected by external events
  // until they are processed by internal event handlers.
  EventCache m_cache[EVENTS_INT];

  // Map each person to his or her last-known age and sexual activity
  // level
  Registry m_registry;

  void consistent() const;

};

#endif // PARTNERSHIP_H
