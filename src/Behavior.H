#ifndef BEHAVIOR_H
#define BEHAVIOR_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <QuantizedSampler.H>
#include <Message.H>
#include <ModelTypes.H>
#include <Person.H>
#include "adevs.h"

// Implements a module simulating behavior change over time,
// implemented as transitions between sexual activity levels. We
// assume women cease sex work according to a fixed rate. Women in
// lower sexual activity levels initiate sex work at a dynamic rate
// intended to maintain a constant proportion of women engaged in sex
// work
//
// Listens for DEBUT, DEATH_NAT, PROGRESS, ART_INIT, ART_CHANGE and
// DEATH_HIV events.
//
//   DEBUT, DEATH_NAT, DEATH_HIV: reschedule next sex work initiation
//   event, as these events may cause the proportion of women engaged
//   in sex work to deviate from its target.
//
//   PROGRESS, ART_INIT, ART_CHANGE: recalculate the propensity of
//   women to initiate or cease sex work (women with unsuppressed HIV
//   in the AIDS stage are less likely to initiate and more likely to
//   cease sex work).
//
// Emits RISK_CHANGE events. RISK_CHANGE messages indicate a change in
// partner change rate and the type of partnerships women engage in.
// These events' messages don't store any information in the person2
// field.

class Behavior : public adevs::Atomic<Message> {
public:
  Behavior();
  virtual ~Behavior();

  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    m_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_rng, seed);
    typename Container::const_iterator pi;
    for (pi = people.begin(); pi != people.end(); ++pi) {
      // This module only supports changes in sexual behavior level
      // for women
      if ((*pi)->sex() == FEMALE) insert(*pi);
    }
    m_tnext = tnext();
  }

  // adevs::Atomic<T> interface
  virtual double ta() {return m_tnext - m_tcurr;}
  virtual void delta_int();
  virtual void delta_ext(double dt, adevs::Bag<Message> const& msgs);
  virtual void delta_conf(adevs::Bag<Message> const& msgs);
  virtual void output_func(adevs::Bag<Message>& msgs);
  virtual void gc_output(adevs::Bag<Message>& garbage) {}

private:
  // Internal events
  // CSW_INIT_INT: woman initiates sex work
  // CSW_EXIT_INT: woman ceases sex work
  enum Event {CSW_INIT_INT=0, CSW_EXIT_INT=1, EVENTS_INT=2};

  // Helper function that returns the next event time. Sets m_enext
  // as a side-effect
  double tnext();

  // Helper function that returns the person experiencing the next
  // event
  Person* pnext() const;

  // helper functions that insert or remove people from tracking
  void insert(Person* const p);
  void remove(Person* const p);
  void update(Person* const p);

  // Internal event handlers
  void csw_init_int(Person* const p);
  void csw_exit_int(Person* const p);

  // External event handlers
  void debut_ext(Person* const p);
  void death_ext(Person* const p);
  void progress_ext(Person* const p);
  void treat_ext(Person* const p);

  double wait_csw_init();
  inline double wait_csw_exit() {
    return gsl_ran_exponential(m_rng, 1.0 / m_women_csw.sum());
  }

  // Helper function: per-capita rate of sex-work initiation
  inline double rate_init(Person* const person) const {
    return params.prop_risk_init[person->stage()][person->treat()];
  }

  inline double rate_exit(Person* const person) const {
    return params.prop_risk_exit[person->stage()][person->treat()] * params.rate_csw_exit;
  }

  gsl_rng* m_rng;

  double m_tcurr; // internal clock
  double m_tnext; // next event time
  Person* m_pnext; // person affected by next event
  Event m_enext; // cache the next event type

  // Associates each woman who is not a sex worker with her relative
  // propensity to initate sex work.
  DynamicSampler<Person*> m_women_non;

  // Associates each female sex worker with the rate she ceases sex
  // work.
  DynamicSampler<Person*> m_women_csw;
};

#endif // BEHAVIOR_H
