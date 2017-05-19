#ifndef RESISTANCE_H
#define RESISTANCE_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <DynamicQueue.H>
#include <Message.H>
#include <ModelTypes.H>
#include <Person.H>
#include "adevs.h"

// Implements a module tracking resistance state
//
// Listens for DEATH_NAT, DEATH_HIV, INFECT, ART_* and PREP_*
// events
//   DEATH_NAT, DEATH_HIV: unschedule any resistance-related events and
//   stop tracking the individual
//
//   INFECT: schedule resistance acquisition on PrEP or reversion with
//   transmitted resistance off PrEP
//
//   PREP_INIT, PREP_CHANGE: schedule resistance acquisition if on
//   PrEP or resistance reversion if PrEP pressure removed
//
//   ART_INIT, ART_CHANGE: schedule resistance acquisition if on ART
//   or resistance reversion of ART pressure removed
//
// Emits RESISTANCE events. RESISTANCE messages indicate a change in
// resistance status (emergence or reversion of resistance) that may
// impact transmission or progression rates. These events' messages
// don't store any information in the person2 field.
class Resistance : public adevs::Atomic<Message> {
public:
  Resistance();
  virtual ~Resistance();

  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    m_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_rng, seed);
    typename Container::const_iterator pi;
    for (pi = people.begin(); pi != people.end(); ++pi) {
      // Potentially insert all infected individuals. The module may
      // not track some of these individuals, for example, if they
      // have wild-type virus and are not taking ARVs
      if ((*pi)->infected()) insert(*pi);
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
  enum Event {EMERGE=0, REVERT=1, RESIST_EVENTS=2};

  typedef DynamicQueue<Person*, double, std::greater<double> > Queue;  

  // helper function that returns the next event time and populates
  // m_enext as a side-effect
  double tnext();
  // helper function that returns the next person experiencing an event
  inline Person* pnext() const {return m_queue[m_enext].top().first;}

  // helper functions that insert or remove people from tracking
  void insert(Person* p);
  void remove(Person* p);

  // Internal event handlers
  void emerge_int(Person* p);
  void revert_int(Person* p);

  // External event handlers
  inline void death_ext(Person* p) {remove(p);}
  inline void infect_ext(Person* p) {insert(p);}
  void treat_ext(Person* p);
  void prep_ext(Person* p);

  double wait_revert(Person* p);
  double wait_emerge(Person* p);

  gsl_rng* m_rng;

  double m_tcurr; // internal clock
  double m_tnext; // next event time
  Person* m_pnext; // caches the next person to change resistance state 
  Event m_enext; // cache the next event type

  // Internal event queues
  Queue m_queue[RESIST_EVENTS];

  // Lookup table that stores HIV variants that remain after reversion
  //
  // Only virus with majority resistance can revert, but as above it
  // is cheaper and more flexible to use a lookup table compared to
  // coding rules
  Virus m_revert[VARIANTS];
};

#endif // RESISTANCE_H
