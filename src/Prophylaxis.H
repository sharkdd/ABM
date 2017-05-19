#ifndef PROPHYLAXIS_H
#define PROPHYLAXIS_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <DynamicQueue.H>
#include <Message.H>
#include <ModelTypes.H>
#include <Person.H>
#include "adevs.h"

// Implements a module for pre-exposure prophylaxis (PrEP) simulation
// 
// Listens for DEATH_NAT, DEATH_HIV and PROGRESS events
//   DEATH_NAT, DEATH_HIV: unschedule any PrEP events scheduled
//   for the deceased individual
//
//   PROGRESS: at progression to detectable infection an individual
//   may be removed from PrEP by surveillance within the PrEP program
//
// Emits PREP_INIT and PREP_CHANGE events.  PREP_CHANGE events
// describe changes in PrEP status, including loss to follow-up and
// program completion, that may impact HIV transmission rates.  These
// events' messages don't store any information in person2 field.
class Prophylaxis : public adevs::Atomic<Message> {
public:
  Prophylaxis();
  virtual ~Prophylaxis();

  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    m_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_rng, seed);
    typename Container::const_iterator pi;
    for (pi = people.begin(); pi != people.end(); ++pi) {
      // Only individuals without detectable HIV infection are
      // eligible to initiate PrEP
      if ((*pi)->stage() < ACUTE_DETECT) insert(*pi);
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

  // The prophylaxis module manages these three internal
  // events. Method implementations may rely on this enum's ordering.
  //
  // FINISH   - complete the course of injections and leave PrEP
  // INJECT   - scheduled time of injection
  // INITIATE - PrEP initiation
  enum Event {FINISH=0, INJECT=1, INITIATE=2, PREP_EVENTS=3};

  typedef DynamicQueue<Person*, double, std::greater<double> > Queue;

  // helper function that returns the next event time and populates
  // m_enext as a side-effect
  double tnext();
  // helper function that returns the next person experiencing an event
  inline Person* pnext() const {return m_queue[m_enext].top().first;}

  // helper functions that add or remove people from tracking
  void insert(Person* p);
  void remove(Person* p);

  // internal event handlers
  void inject_int(Person* p);
  void finish_int(Person* p);
  void initiate_int(Person* p);

  // external event handlers
  inline void death_ext(Person* p) {remove(p);}
  void debut_ext(Person* p);
  void progress_ext(Person* p);

  // Get an individual's time of PrEP initiation
  double time_initiate(Person* p);

  // Get an individual's intended duration of PrEP use
  double wait_finish(Person* p);

  // Get the time until the next injection
  double wait_inject(Person* p);

  gsl_rng* m_rng;

  double m_tcurr; // internal clock
  double m_tnext; // next event time
  Person* m_pnext; // cache the next individual to change PrEP state
  Event m_enext; // cache the next event type

  PrEP m_successor[PREP_STATES]; // successor state when injection received

  // Internal event queues
  Queue m_queue[PREP_EVENTS];

  // Used to determine the external event type to emit for each
  // internal event type
  static EventType output[PREP_EVENTS];
};

#endif // PROPHYLAXIS_H
