#ifndef TREATMENT_H
#define TREATMENT_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <DynamicQueue.H>
#include <Message.H>
#include <ModelTypes.H>
#include <Person.H>
#include "adevs.h"

// Implements a single treatment regimen module
//
// Listens for DEATH_NAT, DEATH_HIV, PREP_*, PROGRESS and
// RESISTANCE events
//
//   DEATH_NAT, DEATH_HIV: unschedule any treatment-related events
//   scheduled for the deceased individual
//
//   PROGRESS, RESISTANCE: update event rates based on the
//   individual's new state
//
//   PREP_*: infected individuals on PrEP may not start ART until
//   their infection is detected
//
// Emits ART_INIT, ART_CHANGE and DEATH_HIV events. ART events change
// in treatment status that may impact transmission and progression
// rates. These events' messages don't store any information in the
// person2 field.
class Treatment : public adevs::Atomic<Message> {
public:
  Treatment();
  virtual ~Treatment();

  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    m_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_rng, seed);
    typename Container::const_iterator pi;
    for (pi = people.begin(); pi != people.end(); ++pi) {
      if ((*pi)->detectable()) insert(*pi);
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

  // The treatment module manages several internal events. Method
  // implementations may rely on this enum's ordering. It is intended
  // to speed up event scheduling.
  //
  // DROPOUT   - loss to follow-up
  // FAILURE   - treatment failure
  // MORTALITY - death attributable to progression or immune reconstitution
  // ADVANCE   - transition from early to late ART
  // INITIATE  - ART initiation among cases eligible for treatment
  enum Event {DROPOUT=0, FAILURE=1, MORTALITY=2, ADVANCE=3, INITIATE=4, ART_EVENTS=5};

  // Data structure used to schedule treatment-related events
  struct Record {
    Record(double t = infinity, Event e = ART_EVENTS) : time(t), event(e) {}
    bool operator<(Record const& record) const {return time > record.time;}
    double time;
    Event event;
  };

  typedef DynamicQueue<Person*, Record> Queue;

  // helper function that returns the next event time and populates
  // m_enext as a side-effect
  double tnext();
  // helper function that returns the next person experiencing an event
  inline Person* pnext() const {return m_queue.top().first;}

  // helper functions that add or remove people from tracking
  void insert(Person* p);
  void remove(Person* p);

  // External event handlers
  inline void death_ext(Person* p) {remove(p);}
  void progress_ext(Person* p);
  void resist_ext(Person* p);
  void prep_ext(Person* p);

  // Get the ART initiation time for an individual
  double time_initiate(Person* p);

  // Get the next event type e and time t for a person receiving ART
  void next_event(Person* p, Event &e, double &t);

  // Return true if a person has majority acquired resistance to ART,
  // false otherwise
  static inline bool art_resistant(Person* p) {
    return p->virus() == R2 || p->virus() == C2;
  } 

  // Return true if a person has ever acquired resistance to ART
  static inline bool art_failure(Person* p) {
    return art_resistant(p) || p->virus() == WR2 || p->virus() == WC2;
  }

  // Calculate the rate of ART-related mortality for person p.
  inline double rate_mortality(Person* p) {
    // HIV disease progression and mortality are handled by the
    // progression module for persons with acquired resistance to ART
    return (((p->virus() != R2) && (p->virus() != C2))
	    * params.art.rate_mortality[p->treat()] 
	    * params.art.prop_mortality[p->stage()]);
  }

  gsl_rng* m_rng;

  // Module internal clock
  double m_tcurr;

  // Next event time
  double m_tnext;

  // Cache that stores the next person scheduled to experience an
  // ART-related event
  Person* m_pnext;

  // Cache that maintains the next event type
  Event m_enext;

  // Internal event queue
  Queue m_queue;

  // Used to store the resulting ART state of each ART event.
  Treat successor[ART_EVENTS];

  // Used to determine the external event type to emit for each
  // internal event type
  static EventType output[ART_EVENTS];
};

#endif // TREATMENT_H
