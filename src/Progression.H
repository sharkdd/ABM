#ifndef PROGRESSION_H
#define PROGRESSION_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <DynamicQueue.H>
#include <Message.H>
#include <ModelTypes.H>
#include <Person.H>
#include "adevs.h"

// Implements a discrete-staged HIV disease progression module.
//
// Listens for DEATH_NAT, DEATH_HIV, INFECT, ART_* and RESISTANCE
// events
//   DEATH_NAT, DEATH_HIV: stop monitoring the deceased individual
//
//   INFECT: monitor the newly infected individual
//
//   ART_INIT, ART_CHANGE, RESISTANCE: update disease progression
//   rates based on the individual's new state
//
// Emits PROGRESS and DEATH_HIV events. Neither of these events stores
// information in the person2 message field.
class Progression : public adevs::Atomic<Message> {
public:
  Progression();
  virtual ~Progression();

  // Container must store pointers of type Person*.
  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    m_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_rng, seed);
    typename Container::const_iterator pi;
    for (pi = people.begin(); pi != people.end(); ++pi) {
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
  typedef DynamicQueue<Person*, double, std::greater<double> > Queue;

  double tnext();
  Person* pnext();

  void insert(Person* p);
  void remove(Person* p);

  // External event handlers
  void death_ext(Person* p);
  void infect_ext(Person* p);
  void treat_ext(Person* p);
  void resist_ext(Person* p);

  // get the time until a person's next progression event
  double wait_progress(Person* p);

  gsl_rng* m_rng;

  // Progression module clock
  double m_tcurr;

  // Next progression event time
  double m_tnext;

  // Cache that stores the next person scheduled to experience a
  // progression event
  Person* m_pnext;

  Queue m_people;

  // maps each infection stage to its immediate successor.
  static Stage m_successor[STAGES];
};

#endif // PROGRESSION_H
