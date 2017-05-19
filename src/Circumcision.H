#ifndef CIRCUMCISION_H
#define CIRCUMCISION_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <tr1/unordered_map>
#include <vector>
#include <Message.H>
#include <ModelTypes.H>
#include <Person.H>
#include "adevs.h"

// Implements a module for male medical circumcision (MMC) uptake
// among sexually active men.
// 
// Listens for DEBUT, DEATH_NAT and DEATH_HIV events
//
//   DEBUT, DEATH_NAT and DEATH_HIV may change MMC coverage. The
//   MMC uptake rate may need re-evaluated when this happens.
//
// Emits MMC events when a sexually active man is newly circumcised.
class Circumcision : public adevs::Atomic<Message> {
public:
  Circumcision();
  virtual ~Circumcision();

  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    m_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_rng, seed);
    typename Container::const_iterator pi;
    for (pi = people.begin(); pi != people.end(); ++pi) {
      if ((*pi)->sex() == MALE) insert(*pi);
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

  // Data structure that maps each person to a sequential index
  typedef std::tr1::unordered_map<Person*, unsigned int> Lookup;

  // helper function that returns the next event time and populates
  // m_enext as a side-effect
  double tnext();

  // helper function that returns the next person experiencing an event
  Person* pnext() const;

  // helper functions that add or remove people from tracking
  void insert(Person* p);
  void remove(Person* p);

  void circumcise_int(Person* p);

  inline void death_ext(Person* p) {if (p->sex() == MALE) {remove(p);}}
  inline void debut_ext(Person* p) {if (p->sex() == MALE) {insert(p);}}

  // Sample a waiting time until the next circumcision event
  double time_circumcise();

  gsl_rng* m_rng;

  double m_tcurr; // internal clock
  double m_tnext; // next event time
  Person* m_pnext; // cache the next individual to change PrEP state

  // Storage for uncircumcised men. m_lookup uses a hash table to
  // associate uncircumcised men (Person*) with sequential indices
  // into m_sampler. m_sampler stores pointers to the value pairs
  // stored by m_lookup. This pair of data structures lets us support
  // O(1) amortized lookup O(1) uniform random sampling of
  // uncircumcised men
  Lookup m_lookup;
  std::vector<Lookup::value_type*> m_sampler;

  // number of circumcised men
  size_t m_men_mmc;
};

#endif // CIRCUMCISION_H
