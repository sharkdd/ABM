#ifndef TRANSMISSION_H
#define TRANSMISSION_H

#include <tr1/unordered_map>
#include <gsl/gsl_rng.h>
#include <BinaryIndexedTree.H>
#include <Message.H>
#include <ModelTypes.H>
#include <Person.H>
#include "adevs.h"

// Implements a transmission module.
//
// Listens for PARTNER, BREAKUP, DEATH_NAT, DEATH_HIV, PROGRESS,
// ART_*, PREP_* and RESISTANCE events.
//   PARTNER: monitor the newly-formed partnership if it is
//   serodiscordant, otherwise the partnership is ignored unless a
//   partner seroconverts.
//
//   BREAKUP: stop monitoring the couple breaking up
//
//   DEATH_NAT, DEATH_HIV: stop monitoring any serodiscordant
//   partnerships the deceased was involved in
//
//   PROGRESS, ART_*, PREP_*, RESISTANCE: update the transmission rate
//   in any serodiscordant partnerships the individual was involved in
//
// Emits INFECT events. A message with event_type=INFECT stores the
// donor in the person1 field and the recipient in person2
class Transmission : public adevs::Atomic<Message> {
public:
  Transmission();
  virtual ~Transmission();

  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    m_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_rng, seed);
    // Scan the initial population and store records for infected
    // people in partnerships -- infected people outside of
    // partnerships pose no transmission risk in this model
    typename Container::const_iterator pi;
    Person::Partners::const_iterator qi;
    for (pi = people.begin(); pi != people.end(); ++pi) {
      if ((*pi)->infected()) {
	for (qi = (*pi)->partners.begin(); qi != (*pi)->partners.end(); ++qi) {
	  if ((*qi)->susceptible()) insert(Couple(*pi, *qi));
	}
      }
    }
    m_tnext = tnext();
  }

  // adevs::Atomic<T> interface
  virtual double ta();
  virtual void delta_int();
  virtual void delta_ext(double dt, adevs::Bag<Message> const& msgs);
  virtual void delta_conf(adevs::Bag<Message> const& msgs);
  virtual void output_func(adevs::Bag<Message>& msgs);
  virtual void gc_output(adevs::Bag<Message>& garbage) {}

private:

  // serodiscordant partnership record. Couple.first is infected,
  // Couple.second is susceptible
  typedef std::pair<Person*, Person*> Couple;

  typedef std::tr1::unordered_map<Couple, unsigned int> Index;

  // external event handlers
  void partner_ext(Person* p1, Person* p2);
  void breakup_ext(Person* p1, Person* p2);
  inline void progress_ext(Person* p) {update(p);}
  inline void treat_ext(Person* p) {update(p);}
  void prep_ext(Person* p);
  inline void resist_ext(Person* p) {update(p);}
  void mmc_ext(Person* p);

  // update a person's infectiousness
  void update(Person* p);

  // draw the wait time until the next transmission event
  double tnext() const;

  // draw the couple among whom transmission next occurs
  Couple cnext() const;

  void insert(Couple couple);
  void remove(Couple couple);

  // return the transmission rate in a serodiscordant partnership
  double risk(Couple const& couple);

  gsl_rng* m_rng;

  // Transmission engine clock 
  double m_tcurr;

  // next transmission time
  double m_tnext;

  // people involved in the next transmission event
  Couple m_cnext;

  // Maps each serodiscordant partnership to a unique index
  Index m_couple2index;

  // maps each index to a serodiscordant partnership
  std::vector<Couple> m_index2couple;

  // stores transmission rates for each serodiscordant partnership
  BinaryIndexedTree<double> m_transmit;

  // Lookup table mapping the donor virus to the transmitted variant
  Virus m_virus[VARIANTS];
};

#endif // TRANSMISSION_H
