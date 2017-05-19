#ifndef ARVMODEL_H
#define ARVMODEL_H

#include <Bookkeeper.H>
#include <Demographer.H>
#include <PshipClassic.H>
//#include <Partnership.H>
#include <Progression.H>
#include <Transmission.H>
#include <Treatment.H>
#include <Prophylaxis.H>
#include <Resistance.H>
#include <Message.H>
#include <Person.H>
#include "adevs.h"

#error "ARV Model is out of date, not all events are communicated properly"

class Model : public adevs::Network<Message> {
public:
  typedef adevs::Devs<Message> Component;

  Model(unsigned int const collect = 0) : adevs::Network<Message>() {
    bookkeeper.collect(collect);
  }
  virtual ~Model() {}

  // Must be called prior to simulation
  template<typename Container>
  void initialize(Container const& people, unsigned long int seed) {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);

    bookkeeper.initialize(people);
    bookkeeper.setParent(this);

    demographer.initialize(people, gsl_rng_get(rng));
    demographer.setParent(this);

    partnership.initialize(people, gsl_rng_get(rng));
    partnership.setParent(this);

    progression.initialize(people, gsl_rng_get(rng));
    progression.setParent(this);

    transmission.initialize(people, gsl_rng_get(rng));
    transmission.setParent(this);

    treatment.initialize(people, gsl_rng_get(rng));
    treatment.setParent(this);

    prophylaxis.initialize(people, gsl_rng_get(rng));
    prophylaxis.setParent(this);

    resistance.initialize(people, gsl_rng_get(rng));
    resistance.setParent(this);

    gsl_rng_free(rng);
  }

  virtual void getComponents(adevs::Set<Component*>& c) {
    c.insert(&bookkeeper);
    c.insert(&demographer);
    c.insert(&partnership);
    c.insert(&progression);
    c.insert(&transmission);
    c.insert(&treatment);
    c.insert(&prophylaxis);
    c.insert(&resistance);
  }

  virtual void route(Message const& msg, Component* sender, adevs::Bag< adevs::Event<Message> >& output) {
    switch (msg.event_type) {
    case DEBUT:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&partnership, msg)); 
      output.insert(adevs::Event<Message>(&prophylaxis, msg));
      break;
    case AGE:
      output.insert(adevs::Event<Message>(&partnership, msg));
      break;
    case PARTNER:
      output.insert(adevs::Event<Message>(&transmission, msg));
      break;
    case BREAKUP:
      output.insert(adevs::Event<Message>(&transmission, msg));
      break;
    case DEATH_NAT:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&partnership, msg));
      output.insert(adevs::Event<Message>(&progression, msg));
      output.insert(adevs::Event<Message>(&treatment, msg));
      output.insert(adevs::Event<Message>(&resistance, msg));
      output.insert(adevs::Event<Message>(&prophylaxis, msg));
      break;
    case DEATH_HIV:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&partnership, msg));
      output.insert(adevs::Event<Message>(&demographer, msg));
      output.insert(adevs::Event<Message>(&resistance, msg));
      output.insert(adevs::Event<Message>(&prophylaxis, msg));
      if (sender == &progression) {
	output.insert(adevs::Event<Message>(&treatment, msg));
      } else {
	output.insert(adevs::Event<Message>(&progression, msg));
      }
      break;
    case INFECT:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&progression, msg));
      output.insert(adevs::Event<Message>(&resistance, msg));
      break;
    case PROGRESS:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&transmission, msg));
      output.insert(adevs::Event<Message>(&treatment, msg));
      output.insert(adevs::Event<Message>(&prophylaxis, msg));
      break;
    case ART_INIT: case ART_CHANGE:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&transmission, msg));
      output.insert(adevs::Event<Message>(&progression, msg));
      output.insert(adevs::Event<Message>(&resistance, msg));
      break;
    case PREP_INIT: case PREP_CHANGE:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&transmission, msg));
      output.insert(adevs::Event<Message>(&treatment, msg));
      output.insert(adevs::Event<Message>(&resistance, msg));
      break;
    case RESISTANCE:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&transmission, msg));
      output.insert(adevs::Event<Message>(&progression, msg));
      output.insert(adevs::Event<Message>(&treatment, msg));
      break;
    default: break;
    }
  }

private:
  // Modules
  Bookkeeper bookkeeper;
  Demographer demographer;
  Partnership partnership;
  Progression progression;
  Transmission transmission;
  Treatment treatment;
  Prophylaxis prophylaxis;
  Resistance resistance;
};

#endif // ARVMODEL_H
