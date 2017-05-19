#ifndef ARTMODEL_H
#define ARTMODEL_H

#include <Behavior.H>
#include <Bookkeeper.H>
#include <Circumcision.H>
#include <Demographer.H>
#include <Progression.H>
#include <Resistance.H>
#include <Transmission.H>
#include <Treatment.H>
#include <Message.H>
#include <Person.H>
#include "adevs.h"

#ifdef GERMS
#  include <Partnership.H>
#else
#  include <PshipClassic.H>
#endif

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

    behavior.initialize(people, gsl_rng_get(rng));
    behavior.setParent(this);

    progression.initialize(people, gsl_rng_get(rng));
    progression.setParent(this);

    circumcision.initialize(people, gsl_rng_get(rng));
    circumcision.setParent(this);

    transmission.initialize(people, gsl_rng_get(rng));
    transmission.setParent(this);

    treatment.initialize(people, gsl_rng_get(rng));
    treatment.setParent(this);

    resistance.initialize(people, gsl_rng_get(rng));
    resistance.setParent(this);

    gsl_rng_free(rng);
  }

  virtual void getComponents(adevs::Set<Component*>& c) {
    c.insert(&bookkeeper);
    c.insert(&demographer);
    c.insert(&partnership);
    c.insert(&behavior);
    c.insert(&progression);
    c.insert(&circumcision);
    c.insert(&transmission);
    c.insert(&treatment);
    c.insert(&resistance);
  }

  virtual void route(Message const& msg, Component* sender, adevs::Bag< adevs::Event<Message> >& output) {
    // fprintf(stderr, "%s[%d]:%s() %u %p\n", __FILE__, __LINE__, __FUNCTION__, msg.event_type, msg.person1);
    switch (msg.event_type) {
    case DEBUT:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&partnership, msg));
      output.insert(adevs::Event<Message>(&behavior, msg));
      output.insert(adevs::Event<Message>(&circumcision, msg));
      break;
    case AGE:
      output.insert(adevs::Event<Message>(&partnership, msg));
      break;
    case RISK_CHANGE:
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
      output.insert(adevs::Event<Message>(&behavior, msg));
      output.insert(adevs::Event<Message>(&progression, msg));
      output.insert(adevs::Event<Message>(&circumcision, msg));
      output.insert(adevs::Event<Message>(&treatment, msg));
      output.insert(adevs::Event<Message>(&resistance, msg));
      break;
    case DEATH_HIV:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&partnership, msg));
      output.insert(adevs::Event<Message>(&demographer, msg));
      output.insert(adevs::Event<Message>(&behavior, msg));
      output.insert(adevs::Event<Message>(&circumcision, msg));
      output.insert(adevs::Event<Message>(&resistance, msg));
      if (sender == &progression) {
      	output.insert(adevs::Event<Message>(&treatment, msg));
      } else {
	output.insert(adevs::Event<Message>(&progression, msg));
      }
      break;
    case MMC:
      output.insert(adevs::Event<Message>(&transmission, msg));
      break;
    case INFECT:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&progression, msg));
      output.insert(adevs::Event<Message>(&resistance, msg));
      break;
    case PROGRESS:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&behavior, msg));
      output.insert(adevs::Event<Message>(&transmission, msg));
      output.insert(adevs::Event<Message>(&treatment, msg));
      break;
    case ART_INIT: case ART_CHANGE:
      output.insert(adevs::Event<Message>(&bookkeeper, msg));
      output.insert(adevs::Event<Message>(&transmission, msg));
      output.insert(adevs::Event<Message>(&partnership, msg));
      output.insert(adevs::Event<Message>(&progression, msg));
      output.insert(adevs::Event<Message>(&behavior, msg));
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
  Behavior behavior;
  Progression progression;
  Circumcision circumcision;
  Transmission transmission;
  Resistance resistance;
  Treatment treatment;
};

#endif // ARTMODEL_H
