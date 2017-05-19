#include <Person.H>

Person::Person(double const tdebut) : m_debut(tdebut), m_alive(true) {m_count = 0;}

Person::~Person() {
#ifdef HISTORY
  History::iterator hi;
  for (hi = history.begin(); hi != history.end(); ++hi) delete *hi;
  history.clear();
#endif // HISTORY
}

Person* Person::build_person(gsl_rng* rng, const double tdebut) {
  Person* person(new Person(tdebut));
  const Sex sex((gsl_rng_uniform(rng) < params.prop_female) ? FEMALE : MALE);

  person->alive(true);
  person->age_band(0);

  // Initialize behavioral attributes
  person->sex(sex);
  person->activity(params.sample_activity_level(rng, sex));
  person->seek = params.rate_partner[sex][person->activity()][0];
  person->damp = params.damp_concurrency[sex][person->activity()];

  // Determine circumcision status
  person->circumcised((sex == MALE) ? (gsl_rng_uniform(rng) < params.mmc.coverage(tdebut).first) : false);

  // Initialize epidemiological attributes
  person->stage(SUSCEPTIBLE);
  person->treat(ART_NONE);
  person->prep(PREP_NONE);
  person->virus(VARIANTS);

#ifdef HISTORY
  person->history.insert(new DebutRecord(tdebut));
#endif // HISTORY

  return person;
}
