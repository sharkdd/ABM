#include <cmath>
#include <Demographer.H>

Demographer::Demographer()
  : m_rng(NULL), m_tcurr(0.0), m_next(EVENTS_INT), m_debut(NULL) {
  for (size_t ei(0); ei <= EVENTS_INT; ++ei) m_time[ei] = infinity;
  // insert a sentinel record so we never have an empty queue
  m_queue.push(NULL, Record(infinity, EVENTS_INT));
}

Demographer::~Demographer() {
  if (m_rng) gsl_rng_free(m_rng);
}

void Demographer::delta_int() {
  m_tcurr = m_time[m_next];
  switch (m_next) {
  case DEBUT_INT: debut_int(); break;
  case AGING_INT: aging_int(); break;
  case DEATH_INT: death_int(); break;
  case ENTER_INT: enter_int(); break;
  default: break;
  }
  update_schedule();
}

void Demographer::delta_ext(double dt, adevs::Bag<Message> const& msgs) {
  m_tcurr += dt;
  adevs::Bag<Message>::const_iterator mi;
  for (mi = msgs.begin(); mi != msgs.end(); ++mi) {
    switch ((*mi).event_type) {
    case DEATH_NAT: death_ext((*mi).person1); break;
    case DEATH_HIV: death_ext((*mi).person1); break;
    default: break;
    }
  }
  update_schedule();
}

void Demographer::delta_conf(adevs::Bag<Message> const& msgs) {
  delta_int();
  delta_ext(0.0, msgs);
}

void Demographer::output_func(adevs::Bag<Message>& msgs) {
  switch (m_next) {
  case DEBUT_INT: debut_out(msgs); break;
  case AGING_INT: aging_out(msgs); break;
  case DEATH_INT: death_out(msgs); break;
  case ENTER_INT: enter_out(msgs); break;
  default: break;
  }
}

Demographer::Record Demographer::next_event(Person* p, bool const initialize) {
  const unsigned int band(p->age_band() + 1);
#ifdef TYPE_TWO_AGE
  // "type two aging" models aging as a stochastic exponential decay
  // process, so that the waiting time between turning age x and x+d
  // is gamma-distributed with shape d and rate 1. Module
  // initialization and simulation must be handled differently. At
  // module initialization, the average waiting time before aging
  // depends on the person's age within his or her age band. After a
  // simulated debut or aging event, the average waiting time before
  // aging is equal to the age band width (span_age_band).
  //
  // p->debut() + params.span_age_band * band - m_tcurr is the
  // (real-valued) expected time until aging. We take the ceiling here
  // so that this process emulates discrete aging in an ODE model. We
  // could accomplish the same thing by assuming the everyone in the
  // initial population was born on January 1.
  const double w_age(initialize
		     ? ceil(p->debut() + params.span_age_band * band - m_tcurr)
		     : params.span_age_band);
  const double t_age(m_tcurr + gsl_ran_gamma(m_rng, w_age, 1.0));
#else // type I aging
  const double t_age(p->debut() + params.span_age_band * band);
#endif // TYPE_TWO_AGE
  const double t_die(m_tcurr + wait_death(p));
  const double t(std::min(t_age, t_die));
  const Event e(((band == AGES) || (t_die < t_age)) ? DEATH_INT : AGING_INT); 
  return Record(t,e);
}

void Demographer::update_schedule() {
  const Record& r(m_queue.top().second);
  m_time[AGING_INT] = infinity;
  m_time[DEATH_INT] = infinity;
  m_time[r.event] = r.time;
  m_next = EVENTS_INT;
  for (unsigned int ei(0); ei < EVENTS_INT; ++ei) {
    if (m_time[ei] < m_time[m_next]) m_next = static_cast<Event>(ei);
  }
}

void Demographer::insert(Person* p, bool const initialize) {
  m_queue.push(p, next_event(p, initialize));
  //  fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, p, p->ref_num(), m_tcurr);
  p->ref_inc();
}

void Demographer::remove(Person* p) {
  if (m_queue.erase(p)) {
    // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, p, p->ref_num(), m_tcurr);
    p->ref_dec();
  }
}

void Demographer::debut_int() {
  insert(m_debut);
  m_time[DEBUT_INT] = m_tcurr + wait_debut();
  m_debut = NULL;
}

void Demographer::aging_int() {
  Person* p(m_queue.top().first);
  p->age_band(p->age_band() + 1);
  p->seek = params.rate_partner[p->sex()][p->activity()][p->age_band()];
  m_queue.push(p, next_event(p));
}

void Demographer::death_int() {
  Person* p(m_queue.top().first);
  p->history.insert(new ExitRecord(m_tcurr));
  p->alive(false);
  //  fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, p, p->ref_num(), m_tcurr);
  p->ref_dec();
  m_queue.pop();
  m_time[DEBUT_INT] = m_tcurr + wait_debut();
}

void Demographer::enter_int() {
  // Seed cases are handled as though inserted at initialization since
  // they may have arbitrary ages
  std::list<Person*>::iterator it;
  for (it = m_enter.begin(); it != m_enter.end(); ++it) {
    insert(*it, true);
    (*it)->history.insert(new InfectRecord(m_tcurr, NULL, (*it)->virus()));
  }
  m_enter.clear();
  m_time[ENTER_INT] = infinity;
}

void Demographer::death_ext(Person* person) {
  remove(person);
  m_time[DEBUT_INT] = m_tcurr + wait_debut();
}

void Demographer::debut_out(adevs::Bag<Message>& msgs) {
  m_debut = Person::build_person(m_rng, m_time[DEBUT]);
  msgs.insert(Message(DEBUT, m_debut));
}

void Demographer::aging_out(adevs::Bag<Message>& msgs) {
  Person* person(m_queue.top().first);
  msgs.insert(Message(AGE, person));
}

void Demographer::death_out(adevs::Bag<Message>& msgs) {
  Person* person(m_queue.top().first);
  msgs.insert(Message(DEATH_NAT, person));
}

void Demographer::enter_out(adevs::Bag<Message>& msgs) {
  // Seed cases are assumed to be 20 years old (and so debuted 5 years
  // prior to entry).
  const double tdebut(params.year_init - 5.0);
  const Activity ki(HIGH);
  const unsigned int ai(2); // age band 2 (ages 25-29)

  Person* qm; // male seed case
  Person* qw; // female seed case

  qm = new Person(tdebut);
  qm->sex(MALE);
  qm->activity(ki);
  qm->age_band(ai);
  qm->seek = params.rate_partner[MALE][ki][ai];
  qm->damp = params.damp_concurrency[MALE][ki];
  qm->circumcised(false);
  qm->stage(ACUTE_WINDOW);
  qm->treat(ART_NONE);
  qm->prep(PREP_NONE);
  qm->virus(WT);

  qw = new Person(tdebut);
  qw->sex(FEMALE);
  qw->activity(ki);
  qw->age_band(ai);
  qw->seek = params.rate_partner[FEMALE][ki][ai];
  qw->damp = params.damp_concurrency[FEMALE][ki];
  qw->circumcised(false);
  qw->stage(ACUTE_WINDOW);
  qw->treat(ART_NONE);
  qw->prep(PREP_NONE);
  qw->virus(WT);

  qm->history.insert(new DebutRecord(tdebut));
  qw->history.insert(new DebutRecord(tdebut));

  m_enter.insert(m_enter.begin(), qm);
  m_enter.insert(m_enter.begin(), qw);

  msgs.insert(Message(DEBUT, qm));
  msgs.insert(Message(DEBUT, qw));
  msgs.insert(Message(INFECT, NULL, qm));
  msgs.insert(Message(INFECT, NULL, qw));
}
