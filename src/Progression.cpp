#include <Progression.H>

Stage Progression::m_successor[] = {
  STAGES,        // susceptible: no progression, so no valid successor
  ACUTE_DETECT,  // acute infection, window period
  CHRONIC_EARLY, // acute infection with detectable infection
  CHRONIC_LATE,  // early chronic infection
  AIDS,          // late chronic infection
  STAGES};       // AIDS precedes death, so no valid successor stage

Progression::Progression()
  : m_rng(NULL), m_tcurr(0.0), m_tnext(infinity), m_pnext(NULL) {
  // insert a sentinel record so that we never have to deal with an
  // empty queue
  m_people.push(NULL, infinity);
}

Progression::~Progression() {
  if (m_rng) gsl_rng_free(m_rng);
}

void Progression::delta_int() {
  m_tcurr = m_tnext;
  if (m_pnext->stage() == AIDS) {
    m_pnext->alive(false);
    // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, m_pnext, m_pnext->ref_num(), m_tcurr);
    m_pnext->ref_dec();
    m_people.pop();
    m_pnext->history.insert(new HIVDeathRecord(m_tcurr));
  } else {
    double dt;
    m_pnext->stage(m_successor[m_pnext->stage()]);
    dt = wait_progress(m_pnext);
    m_people.push(m_pnext, m_tcurr + dt);
    m_pnext->history.insert(new ProgressRecord(m_tcurr, m_pnext->stage()));
  }
  m_tnext = tnext();
}

void Progression::delta_ext(double dt, adevs::Bag<Message> const& msgs) {
  m_tcurr += dt;
  adevs::Bag<Message>::const_iterator mi;
  for (mi = msgs.begin(); mi != msgs.end(); ++mi) {
    switch ((*mi).event_type) {
    case DEATH_NAT: death_ext((*mi).person1); break;
    case DEATH_HIV: death_ext((*mi).person1); break;
    case INFECT: infect_ext((*mi).person2); break;
    case ART_INIT: treat_ext((*mi).person1); break;
    case ART_CHANGE: treat_ext((*mi).person1); break;
    case RESISTANCE: resist_ext((*mi).person1); break;
    default: break;
    }
  }
  m_tnext = tnext();
}

void Progression::delta_conf(adevs::Bag<Message> const& msgs) {
  delta_int();
  delta_ext(0.0, msgs);
}

void Progression::output_func(adevs::Bag<Message>& msgs) {
  m_pnext = pnext();
  msgs.insert(Message((m_pnext->stage() == AIDS) ? DEATH_HIV : PROGRESS, m_pnext));
}

double Progression::tnext() {
  return m_people.top().second;
}

Person* Progression::pnext() {
  return m_people.top().first;
}

void Progression::insert(Person* p) {
  m_people.push(p, m_tcurr + wait_progress(p));
  //  fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, p, p->ref_num(), m_tcurr);
  p->ref_inc();
}

void Progression::remove(Person* p) {
  if (p->infected()) {
    m_people.erase(p);
    // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, p, p->ref_num(), m_tcurr);
    p->ref_dec();
  }
}

void Progression::death_ext(Person* p) {
  remove(p);
}

void Progression::infect_ext(Person* p) {
  insert(p);
}

void Progression::treat_ext(Person* p) {
  m_people.push(p, m_tcurr + wait_progress(p));
}

void Progression::resist_ext(Person* p) {
  m_people.push(p, m_tcurr + wait_progress(p));
}

double Progression::wait_progress(Person* p) {
  return gsl_ran_exponential(m_rng, params.wait_progress(p));
}
