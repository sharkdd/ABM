#include <Behavior.H>

Behavior::Behavior() 
  : m_rng(NULL), m_tcurr(0.0), m_pnext(NULL) {
  m_tnext = infinity;
}

Behavior::~Behavior() {
  if (m_rng) gsl_rng_free(m_rng);
}

void Behavior::delta_int() {
  m_tcurr = m_tnext;
  switch (m_enext) {
  case CSW_INIT_INT: csw_init_int(m_pnext); break;
  case CSW_EXIT_INT: csw_exit_int(m_pnext); break;
  default: assert(false); break; // should never get here
  }
  m_tnext = tnext();
}

void Behavior::delta_ext(double dt, adevs::Bag<Message> const& msgs) {
  m_tcurr += dt;
  adevs::Bag<Message>::const_iterator mi;
  for (mi = msgs.begin(); mi != msgs.end(); ++mi) {
    switch ((*mi).event_type) {
    case DEBUT: debut_ext((*mi).person1); break;
    case DEATH_NAT: death_ext((*mi).person1); break;
    case DEATH_HIV: death_ext((*mi).person1); break;
    case PROGRESS: progress_ext((*mi).person1); break;
    case ART_INIT: treat_ext((*mi).person1); break;
    case ART_CHANGE: treat_ext((*mi).person1); break;
    default: break;
    }
  }
  m_tnext = tnext();
}

void Behavior::delta_conf(adevs::Bag<Message> const& msgs) {
  delta_ext(0.0, msgs);
  delta_int();
}

void Behavior::output_func(adevs::Bag<Message>& msgs) {
  m_pnext = pnext();
  msgs.insert(Message(RISK_CHANGE, m_pnext));
}

double Behavior::tnext() {
  double winit, wexit;
  winit = wait_csw_init();
  wexit = wait_csw_exit(); 
  if (winit < wexit) {
    m_enext = CSW_INIT_INT;
    return m_tcurr + winit;
  } else {
    m_enext = CSW_EXIT_INT;
    return m_tcurr + wexit;
  }
}

Person* Behavior::pnext() const {
  if (m_enext == CSW_INIT_INT) {
    return m_women_non.sample(m_rng);
  } else if (m_enext == CSW_EXIT_INT) {
    return m_women_csw.sample(m_rng);
  } else {
    fprintf(stderr, "%s[%d]:%s unrecognized event.\n", __FILE__, __LINE__, __FUNCTION__);
    assert(false);
  }
}

void Behavior::insert(Person* const p) {
  assert(p->sex() == FEMALE);
  if (p->activity() == HIGH) {
    m_women_csw.insert(p, rate_exit(p));
  } else {
    m_women_non.insert(p, rate_init(p));
  }
  p->ref_inc();
}

void Behavior::remove(Person* const p) {
  if (p->activity() == HIGH) {
    assert(m_women_csw.member(p));
    m_women_csw.remove(p);
  } else {
    assert(m_women_non.member(p));
    m_women_non.remove(p);
  }
  p->ref_dec();
}

void Behavior::update(Person* const p) {
  assert(p->sex() == FEMALE);
  if (p->activity() == HIGH) {
    assert(m_women_csw.member(p));
    m_women_csw.update(p, rate_exit(p));
  } else {
    assert(m_women_non.member(p));
    m_women_non.update(p, rate_init(p));
  }
}

void Behavior::csw_init_int(Person* const p) {
  assert(m_women_non.member(p));
  assert(p->activity() != HIGH);
  p->activity(HIGH);
  p->seek = params.rate_partner[p->sex()][p->activity()][p->age_band()];
  m_women_non.remove(p);
  m_women_csw.insert(p, rate_exit(p));
  p->history.insert(new BehaviorRecord(m_tcurr, p->activity()));
}

void Behavior::csw_exit_int(Person* const p) {
  assert(p->activity() == HIGH);
  assert(m_women_csw.member(p));
  Activity level;
  do level = params.sample_activity_level(m_rng, p->sex()); while(level == HIGH);
  p->activity(level);
  p->seek = params.rate_partner[p->sex()][p->activity()][p->age_band()];
  m_women_csw.remove(p);
  m_women_non.insert(p, rate_init(p));
  p->history.insert(new BehaviorRecord(m_tcurr, p->activity()));
}

void Behavior::debut_ext(Person* const p) {
  if (p->sex() == FEMALE) insert(p);
}

void Behavior::death_ext(Person* const p) {
  if (p->sex() == FEMALE) remove(p);
}

void Behavior::progress_ext(Person* const p) {
  if (p->sex() == FEMALE) update(p);
}

void Behavior::treat_ext(Person* const p) {
  if (p->sex() == FEMALE) update(p);
}

double Behavior::wait_csw_init() {
  // The overall rate that women initiate sex work is calculated as
  // the difference between the target number and actual number of
  // CSWs plus the rate that women leave sex work.
  size_t num_csw, num_non;
  double rate_exit, rate_init;

  num_csw = m_women_csw.size();
  num_non = m_women_non.size();

  rate_exit = m_women_csw.sum();
  rate_init = std::max(0.0, params.prop_activity(FEMALE, HIGH) * num_non - num_csw + rate_exit);

  //  fprintf(stderr, "%s[%d]:%s(%f) %d %d %f %f %f\n", __FILE__, __LINE__, __FUNCTION__,
  //	  m_tcurr, num_csw, num_non, params.prop_activity(FEMALE, HIGH) * num_non, rate_init, rate_exit);

  return gsl_ran_exponential(m_rng, 1.0 / rate_init);
}
