#include <Prophylaxis.H>

EventType Prophylaxis::output[] = {
  PREP_CHANGE, // PrEP program completion
  PREP_CHANGE, // PrEP injection
  PREP_INIT    // PrEP initiation
};

Prophylaxis::Prophylaxis()
  : m_rng(NULL), m_tcurr(0.0), m_pnext(NULL) {
  for (unsigned int ei(0); ei < PREP_EVENTS; ++ei) {
    m_queue[ei].push(NULL, infinity);
  }

  m_successor[PREP_EARLY] = PREP_LATE;
  m_successor[PREP_LATE ] = PREP_LATE;

  m_tnext = tnext();
}

Prophylaxis::~Prophylaxis() {
  if (m_rng) gsl_rng_free(m_rng);
}

void Prophylaxis::delta_int() {
  m_tcurr = m_tnext;
  m_queue[m_enext].pop();
  switch (m_enext) {
  case FINISH: finish_int(m_pnext); break;
  case INJECT: inject_int(m_pnext); break;
  case INITIATE: initiate_int(m_pnext); break;
  default: assert(false); break;
  }
  m_tnext = tnext();
}

void Prophylaxis::delta_ext(double dt, adevs::Bag<Message> const& msgs) {
  m_tcurr += dt;
  adevs::Bag<Message>::const_iterator mi;
  for (mi = msgs.begin(); mi != msgs.end(); ++mi) {
    switch ((*mi).event_type) {
    case DEBUT: debut_ext((*mi).person1); break;
    case DEATH_NAT: death_ext((*mi).person1); break;
    case DEATH_HIV: death_ext((*mi).person1); break;
    case PROGRESS: progress_ext((*mi).person1); break;
    default: break;
    }
  }
  m_tnext = tnext();  
}

void Prophylaxis::delta_conf(adevs::Bag<Message> const& msgs) {
  delta_ext(0.0, msgs);
  delta_int();
}

void Prophylaxis::output_func(adevs::Bag<Message>& msgs) {
  m_pnext = pnext();
  msgs.insert(Message(output[m_enext], m_pnext));
}

double Prophylaxis::tnext() {
  m_enext = PREP_EVENTS;
  m_tnext = infinity;
  for (unsigned int ei(0); ei < PREP_EVENTS; ++ei) {
    if (m_queue[ei].top().second < m_tnext) {
      m_enext = static_cast<Event>(ei);
      m_tnext = m_queue[ei].top().second;
    }
  }
  return m_tnext;
}

void Prophylaxis::insert(Person* p) {
  if (!p->detectable()) {
    if (p->prep() == PREP_NONE) {
      m_queue[INITIATE].push(p, time_initiate(p));
    } else {
      m_queue[FINISH].push(p, m_tcurr + wait_finish(p));
      m_queue[INJECT].push(p, m_tcurr + wait_inject(p));
    }
    p->ref_inc();
  } else {
    if (p->prep() != PREP_NONE) {
      m_queue[INJECT].push(p, m_tcurr + wait_inject(p));
      m_queue[FINISH].push(p, m_tcurr + wait_finish(p));
      p->ref_inc();
    }
  }
}

void Prophylaxis::remove(Person* p) {
  unsigned int refs(0);
  for (unsigned int ei(0); ei < PREP_EVENTS; ++ei) refs += m_queue[ei].erase(p);
  if (refs) {
    p->ref_dec();
  }
}

void Prophylaxis::finish_int(Person* p) {
  p->prep(PREP_NONE);
  if (!p->detectable()) {
    m_queue[INITIATE].push(p, time_initiate(p));
  } else {
    p->ref_dec();
  }
  m_queue[INJECT].erase(p);
  p->history.insert(new ProphylaxisRecord(m_tcurr, p->prep()));
}

void Prophylaxis::inject_int(Person* p) {
  if (gsl_ran_bernoulli(m_rng, params.prep.prop_missed[p->prep()])) {
    // LTFU when an injection is missed
    p->prep(PREP_NONE);
    m_queue[FINISH].erase(p);
    if (!p->detectable()) {
      m_queue[INITIATE].push(p, time_initiate(p));
    } else {
      p->ref_dec();
    }
  } else {
    // Appointment met. The individual is removed from PrEP if HIV
    // infection is detected, otherwise the person's dose status is
    // updated.
    if (p->detectable()) {
      p->prep(PREP_NONE);
      m_queue[FINISH].erase(p);
      p->ref_dec();
    } else {
      p->prep(m_successor[p->prep()]);
      m_queue[INJECT].push(p, m_tcurr + wait_inject(p));
    }
  }
  p->history.insert(new ProphylaxisRecord(m_tcurr, p->prep()));
}

void Prophylaxis::initiate_int(Person* p) {
  assert(!p->detectable());
  p->prep(PREP_EARLY);
  m_queue[FINISH].push(p, m_tcurr + wait_finish(p));
  m_queue[INJECT].push(p, m_tcurr + wait_inject(p));
  p->history.insert(new ProphylaxisRecord(m_tcurr, p->prep()));
}

void Prophylaxis::debut_ext(Person* p) {
  double t;
  if (gsl_ran_bernoulli(m_rng, params.prep.prop_debut(m_tcurr, p->sex(), p->activity()))) {
    t = m_tcurr;
  } else {
    t = time_initiate(p);
  }
  m_queue[INITIATE].push(p,t);
  p->ref_inc();
}

void Prophylaxis::progress_ext(Person* p) {
  // Individuals who progress to detectable infection can not be
  // enrolled in the PrEP program. Individuals who progress to
  // detectable infection will be removed from PrEP the next time they
  // show up for an injection.
  if (p->detectable() && (p->prep() == PREP_NONE)) {
    if (m_queue[INITIATE].erase(p)) {
      p->ref_dec();
    }
  }
}

double Prophylaxis::time_initiate(Person* p) {
  // PrEP uptake occurs according to a non-homogeneous Poisson
  // process. This is simulated using the thinning algorithm from pg
  // 78 of Ross, Simulation 3rd Ed. This works by rejection sampling,
  // so we need the first year PrEP is available to avoid sampling
  // absurd times.
  double t(std::max(m_tcurr, params.prep.year_available));
  double U;
  while (t < params.year_done) {
    t += gsl_ran_exponential(m_rng, 1 / params.prep.rate_uptake_max);
    U = gsl_rng_uniform(m_rng) * params.prep.rate_uptake_max;
    if (U < (*(params.prep.rate_uptake))(t, p->sex(), p->activity(), p->age_band())) return t;
  }
  return infinity;
}

double Prophylaxis::wait_finish(Person* p) {
  return gsl_ran_exponential(m_rng, params.prep.wait_finish);
}

double Prophylaxis::wait_inject(Person* p) {
  return gsl_ran_exponential(m_rng, params.prep.wait_inject);
}
