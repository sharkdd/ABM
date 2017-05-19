#include <Treatment.H>

EventType Treatment::output[] = {
  ART_CHANGE, // dropout
  ART_CHANGE, // treatment failure due to non-adherence
  DEATH_HIV,  // HIV mortality on ART
  ART_CHANGE, // transition from early to late ART
  ART_INIT};  // ART initiation

Treatment::Treatment()
  : m_rng(NULL), m_tcurr(0.0), m_pnext(NULL) {

  // populate the event queue with a sentinel record
  m_queue.push(NULL, Record(infinity, ART_EVENTS));

  // initialize the successor function mapping events to the resulting
  // ART state. An invalid successor is specified for MORTALITY for
  // debugging purposes because the previous state should be retained.
  successor[DROPOUT] = ART_NONE;
  successor[FAILURE] = ART_FAIL;
  successor[ADVANCE] = ART_LATE;
  successor[INITIATE] = ART_EARLY;
  successor[MORTALITY] = ART_STATES;

  m_tnext = tnext();
}

Treatment::~Treatment() {
  if (m_rng) gsl_rng_free(m_rng);
}

void Treatment::delta_int() {
  m_tcurr = m_tnext;

  Event e;
  double t;

  // remove the affected person from the current event queue
  m_queue.pop();

  if (m_enext != MORTALITY) {
    // Update the individual's ART status (does not change upon death)
    m_pnext->treat(successor[m_enext]);

    // Check that we never have persons on both ART and PrEP
    assert(m_pnext->prep() == PREP_NONE || m_pnext->treat() == ART_NONE);

    // Special case: patients who restart ART after virologic failure
    // due to acquired resistance re-enter late treatment
    if (m_enext == INITIATE && art_failure(m_pnext)) {
      m_pnext->treat(ART_LATE);
    }

    // schedule the affected person's next event
    next_event(m_pnext, e, t);
    m_queue.push(m_pnext, Record(t, e));
    m_pnext->history.insert(new TreatmentRecord(m_tcurr, m_pnext->treat()));
  } else {
    m_pnext->alive(false);
    m_pnext->ref_dec();
    m_pnext->history.insert(new HIVDeathRecord(m_tcurr));
  }

  // determine the module's next event 
  m_tnext = tnext();
}

void Treatment::delta_ext(double dt, adevs::Bag<Message> const& msgs) {
  m_tcurr += dt;
  adevs::Bag<Message>::const_iterator mi;
  for (mi = msgs.begin(); mi != msgs.end(); ++mi) {
    switch ((*mi).event_type) {
    case DEATH_NAT: death_ext((*mi).person1); break;
    case DEATH_HIV: death_ext((*mi).person1); break;
    case PROGRESS: progress_ext((*mi).person1); break;
    case RESISTANCE: resist_ext((*mi).person1); break;
    case PREP_INIT: prep_ext((*mi).person1); break;
    case PREP_CHANGE: prep_ext((*mi).person1); break;
    default: break;
    }
  }
  m_tnext = tnext();
}

void Treatment::delta_conf(adevs::Bag<Message> const& msgs) {
  delta_ext(0.0, msgs);
  delta_int();
}

void Treatment::output_func(adevs::Bag<Message>& msgs) {
  m_pnext = pnext();
  msgs.insert(Message(output[m_enext], m_pnext));
}

double Treatment::tnext() {
  m_tnext = m_queue.top().second.time;
  m_enext = m_queue.top().second.event;
  return m_tnext;
}

void Treatment::insert(Person* p) {
  if (p->detectable() && (p->prep() == PREP_NONE)) {
    double t;
    Event e;
    next_event(p, e, t);
    m_queue.push(p, Record(t, e));
    p->ref_inc();
  }
}

void Treatment::remove(Person* p) {
  // Only persons with detectable HIV infection are eligible for
  // treatment-related events.
  if (p->detectable() && m_queue.erase(p)) p->ref_dec();
}

void Treatment::progress_ext(Person* p) {
  // remove and reinsert avoids double-counting when maintaining
  // Person reference counts.
  remove(p);
  insert(p);
}

void Treatment::resist_ext(Person* p) {
  if (p->treat() == ART_EARLY && art_failure(p)) {
    // Patients on ART who have virus with acquired ART resistance are
    // categorized as late treatment
    m_queue.push(p, Record(m_tcurr, ADVANCE));
  } else {
    // Resistance acquisition changes the rate of treatment-related
    // events. Calling insert updates that rate.
    //
    // remove and reinsert avoids double-counting when maintaining
    // Person reference counts.
    remove(p);
    insert(p);
  }
}

void Treatment::prep_ext(Person* p) {
  // PrEP users with breakthrough infection are not tracked by the
  // Treatment module until they stop PrEP
  if (p->detectable() && (p->prep() == PREP_NONE)) {
    m_queue.push(p, Record(time_initiate(p), INITIATE));
    p->ref_inc();
  }
}

double Treatment::time_initiate(Person* p) {
  // Infected individuals on PrEP are assumed unaware of their
  // infected status, and therefore do not initiate treatment
  if (p->prep() != PREP_NONE) return infinity;

  // Treatment uptake occurs according to a non-homogeneous Poisson
  // process. This is simulated using the thinning algorithm from pg
  // 78 of Ross, Simulation 3rd Ed. This works by rejection sampling,
  // so we need the first year of ART availability to avoid sampling
  // times when treatment is unavailable.
  double t(std::max(m_tcurr, params.art.year_available[p->stage()]));
  double U;
  while (t < params.year_done) {
    t += gsl_ran_exponential(m_rng, 1 / params.art.rate_uptake_max[p->stage()]);
    U = gsl_rng_uniform(m_rng) * params.art.rate_uptake_max[p->stage()];
    if (U <= params.art.uptake(t,p->stage())) return t;
  }
  return infinity;
}

void Treatment::next_event(Person* p, Event &e, double &t) {
  if (p->treat() == ART_NONE) {
    e = INITIATE;
    t = time_initiate(p);
  } else {
    double rate[ART_EVENTS], sum(0.0);
    rate[INITIATE] = 0.0;

    sum += rate[ADVANCE] = (p->treat() == ART_EARLY) * params.art.rate_advance;
    sum += rate[DROPOUT] = params.art.rate_dropout[p->treat()];

    // virologic failure due to non-adherence is ignored in case of
    // acquired resistance
    sum += rate[FAILURE] = (p->virus() != R2 && p->virus() != C2) * params.art.rate_nonadhere[p->treat()];
    sum += rate[MORTALITY] = rate_mortality(p);

    // next event time
    t = m_tcurr + gsl_ran_exponential(m_rng, 1 / sum);

    // select the next event using Gillespie's direct method
    unsigned int i;
    const double U(gsl_rng_uniform(m_rng) * sum);
    for (i = 0; i < ART_EVENTS && U > rate[i]; ++i, rate[i] += rate[i-1]);
    e = static_cast<Event>(i);    
  }

}
