#include <Resistance.H>

Resistance::Resistance() 
  : m_rng(NULL), m_tcurr(0.0), m_pnext(NULL) {
  for (unsigned int ei(0); ei < RESIST_EVENTS; ++ei) {
    m_queue[ei].push(NULL, infinity);
  }
  m_tnext = tnext();

  // record the HIV variants that emerge when majority resistant
  // variants emerge. 
  m_revert[ R1] = WR1;
  m_revert[ C1] = WC1;
  m_revert[ Q1] = WQ1;
  m_revert[ R2] = WR2;
  m_revert[ C2] = WC2;
  m_revert[ Q2] = WQ2;

  // Wild-type and minority resistant variants do not revert, so their
  // successor variants are marked with the invalid value VARIANTS
  m_revert[ WT] = VARIANTS;
  m_revert[WR1] = VARIANTS;
  m_revert[WC1] = VARIANTS;
  m_revert[WQ1] = VARIANTS;
  m_revert[WR2] = VARIANTS;
  m_revert[WC2] = VARIANTS;
  m_revert[WQ2] = VARIANTS;
}

Resistance::~Resistance() {
  if (m_rng) gsl_rng_free(m_rng);
}

void Resistance::delta_int() {
  m_tcurr = m_tnext;
  m_queue[m_enext].pop();
  switch (m_enext) {
  case EMERGE: emerge_int(m_pnext); break;
  case REVERT: revert_int(m_pnext); break;
  default: assert(false); break; // should never get here
  }
  m_tnext = tnext();

  // check if resistance emerges or reverts inappropriately
  assert(m_pnext->virus() != VARIANTS);
}

void Resistance::delta_ext(double dt, adevs::Bag<Message> const& msgs) {
  m_tcurr += dt;
  adevs::Bag<Message>::const_iterator mi;
  for (mi = msgs.begin(); mi != msgs.end(); ++mi) {
    switch ((*mi).event_type) {
    case DEATH_NAT: death_ext((*mi).person1); break;
    case DEATH_HIV: death_ext((*mi).person1); break;
    case INFECT: infect_ext((*mi).person2); break;
    case ART_INIT: treat_ext((*mi).person1); break;
    case ART_CHANGE: treat_ext((*mi).person1); break;
    case PREP_INIT: prep_ext((*mi).person1); break;
    case PREP_CHANGE: prep_ext((*mi).person1); break;
    default: break;
    }
  }
  m_tnext = tnext();
}

void Resistance::delta_conf(adevs::Bag<Message> const& msgs) {
  delta_ext(0.0, msgs);
  delta_int();
}

void Resistance::output_func(adevs::Bag<Message>& msgs) {
  m_pnext = pnext();
  msgs.insert(Message(RESISTANCE, m_pnext));
}

double Resistance::tnext() {
  m_enext = static_cast<Event>(0);
  m_tnext = m_queue[m_enext].top().second;
  for (unsigned int ei(1); ei < RESIST_EVENTS; ++ei) {
    if (m_queue[ei].top().second < m_tnext) {
      m_enext = static_cast<Event>(ei);
      m_tnext = m_queue[ei].top().second;
    }
  }
  assert(m_tnext >= m_tcurr);
  return m_tnext;
}

void Resistance::insert(Person* p) {
  // This supports either reversion or emergence of resistance in any
  // ARV state with any virus. For example, this would allow
  // experiments in which resistance may either emerge or revert while
  // non-adherent to ART. It is the experiment designer's
  // responsibility to assign parameter values (params.wait_emerge and
  // params.wait_revert) that implement model specifications, for
  // example, by setting wait_emerge=infinity in situations where
  // resistance can not emerge.
  const double tr(m_tcurr + wait_revert(p));
  const double te(m_tcurr + wait_emerge(p));
  double t;
  Event e;
  if (tr < te) {
    t = tr;
    e = REVERT;
  } else {
    t = te;
    e = EMERGE;
  }
  m_queue[e].push(p,t);
  p->ref_inc();
}

void Resistance::remove(Person* p) {
  unsigned int refs(0);
  for (unsigned int ei(0); ei < RESIST_EVENTS; ++ei) refs += m_queue[ei].erase(p);
  if (refs) p->ref_dec();
}

void Resistance::emerge_int(Person* p) {
  Virus virus;
  if (p->treat() != ART_NONE) {
    if (p->virus() == WT || p->virus() == R1 || p->virus() == WR1) {
      virus = (gsl_ran_bernoulli(m_rng, params.art.prop_cross_resist) ? C2 : R2);
    } else if (p->virus() == WR2) {
      virus = R2;
    } else {
      virus = C2;
    }
  } else if (p->prep() != PREP_NONE) {
    virus = Q2;
  } else {
    // should never get here, since resistance only emerges with ART
    // or PrEP pressure in our model. If this assertion fails,
    // params.wait_emerge may incorrectly be finite without ARVs.
    virus = VARIANTS;
    assert(false);
  }
  p->ref_dec();
  p->virus(virus);
  p->history.insert(new ResistanceRecord(m_tcurr, p->virus()));
}

void Resistance::revert_int(Person* p) {
  p->virus(m_revert[p->virus()]);
  p->history.insert(new ResistanceRecord(m_tcurr, p->virus()));

  // This module allows reversion of resistant variants under drug
  // pressure; for example, PrEP does not maintain ART resistance if
  // that variant lacks PrEP cross resistance. In this situation, a
  // person on PrEP is at risk of both wild-type reversion and PrEP
  // resistance acquisition. If reversion to wild-type occurs, he or
  // she remains at risk for acquiring PrEP resistance
  const double t(m_tcurr + wait_emerge(p));
  if (t < infinity) {
    m_queue[EMERGE].push(p,t);
  } else {
    p->ref_dec();
  }
}

void Resistance::treat_ext(Person* p) {
  remove(p);
  insert(p);
}

void Resistance::prep_ext(Person* p) {
  if (p->susceptible()) return; // Ignore susceptibles
  remove(p);
  insert(p);
}

double Resistance::wait_revert(Person* p) {
  return gsl_ran_exponential(m_rng, params.wait_revert[p->virus()][p->treat()][p->prep()]);
}

double Resistance::wait_emerge(Person* p) {
  return gsl_ran_exponential(m_rng, params.wait_emerge[p->virus()][p->treat()][p->prep()]);
}
