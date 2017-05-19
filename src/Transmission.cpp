#include <gsl/gsl_randist.h>
#include <Transmission.H>

Transmission::Transmission() : m_rng(NULL), m_tcurr(0.0), m_tnext(infinity) {
  m_virus[WT ] = WT;
  m_virus[R1 ] = R1;
  m_virus[R2 ] = R1;
  m_virus[C1 ] = C1;
  m_virus[C2 ] = C1;
  m_virus[Q1 ] = Q1;
  m_virus[Q2 ] = Q1;
  m_virus[WR1] = WT;
  m_virus[WR2] = WT;
  m_virus[WC1] = WT;
  m_virus[WC2] = WT;
  m_virus[WQ1] = WT;
  m_virus[WQ2] = WT;
}

Transmission::~Transmission() {
  if (m_rng) gsl_rng_free(m_rng);
}

double Transmission::ta() {return m_tnext - m_tcurr;}

void Transmission::delta_int() {
  //  printf("%10.6f %8lu %10.6f\n", m_tcurr, m_transmit.size(), m_tnext - m_tcurr);
  m_tcurr = m_tnext;

  // upon becoming infected, all of the previously negative partner's
  // concordant negative partnerships become serodiscordant, and
  // serodiscordant partnerships become concordant positive
  Person* p(m_cnext.second);
  p->stage(ACUTE_WINDOW);
  p->virus(m_virus[m_cnext.first->virus()]);
  p->history.insert(new InfectRecord(m_tcurr, m_cnext.first, p->virus()));
  Person::Partners::const_iterator qi;
  for (qi = p->partners.begin(); qi != p->partners.end(); ++qi) {
    if ((*qi)->susceptible()) insert(Couple(p, *qi));
    else remove(Couple(*qi, p));
  }
  m_tnext = m_tcurr + tnext();
}

void Transmission::delta_ext(double dt, adevs::Bag<Message> const& msgs) {
  m_tcurr += dt;
  adevs::Bag<Message>::const_iterator mi;
  for (mi = msgs.begin(); mi != msgs.end(); ++mi) {
    switch ((*mi).event_type) {
    case PARTNER: partner_ext((*mi).person1, (*mi).person2); break;
    case BREAKUP: breakup_ext((*mi).person1, (*mi).person2); break;
    case PROGRESS: progress_ext((*mi).person1); break;
    case ART_INIT: treat_ext((*mi).person1); break;
    case ART_CHANGE: treat_ext((*mi).person1); break;
    case PREP_INIT: prep_ext((*mi).person1); break;
    case PREP_CHANGE: prep_ext((*mi).person1); break;
    case MMC: mmc_ext((*mi).person1); break;
    case RESISTANCE: resist_ext((*mi).person1); break;
    default: break;
    }
  }
  m_tnext = m_tcurr + tnext();
}

void Transmission::delta_conf(adevs::Bag<Message> const& msgs) {
  delta_int();
  delta_ext(0.0, msgs);
}

void Transmission::output_func(adevs::Bag<Message>& msgs) {
  m_cnext = cnext();
  msgs.insert(Message(INFECT, m_cnext.first, m_cnext.second));
}

void Transmission::partner_ext(Person* p1, Person* p2) {
  // This assumes the transmission rates in a partnership are
  // independent of each persons' number of partners.
  //
  // If transmission rates depend on numbers of partners, then we need
  // to update transmission rates to susceptible partners of any
  // infected people among p1 and p2. Suppose, for example, that each
  // person has a fixed number of acts per year, but those acts are
  // divided evenly among all partners. Then a new partnership among
  // two infected people will reduce the rate of transmission to any
  // susceptible partners they have.
  if (p1->infected()) {
    if (p2->susceptible()) insert(Couple(p1,p2));
  } else {
    if (p2->infected()) insert(Couple(p2,p1));
  }
}

void Transmission::breakup_ext(Person* p1, Person* p2) {
  // As in partner_ext, if transmission rates depend on the number of
  // partners, then dissolution of concordant positive partnerships
  // requires re-evaluation of transmission risks in other
  // partnerships involving p1 or p2.
  if (p1->infected()) {
    if (p2->susceptible()) remove(Couple(p1,p2));
  } else {
    if (p2->infected()) remove(Couple(p2,p1));
  }
}

void Transmission::prep_ext(Person* p) {
  // The rate of acquisition from each seropositive partner may change
  // after a susceptible individual changes PrEP status
  if (p->susceptible()) {
    Couple couple(NULL, p);
    Person::Partners::const_iterator qi;
    for (qi = p->partners.begin(); qi != p->partners.end(); ++qi) {
      if ((*qi)->infected()) {
	couple.first = *qi;
	Index::iterator i(m_couple2index.find(couple));
	assert(i != m_couple2index.end());
	m_transmit.set(i->second, risk(couple));
      }
    }
  }
}

void Transmission::mmc_ext(Person* p) {
  // The rate of acquisition from each seropositive partner may change
  // after a susceptible individual is circumcised
  if (p->susceptible()) {
    Couple couple(NULL, p);
    Person::Partners::const_iterator qi;
    for (qi = p->partners.begin(); qi != p->partners.end(); ++qi) {
      if ((*qi)->infected()) {
	couple.first = *qi;
	Index::iterator i(m_couple2index.find(couple));
	assert(i != m_couple2index.end());
	m_transmit.set(i->second, risk(couple));
      }
    }
  }
}

void Transmission::update(Person* p) {
  // The rate of transmission to each seronegative partner may change
  // after a change in infectiousness
  Couple couple(p, NULL);
  Person::Partners::const_iterator qi;
  for (qi = p->partners.begin(); qi != p->partners.end(); ++qi) {
    if ((*qi)->susceptible()) {
      couple.second = *qi;
      Index::iterator i(m_couple2index.find(couple));
      assert(i != m_couple2index.end());
      m_transmit.set(i->second, risk(couple));
    }
  }
}

double Transmission::tnext() const {
  const double rate(m_transmit.sum());
  return (rate > 0.0) ? gsl_ran_exponential(m_rng, 1.0 / rate) : infinity;
}

Transmission::Couple Transmission::cnext() const {
  const double U(gsl_rng_uniform(m_rng) * m_transmit.sum());
  return m_index2couple[m_transmit.search(U)];
}

void Transmission::insert(Couple couple) {
  if (m_couple2index.find(couple) == m_couple2index.end()) {
    size_t index(m_index2couple.size());
    m_couple2index[couple] = index;
    m_index2couple.push_back(couple);
    m_transmit.insert(risk(couple));

    // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, couple.first, couple.first->ref_num(), m_tcurr);
    // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, couple.second, couple.second->ref_num(), m_tcurr);
    couple.first->ref_inc();
    couple.second->ref_inc();
  }
}

void Transmission::remove(Couple couple) {
  Index::iterator ci(m_couple2index.find(couple));
  size_t index;
  Couple cswap;
  if (ci != m_couple2index.end()) {
    index = m_transmit.remove(ci->second);
    cswap = m_index2couple[index];
    m_index2couple[ci->second] = cswap;
    m_index2couple.resize(m_index2couple.size()-1);
    m_couple2index[cswap] = ci->second;
    m_couple2index.erase(ci);

    // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, couple.first, couple.first->ref_num(), m_tcurr);
    // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, couple.second, couple.second->ref_num(), m_tcurr);
    couple.first->ref_dec();
    couple.second->ref_dec();
  }
}

double Transmission::risk(Couple const& couple) {
  return params.rate_transmit(couple.first, couple.second, m_tcurr);
}
