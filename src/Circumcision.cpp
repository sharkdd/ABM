#include <Circumcision.H>

Circumcision::Circumcision()
  : m_rng(NULL), m_tcurr(0.0), m_pnext(NULL), m_men_mmc(0) {
  m_tnext = infinity;
}

Circumcision::~Circumcision() {
  if (m_rng) gsl_rng_free(m_rng);
}

void Circumcision::delta_int() {
  m_tcurr = m_tnext;
  circumcise_int(m_pnext);
  m_tnext = tnext();
}

void Circumcision::delta_ext(double dt, adevs::Bag<Message> const& msgs) {
  m_tcurr += dt;
  adevs::Bag<Message>::const_iterator mi;
  for (mi = msgs.begin(); mi != msgs.end(); ++mi) {
    switch ((*mi).event_type) {
    case DEBUT: debut_ext((*mi).person1); break;
    case DEATH_NAT: death_ext((*mi).person1); break;
    case DEATH_HIV: death_ext((*mi).person1); break;
    default: break;
    }
  }
  m_tnext = tnext();  
}

void Circumcision::delta_conf(adevs::Bag<Message> const& msgs) {
  // Handle external events first so that 1) new debut's are eligible
  // and 2) recently deceased individuals are ineligible for MMC
  delta_ext(0.0, msgs);
  delta_int();
}

void Circumcision::output_func(adevs::Bag<Message>& msgs) {
  m_pnext = pnext();
  msgs.insert(Message(MMC, m_pnext));
}

double Circumcision::tnext() {
  // #warning "MMC uptake ignores temporal variation in the coverage target"
  //   const double target(params.mmc.cover_target(m_tcurr));
  //   const double change(params.mmc.cover_change(m_tcurr));
  //   const size_t nu(m_sampler.size());
  //   const size_t nc(m_men_mmc);
  //   const double rate(std::max(0.0, target * (nu + nc) - nc + change * (nu + nc)));
  //   return m_tcurr + gsl_ran_exponential(m_rng, 1.0 / rate);
  return time_circumcise();
}

Person* Circumcision::pnext() const {
  unsigned long int k;
  k = gsl_rng_uniform_int(m_rng, m_sampler.size());
  return m_sampler[k]->first;
}

void Circumcision::insert(Person* p) {
  assert(p->sex() == MALE);
  if (p->circumcised()) {
    ++m_men_mmc;
  } else {
    // The person will be placed at the end of m_sampler, so position
    // will be equal to the size of m_sampler before insertion. We 1)
    // place the (person, position) pair in the hash table (m_lookup),
    // then 2) insert a pointer to the hash table value in m_sampler.
    Lookup::value_type value(p, m_sampler.size());
    std::pair<Lookup::iterator, bool> retval;
    retval = m_lookup.insert(value);
    assert(retval.second);
    m_sampler.push_back(&(*retval.first));
    p->ref_inc();
  }
}

void Circumcision::remove(Person* p) {
  if (p->circumcised()) {
    --m_men_mmc;
  } else {
    size_t j, k;

    Lookup::iterator it(m_lookup.find(p));
    assert(it != m_lookup.end());

    j = it->second;
    k = m_sampler.size() - 1;

    // Replace the record for p in m_sampler with the last record in
    // m_sampler. Link the corresponding lookup table record with that
    // updated sampler record. Then shrink the m_sampler to eliminate
    // the now-obsolete last record.
    m_sampler[j] = m_sampler[k];
    m_sampler[j]->second = j;
    m_sampler.resize(m_sampler.size() - 1);

    // Erase p's record from the lookup table and decrease p's
    // reference count.
    m_lookup.erase(it);
    p->ref_dec();
  }
}

void Circumcision::circumcise_int(Person* p) {
  // Circumcising a man removes him from the module, since to
  // calculate MMC uptake rates, we need to know how many men are
  // circumcised, but we do not need to know who they are.
  remove(p);
  p->history.insert(new CircumcisionRecord(m_tcurr));
  p->circumcised(true);
  ++m_men_mmc;
}

double Circumcision::time_circumcise() {
  // This may be well-approximated by a homogeneous Poisson process
  // with rate
  //
  //   target = coverage(m_tcurr).first;
  //   change = coverage(m_tcurr).second;
  //   rate = std::max(0.0, target * (nu + nc) - nc + change * (nu + nc));
  //
  // where nu and nc are as defined below

  const size_t nu(m_sampler.size());
  const size_t nc(m_men_mmc);
  double t(std::max(m_tcurr, params.mmc.year_init)), U, rate, rmax;
  double target, change;
  std::pair<double,double> cvg;

  rmax = std::max(0.0, params.mmc.target_max * (nu + nc) - nc + params.mmc.change_max * (nu + nc));
  while (t < params.mmc.year_done) {
    cvg = params.mmc.coverage(t);
    target = cvg.first;
    change = cvg.second;
    rate = std::max(0.0, target * (nu + nc) - nc + change * (nu + nc));
    t += gsl_ran_exponential(m_rng, 1.0 / rmax);
    U = gsl_rng_uniform(m_rng) * rmax;
    if (U < rate) return t;
  }
  return infinity;
}
