#include <algorithm>
#include <Partnership.H>

// The GEOMETRIC compiler flag selects between two different ways of
// calculating the partnership formation rate. When GEOMETRIC is not
// set, pairs form partnerships based on the arithmetic mean of their
// seek rates, as in the original GERMS model. Pairs form partnerships
// based on the geometric mean when GEOMETRIC is set, which closely
// approximates garnett1994mmb-type mixing and balancing.

Partnership::Partnership()
  : m_rng(NULL), m_tcurr(0.0), m_next(EVENTS_INT) {
  for (size_t ei(0); ei <= EVENTS_INT; ++ei) m_time[ei] = infinity;
  // insert a sentinel record to avoid handling an empty breakup queue
  m_breakup.push(Match(NULL, NULL), infinity);
}

Partnership::~Partnership() {
  if (m_rng) gsl_rng_free(m_rng);
}

void Partnership::insert(Person* const person) {
  const double seek_base(Partnership::seek(person));
  int index[6], n;
  double weight[6];
  bins(person, index, weight, n);
  for (int b(0); b < n; ++b) {
    m_bin[index[b]].insert(person, weight[b] * seek_base);
  }
  // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, person, person->ref_num(), m_tcurr);
  person->ref_inc();
}

void Partnership::remove(Person* const person) {
  int index[6], n;
  double weight[6];
  bins(person, index, weight, n);
  for (int b(0); b < n; ++b) m_bin[index[b]].remove(person);
  // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, person, person->ref_num(), m_tcurr);
  person->ref_dec();
}

void Partnership::update(Person* const person) {
  const double seek_base(Partnership::seek(person));
  int index[6], n;
  double weight[6];
  bins(person, index, weight, n);
  for (int b(0); b < n; ++b) {
    m_bin[index[b]].update(person, weight[b] * seek_base);
  }
}

double Partnership::seek(Person const* const person) {
  size_t ni(person->partners.size());
  double ri(person->seek * params.prop_seek_decr[person->stage()][person->treat()]);
  double di(person->damp);
  return ri * gsl_pow_int(di, ni);
}

void Partnership::bins(Person const* const person, int *index, double *weight, int &n) {
  const AgeBand a(static_cast<AgeBand>(person->age_band()));
  const Activity k(static_cast<Activity>(person->activity()));

  // all sexually active individuals are guaranteed to seek partners
  // in the random mixing bin and one each of same-age, same-level and
  // same-age + same-level preferential mixing bins.
  index[0] = BIN_RAND;
  index[1] = BIN_ACT + k;
  index[2] = BIN_SAGE + a;
  index[3] = BIN_SAGE_ACT + a * LEVELS + k;

  // all but the youngest men and oldest women may preferentially seek
  // disparate-age partnerships
  if (person->sex() == MALE) {
    if (a > 0) {
      index[4] = BIN_DAGE + a - 1;
      index[5] = BIN_DAGE_ACT + (a - 1) * LEVELS + k;
      n = 6;
    } else {
      n = 4;
    }
  } else {
    if (a + 1 < AGES) {
      index[4] = BIN_DAGE + a;
      index[5] = BIN_DAGE_ACT + a * LEVELS + k;
      n = 6;
    } else {
      n = 4;
    }
  }

  // index[0] random mixing
  // index[1] risk-assortative, age-random; 
  // index[2] risk-random, same-age-assortative
  // index[3] risk-assortative, same-age-assortative
  // index[4] risk-random, disparate-age-assortative
  // index[5] risk-assortative, disparate-age-assortative
  weight[0] = (1.0 - params.assort_act) * (1.0 - params.assort_age); 
  weight[1] = params.assort_act * (1.0 - params.assort_age); 
  if (n == 4) {
    weight[2] = (1.0 - params.assort_act) * params.assort_age;
    weight[3] = params.assort_act * params.assort_age;
  } else if (n == 6) {
    weight[2] = (1.0 - params.assort_act) * params.assort_age * (1.0 - params.assort_dif); 
    weight[3] = params.assort_act * params.assort_age * (1.0 - params.assort_dif); 
    weight[4] = (1.0 - params.assort_act) * params.assort_age * params.assort_dif;
    weight[5] = params.assort_act * params.assort_age * params.assort_dif;
  } else {
    assert(false);
  }

#warning "Debugging code in bins calculation"
  // DEBUGGING CODE
  double sum(std::accumulate(weight, weight + n, 0.0));
  assert(fabs(sum - 1.0) < 1e-12);
  // END DEBUGGING CODE
}

void Partnership::delta_int() {
  m_tcurr = m_time[m_next];
  switch (m_next) {
  case PARTNER_INT: partner_int(); break;
  case BREAKUP_INT: breakup_int(); break;
  case DEATH_INT: death_int(); break;
  case DEBUT_INT: debut_int(); break;
  case AGE_INT: age_int(); break;
  case RISK_INT: risk_int(); break;
  case SEEK_INT: seek_int(); break;
  default: break;
  }
  m_time[PARTNER_INT] = tnext();
  m_time[BREAKUP_INT] = m_breakup.top().second;
  m_next = EVENTS_INT;

  for (size_t ei(0); ei < EVENTS_INT; ++ei)
    if (m_time[ei] < m_time[m_next])
      m_next = static_cast<Event>(ei);
}

void Partnership::delta_ext(double dt, Bag<Message> const& msgs) {
  // Dispatches handlers for external events. The handlers schedule
  // immediate internal transitions for these events, so it is not
  // necessary to find the next event time here
  m_tcurr += dt;
  Bag<Message>::const_iterator mi;
  for (mi = msgs.begin(); mi != msgs.end(); ++mi) {
    switch ((*mi).event_type) {
    case AGE: schedule_ext((*mi).person1, AGE_INT); break;
    case DEBUT: schedule_ext((*mi).person1, DEBUT_INT); break;
    case DEATH_NAT: schedule_ext((*mi).person1, DEATH_INT); break;
    case DEATH_HIV: schedule_ext((*mi).person1, DEATH_INT); break;
    case RISK_CHANGE: schedule_ext((*mi).person1, RISK_INT); break;
    case PROGRESS: schedule_ext((*mi).person1, SEEK_INT); break;
    case ART_INIT: schedule_ext((*mi).person1, SEEK_INT); break;
    case ART_CHANGE: schedule_ext((*mi).person1, SEEK_INT); break;
    default: break;
    }
  }
}

void Partnership::delta_conf(Bag<Message> const& msgs) {
  delta_int();
  delta_ext(0.0, msgs);
}

void Partnership::output_func(Bag<Message>& msgs) {
  switch (m_next) {
  case PARTNER_INT: partner_out(msgs); break;
  case BREAKUP_INT: breakup_out(msgs); break;
  case DEATH_INT: death_out(msgs); break;
  default: break;
  }
}

double Partnership::rate() {
  double rate(0.0);
  for (int k(0); k < BINS; ++k) {
    m_rate[k] = m_bin[k].rate();
    rate += m_rate[k];
  }
  return rate;
}

void Partnership::partner_int() {
  //   // BGN DEBUG: capture first partnership formation event
  // #warning "Partnership module terminates simulation early for debugging purposes"
  //   assert(m_mnext.first && m_mnext.second);
  //   fprintf(stderr, "%f %d %d %d %d\n", m_tcurr,
  // 	  m_mnext.second->age_band(), // female partner
  // 	  m_mnext.second->activity(),
  // 	  m_mnext.first->age_band(), // male partner
  // 	  m_mnext.first->activity());
  //   exit(1);
  //   // END DEBUG: capture first partnership formation event

  if (m_mnext.first && m_mnext.second) {
    const Activity ki(m_mnext.first->activity());
    const Activity kj(m_mnext.second->activity());
    const double mean(1.0 / params.rate_breakup[ki][kj]);

    m_breakup.push(m_mnext, m_tcurr + gsl_ran_exponential(m_rng, mean));

#ifdef PARTNER_HISTORY
    // When enabled, all partnership formation events are stored in
    // each partners' history. Each person may experience a large
    // number of these events and most are irrelevant for HIV
    // transmission, so we want these disabled in large-scale
    // simulations
    m_mnext.first->history.insert(new PartnerRecord(m_tcurr, m_mnext.second));
    m_mnext.second->history.insert(new PartnerRecord(m_tcurr, m_mnext.first));
#endif // PARTNER_HISTORY

    // Update each persons' partner set
    m_mnext.first->partners.insert(m_mnext.second);
    m_mnext.second->partners.insert(m_mnext.first);

    // Update each persons' partnership formation rate
    update(m_mnext.first);
    update(m_mnext.second);
  }
}

void Partnership::breakup_int() {
  const Match m(m_breakup.top().first);
  m_breakup.pop();

#ifdef PARTNER_HISTORY
  // When enabled, all partnership dissolution events are stored in
  // each partners' history. Each person may experience a large
  // number of these events and most are irrelevant, so we do not
  // want these enabled in production
  m.first->history.insert(new BreakupRecord(m_tcurr, m.second));
  m.second->history.insert(new BreakupRecord(m_tcurr, m.first));
#endif // PARTNER_HISTORY

  // Update each persons' partner set
  m.first->partners.erase(m.second);
  m.second->partners.erase(m.first);

  // Update each persons' partnership formation rate
  update(m.first);
  update(m.second);
}

void Partnership::debut_int() {
  EventCache::const_iterator di;
  for (di = m_cache[DEBUT_INT].begin(); di != m_cache[DEBUT_INT].end(); ++di) {
    insert(*di);
  }
  m_cache[DEBUT_INT].clear();
  m_time[DEBUT_INT] = infinity;
}

void Partnership::death_int() {
  // Dissolve partnerships with each deceased person.
  EventCache::const_iterator di;
  Person::Partners::iterator pi;
  for (di = m_cache[DEATH_INT].begin(); di != m_cache[DEATH_INT].end(); ++di) {
    for (pi = (*di)->partners.begin(); pi != (*di)->partners.end(); ++pi) {
#ifdef PARTNER_HISTORY
      // Store breakup events indicating mortality-related dissolution
      (*di)->history.insert(new BreakupRecord(m_tcurr, *pi));
      (*pi)->history.insert(new BreakupRecord(m_tcurr, *di));
#endif // PARTNER_HISTORY
      (*pi)->partners.erase(*di);
      update(*pi);
    }

    // Deschedule breakup records for individuals who have died
    if ((*di)->sex() == MALE) {
      for (pi = (*di)->partners.begin(); pi != (*di)->partners.end(); ++pi) {
	m_breakup.erase(Match(*di, *pi));
      }
    } else {
      for (pi = (*di)->partners.begin(); pi != (*di)->partners.end(); ++pi) {
	m_breakup.erase(Match(*pi, *di));
      }
    }

    (*di)->partners.clear();
    remove(*di);
  }
  m_cache[DEATH_INT].clear();
  m_time[DEATH_INT] = infinity;
}

void Partnership::age_int() {
  Person* person;
  EventCache::const_iterator pi;
  int index[6], n;
  double weight[6], seek_base;

  for (pi = m_cache[AGE_INT].begin(); pi != m_cache[AGE_INT].end(); ++pi) {
    person = *pi;
    seek_base = Partnership::seek(person);

    // This implementation is currently overkill, but avoids problems
    // in the corner case of an individual who ages and changes sexual
    // activity levels simultaneously.

    // 1. Remove from all bins (we could exploit that these people
    //    should be present only in 4 or 6 bins).
    for (int b(0); b < BINS; ++b) {
      if (m_bin[b].member(person)) m_bin[b].remove(person);
    }

    // 2. Insert in appropriate bins
    bins(person, index, weight, n);
    for (int b(0); b < n; ++b) m_bin[index[b]].insert(person, weight[b] * seek_base);
  }
  m_cache[AGE_INT].clear();
  m_time[AGE_INT] = infinity;
}

void Partnership::risk_int() {
  Person* person;
  EventCache::const_iterator pi;
  int index[6], n;
  double weight[6], seek_base;
  for (pi = m_cache[RISK_INT].begin(); pi != m_cache[RISK_INT].end(); ++pi) {
    person = *pi;
    seek_base = Partnership::seek(person);

    // 1. Remove person from previous risk-assortative mixing bins.
    //    The partnership module does not know the person's previous
    //    sexual activity level, so we check each bin to be safe. (We
    //    could exploit that these people should be present in only 4
    //    or 6 bins).
    for (int b(0); b < BINS; ++b) {
      if (m_bin[b].member(person)) m_bin[b].remove(person);
    }

    // 2. Reinsert person
    bins(person, index, weight, n);
    for (int b(0); b < n; ++b) m_bin[index[b]].insert(person, weight[b] * seek_base);
  }
  m_cache[RISK_INT].clear();
  m_time[RISK_INT] = infinity;
}

void Partnership::seek_int() {
  EventCache::const_iterator pi;
  for (pi = m_cache[SEEK_INT].begin(); pi != m_cache[SEEK_INT].end(); ++pi) {
    update(*pi);
  }
  m_cache[SEEK_INT].clear();
  m_time[SEEK_INT] = infinity;
}

void Partnership::partner_out(Bag<Message>& msgs) {
  m_mnext = mnext();

  // See the comments on tnext() and mnext(). Partnership seek rate
  // computation has been optimized by pretending that all
  // partnerships are eligible to form, including existing
  // partnerships. A candidate partnership (m_mnext) is selected at
  // the next formation time, but if its members are already partners,
  // the match is invalid and the member pointers are nullified.
  //
  // To reduce unnecessary messaging, output is generated only when
  // serodiscordant partnerships are formed
  if (m_mnext.first && m_mnext.second && (m_mnext.first->infected() ^ m_mnext.second->infected())) {
    msgs.insert(Message(PARTNER, m_mnext.first, m_mnext.second));
  }
}

void Partnership::breakup_out(Bag<Message>& msgs) {
  Match const& m(m_breakup.top().first);
  // To reduce unnecessary messaging, output is generated only when a
  // serodiscordant partnership has ended
  if (m.first->infected() ^ m.second->infected()) {
    msgs.insert(Message(BREAKUP, m.first, m.second));
  }
}

void Partnership::death_out(Bag<Message>& msgs) {
  EventCache::const_iterator di;
  Person::Partners::iterator pi;
  for (di = m_cache[DEATH_INT].begin(); di != m_cache[DEATH_INT].end(); ++di) {
    for (pi = (*di)->partners.begin(); pi != (*di)->partners.end(); ++pi) {
      // To reduce unnecessary messaging, output is generated only if
      // the partnership that ended was serodiscordant
      if ((*di)->infected() ^ (*pi)->infected()) {
	msgs.insert(Message(BREAKUP, *di, *pi));
      }
    }
  }
}

void Partnership::schedule_ext(Person* p, Event e) {
  m_cache[e].insert(m_cache[e].begin(), p);
  m_time[e] = m_tcurr;
  m_next = e;
}

double Partnership::tnext() {
  // The overall rate of partnership formation is the sum of rates in
  // each bin. Bin rates are weighted by the degree of assortativity
  // to account for the proportion of time persons spend seeking
  // partnerships proportionally or assortatively

  double rate(0.0);
  for (int k(0); k < BINS; ++k) {
    m_rate[k] = m_bin[k].rate();
    rate += m_rate[k];
  }

  return m_tcurr + gsl_ran_exponential(m_rng, 1 / rate);
}

Partnership::Match Partnership::mnext() const {
  // partnerships come from bins proportional to the contribution of
  // each bin to the overall rate
  double psum[BINS];
  psum[0] = m_rate[0];
  for (int ai(1); ai < BINS; ++ai) psum[ai] = m_rate[ai] + psum[ai-1];
  const double U(psum[BINS - 1] * gsl_rng_uniform(m_rng));
  for (int ai(0); ai < BINS; ++ai)
    if (U < psum[ai]) return m_bin[ai].sample(m_rng);
  assert(false); // should never get here
}

Partnership::Bin::Bin() {}
Partnership::Bin::~Bin() {}

void Partnership::Bin::insert(Person* const person, double const seek_rate) {
  m_valid = false;
  seek[person->sex()].insert(person, seek_rate);
}

void Partnership::Bin::remove(Person* const person) {
  m_valid = false;
  seek[person->sex()].remove(person);
}

void Partnership::Bin::update(Person* const person, double const seek_rate) {
  m_valid = false;
  seek[person->sex()].update(person, seek_rate);
}

bool Partnership::Bin::member(Person* const person) {
  return seek[person->sex()].member(person);
}

double Partnership::Bin::rate() {
  if (!m_valid) {
    // New partnerships form at rate sqrt(Sm * Sf) - sqrt(Sp) in a bin,
    // where Sm and Sf are the sums of male and female partnership seek
    // rates and Sp is the contribution to Sm * Sf from members of
    // existing partnerships. This has been optimized by ignoring the Sp
    // term and allowing selection of invalid pairings at the when a
    // partnership formation event occurs. These invalid pairings are
    // ignored. This is provably equivalent to sampling partnerships at
    // the rate above, but is more space-efficient
    const double Sm(seek[ MALE ].sum());
    const double Sf(seek[FEMALE].sum());
#ifdef GEOMETRIC
    m_rate = sqrt(Sm * Sf);
#else
    m_rate = 0.5 * (Sm + Sf);
#endif
    m_valid = true;
  }
  return m_rate;
}

Partnership::Match Partnership::Bin::sample(gsl_rng* rng) const {
  // tnext() was optimized by assuming that the number of existing
  // partnerships is much smaller than the number of potential
  // pairings. We took advantage of this assumption by summing
  // partnership formation rates among all possible pairings without
  // subtracting the sum of formation rates among existing
  // partnerships. For the simulation to remain correct, we have to
  // select the existing partnership among all possible pairings; if
  // the selected pair are already partners, we mark that match as
  // invalid and continue on.
  //
  // Notice that rejection sampling (resampling until we get a valid
  // pair) results in an incorrect algorithm, since new partnerships
  // would form too rapidly.
  Match m;
  m.first  = seek[ MALE ].sample(rng);
  m.second = seek[FEMALE].sample(rng);
  if (m.first->partners.find(m.second) != m.first->partners.end()) {
    m.first = m.second = NULL;
  }
  return m;
}
