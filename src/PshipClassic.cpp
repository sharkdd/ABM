#include <algorithm>
#include <PshipClassic.H>

Partnership::Partnership()
  : m_rng(NULL), m_tcurr(0.0), m_next(EVENTS_INT) {
  for (size_t ei(0); ei <= EVENTS_INT; ++ei) m_time[ei] = infinity;

  // initialize supply caches
  for (int ai(0); ai < AGES; ++ai) {
    std::fill(m_supply[ MALE ][ai], m_supply[ MALE ][ai] + LEVELS, 0.0);
    std::fill(m_supply[FEMALE][ai], m_supply[FEMALE][ai] + LEVELS, 0.0);
  }
  std::fill(m_supply_age[ MALE ], m_supply_age[ MALE ] + AGES, 0.0);
  std::fill(m_supply_age[FEMALE], m_supply_age[FEMALE] + AGES, 0.0);
  m_supply_sex[ MALE ] = m_supply_sex[FEMALE] = 0.0;

  // insert a sentinel record to avoid handling an empty breakup queue
  m_breakup.push(Match(NULL, NULL), infinity);
}

Partnership::~Partnership() {
  if (m_rng) gsl_rng_free(m_rng);
}

void Partnership::insert(Person* const person) {
  const Sex gj(person->sex());
  const int kj(person->activity()), aj(person->age_band());
  m_seek[gj][aj][kj].insert(person, Partnership::seek(person));

  // update supply caches
  m_supply[gj][aj][kj] = m_seek[gj][aj][kj].sum();
  m_supply_age[gj][aj] = std::accumulate(m_supply[gj][aj], m_supply[gj][aj] + LEVELS, 0.0);
  m_supply_sex[gj] = std::accumulate(m_supply_age[gj], m_supply_age[gj] + AGES, 0.0);

  // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, person, person->ref_num(), m_tcurr);
  m_registry[person] = Record(aj, kj);
  person->ref_inc();
}

void Partnership::remove(Person* const person) {
  const Sex gj(person->sex());
  const int kj(person->activity()), aj(person->age_band());
  m_seek[gj][aj][kj].remove(person);

  // update supply caches
  m_supply[gj][aj][kj] = m_seek[gj][aj][kj].sum();
  m_supply_age[gj][aj] = std::accumulate(m_supply[gj][aj], m_supply[gj][aj] + LEVELS, 0.0);
  m_supply_sex[gj] = std::accumulate(m_supply_age[gj], m_supply_age[gj] + AGES, 0.0);

  // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, person, person->ref_num(), m_tcurr);
  m_registry.erase(person);
  person->ref_dec();
}

void Partnership::update(Person* const person) {
  const Sex gj(person->sex());
  const int kj(person->activity()), aj(person->age_band());

  // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, person, person->ref_num(), m_tcurr);

  assert(m_seek[gj][aj][kj].member(person));

  m_seek[gj][aj][kj].update(person, Partnership::seek(person));

  // update supply caches
  m_supply[gj][aj][kj] = m_seek[gj][aj][kj].sum();
  m_supply_age[gj][aj] = std::accumulate(m_supply[gj][aj], m_supply[gj][aj] + LEVELS, 0.0);
  m_supply_sex[gj] = std::accumulate(m_supply_age[gj], m_supply_age[gj] + AGES, 0.0);
}

void Partnership::reassign(Person* const person) {
  Registry::iterator it(m_registry.find(person));
  assert(it != m_registry.end());
  const Sex gj(person->sex());
  const int kcurr(person->activity()), acurr(person->age_band());
  const int kprev(it->second.activity), aprev(it->second.band);

  // fprintf(stderr, "%s[%d]:%s() %p %u %10.4f\n", __FILE__, __LINE__, __FUNCTION__, person, person->ref_num(), m_tcurr);

  if (kcurr != kprev || acurr != aprev) {  
    m_seek[gj][aprev][kprev].remove(person);
    m_seek[gj][acurr][kcurr].insert(person, Partnership::seek(person));

    m_supply[gj][aprev][kprev] = m_seek[gj][aprev][kprev].sum();
    m_supply[gj][acurr][kcurr] = m_seek[gj][acurr][kcurr].sum();
    m_supply_age[gj][aprev] = std::accumulate(m_supply[gj][aprev], m_supply[gj][aprev] + LEVELS, 0.0);
    m_supply_age[gj][acurr] = std::accumulate(m_supply[gj][acurr], m_supply[gj][acurr] + LEVELS, 0.0);
    m_supply_sex[gj] = std::accumulate(m_supply_age[gj], m_supply_age[gj] + AGES, 0.0);

    it->second.band = acurr;
    it->second.activity = kcurr;
  } else {
    assert(m_seek[gj][acurr][kcurr].member(person));
    m_seek[gj][acurr][kcurr].update(person, Partnership::seek(person));

    m_supply[gj][acurr][kcurr] = m_seek[gj][acurr][kcurr].sum();
    m_supply_age[gj][acurr] = std::accumulate(m_supply[gj][acurr], m_supply[gj][acurr] + LEVELS, 0.0);
    m_supply_sex[gj] = std::accumulate(m_supply_age[gj], m_supply_age[gj] + AGES, 0.0);
  }
}

double Partnership::seek(Person const* const person) {
  size_t ni(person->partners.size());
  double ri(person->seek * params.prop_seek_decr[person->stage()][person->treat()]);
  double di(person->damp);
  return ri * gsl_pow_int(di, ni);
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

void Partnership::calculate_mixing() {
  const double eps(std::numeric_limits<double>::denorm_min());

  const int W(FEMALE), M(MALE);
  double mixw, mixm;
  double uk, us, ud;
  for (int wa(0); wa < AGES; ++wa) {
    for (int wk(0); wk < LEVELS; ++wk) {
      for (int ma(0); ma < AGES; ++ma) {
	for (int mk(0); mk < LEVELS; ++mk) {
	  // men
	  us = (1 - params.assort_age) * m_supply_age[W][wa] / (eps + m_supply_sex[W]) + params.assort_age * (wa == ma);
	  uk = (1 - params.assort_act) * m_supply[W][wa][wk] / (eps + m_supply_age[W][wa]) + params.assort_act * (mk == wk);
	  if (ma == wa && ma > 0) {
	    mixm = (1 - params.assort_dif) * us * uk;
	  } else if (ma == wa + 1) {
	    ud = (1 - params.assort_age) * m_supply_age[W][wa+1] / (eps + m_supply_sex[W]) + params.assort_age;
	    mixm = (us + params.assort_dif * ud) * uk;
	  } else {
	    mixm = us * uk;
	  }

	  // women
	  us = (1 - params.assort_age) * m_supply_age[M][ma] / (eps + m_supply_sex[M]) + params.assort_age * (wa == ma);
	  uk = (1 - params.assort_act) * m_supply[M][ma][mk] / (eps + m_supply_age[M][ma]) + params.assort_act * (mk == wk);
	  if (ma == wa && wa + 1 < AGES) {
	    mixw = (1 - params.assort_dif) * us * uk;
	  } else if (ma == wa + 1) {
	    ud = (1 - params.assort_age) * m_supply_age[M][ma-1] / (eps + m_supply_sex[M]) + params.assort_age;
	    mixw = (us + params.assort_dif * ud) * uk;
	  } else {
	    mixw = us * uk;
	  }

	  m_mix[ma][mk][wa][wk] = mixm * mixw;
	}
      }
    }
  }
}

void Partnership::calculate_mixing_new() {
  const double eps(std::numeric_limits<double>::denorm_min());
  const int W(FEMALE), M(MALE);

  double mix_age[SEXES][AGES][AGES];
  double mix_act[SEXES][LEVELS][LEVELS][AGES];
  double mixw, mixm;

  for (int wa(0); wa < AGES; ++wa) {
    for (int ma(0); ma < AGES; ++ma) {
      mix_age[M][ma][wa] = (1 - params.assort_age) * m_supply_age[W][wa] / (eps + m_supply_sex[W]) + params.assort_age * (wa == ma);
      mix_age[W][wa][ma] = (1 - params.assort_age) * m_supply_age[M][ma] / (eps + m_supply_sex[M]) + params.assort_age * (wa == ma);
    }
  }

  for (int wk(0); wk < LEVELS; ++wk) {
    for (int mk(0); mk < LEVELS; ++mk) {
      for (int ai(0); ai < AGES; ++ai) {
	mix_act[M][mk][wk][ai] = (1 - params.assort_act) * m_supply[W][ai][wk] / (eps + m_supply_age[W][ai]) + params.assort_act * (mk == wk);
	mix_act[W][wk][mk][ai] = (1 - params.assort_act) * m_supply[M][ai][mk] / (eps + m_supply_age[M][ai]) + params.assort_act * (mk == wk);
      }
    }
  }

  for (int ma(0); ma < AGES; ++ma) {
    for (int mk(0); mk < LEVELS; ++mk) {
      for (int wa(0); wa < AGES; ++wa) {
	for (int wk(0); wk < LEVELS; ++wk) {
	  if (ma == wa && ma > 0) {
	    mixm = mix_act[M][mk][wk][wa] * (1 - params.assort_dif) * mix_age[M][ma][wa];
	  } else if (ma == wa + 1) {
	    mixm = mix_act[M][mk][wk][wa] * (mix_age[M][ma][wa] + params.assort_dif * mix_age[M][ma][wa+1]);
	  } else {
	    mixm = mix_act[M][mk][wk][wa] * mix_age[M][ma][wa];
	  }

	  if (ma == wa && wa + 1 < AGES) {
	    mixw = mix_act[W][wk][mk][ma] * (1 - params.assort_dif) * mix_age[W][wa][ma];
	  } else if (ma == wa + 1) {
	    mixw = mix_act[W][wk][mk][ma] * (mix_age[W][wa][ma] + params.assort_dif * mix_age[W][wa][ma-1]);
	  } else {
	    mixw = mix_act[W][wk][mk][ma] * mix_age[W][wa][ma];
	  }

	  m_mix[ma][mk][wa][wk] = mixm * mixw;
	}
      }
    }
  }
}

double Partnership::rate() {
  calculate_mixing_new();
  double rsum(0.0);

  for (int ma(0); ma < AGES; ++ma) {
    for (int mk(0); mk < LEVELS; ++mk) {
      for (int wa(0); wa < AGES; ++wa) {
	for (int wk(0); wk < LEVELS; ++wk) {
	  m_rate[ma][mk][wa][wk] = sqrt(m_mix[ma][mk][wa][wk] * m_supply[ MALE ][ma][mk] * m_supply[FEMALE][wa][wk]);
	  rsum += m_rate[ma][mk][wa][wk];
	}
      }
    }
  }
  return rsum;
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
  EventCache::const_iterator pi;
  for (pi = m_cache[AGE_INT].begin(); pi != m_cache[AGE_INT].end(); ++pi) {
    reassign(*pi);
  }
  m_cache[AGE_INT].clear();
  m_time[AGE_INT] = infinity;
}

void Partnership::risk_int() {
  EventCache::const_iterator pi;
  for (pi = m_cache[RISK_INT].begin(); pi != m_cache[RISK_INT].end(); ++pi) {
    reassign(*pi);
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
  return m_tcurr + gsl_ran_exponential(m_rng, 1.0 / rate());
}

Partnership::Match Partnership::mnext() const {
  const int n(AGES*AGES*LEVELS*LEVELS);
  Match m;
  double psum[n], grpr[n], U;
  unsigned int agem[n], actm[n], agew[n], actw[n];
  int i(0);
  for (int ma(0); ma < AGES; ++ma) {
    for (int mk(0); mk < LEVELS; ++mk) {
      for (int wa(0); wa < AGES; ++wa) {
	for (int wk(0); wk < LEVELS; ++wk) {
	  grpr[i] = m_rate[ma][mk][wa][wk];
	  agem[i] = ma;
	  actm[i] = mk;
	  agew[i] = wa;
	  actw[i] = wk;
	  ++i;
	}
      }
    }
  }

  psum[0] = grpr[0];
  for (int j(1); j < n; ++j) psum[j] = grpr[j] + psum[j-1];
  U = psum[n-1] * gsl_rng_uniform(m_rng);
  for (int j(0); j < n; ++j) {
    if (U < psum[j]) {
      // tnext() is optimized by assuming that the number of existing
      // partnerships is much smaller than the number of potential
      // pairings. We exploit this assumption by summing partnership
      // formation rates among all possible pairings without
      // subtracting the formation rates in existing partnerships. For
      // correct simulation, existing partnerships can be sampled and
      // rejected. Rejected pairings are marked by setting both
      // members to NULL.
      //
      // Note that if we instead resampled until we get a valid pair,
      // the algorithm would be incorrect because new partnerships
      // would form too rapidly.

      assert(m_seek[ MALE ][agem[j]][actm[j]].size() > 0);
      assert(m_seek[FEMALE][agew[j]][actw[j]].size() > 0);

      m.first  = m_seek[ MALE ][agem[j]][actm[j]].sample(m_rng);
      m.second = m_seek[FEMALE][agew[j]][actw[j]].sample(m_rng);
      if (m.first->partners.find(m.second) != m.first->partners.end()) {
	m.first = m.second = NULL;
      }
      return m;
    }
  }
  assert(false); // should never get here
}

void Partnership::consistent() const {
  for (int ai(0); ai < AGES; ++ai) {
    for (int ki(0); ki < LEVELS; ++ki) {
      if (m_seek[ MALE ][ai][ki].size() == 0) {assert(m_supply[ MALE ][ai][ki] == 0.0);}
      if (m_seek[FEMALE][ai][ki].size() == 0) {assert(m_supply[FEMALE][ai][ki] == 0.0);}
    }
  }
}
