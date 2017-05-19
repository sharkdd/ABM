#include <numeric>
#include <Bookkeeper.H>

Bookkeeper::Bookkeeper(double interval, size_t ncollect)
  : m_tcurr(0.0), m_tnext(0.0), m_twait(interval), m_count(0),
    m_infected(0), m_sum_tdr(0), m_death_hiv(0), m_collect(ncollect), m_deaths(0),
    m_num_art(0), m_num_prep(0), m_t_art(0.0), m_t_prep(0.0), m_py_art(0.0), m_py_prep(0.0),
    m_report_pships(false) {
  for (unsigned int bi(0); bi < AGES; ++bi) {
    m_incident[MALE][bi] = m_incident[FEMALE][bi] = 0;
  }
}

Bookkeeper::~Bookkeeper() {
#ifdef PARTNER_HISTORY
  print_pship();
#endif // PARTNER_HISTORY
  for (Live::iterator li(live.begin()); li != live.end(); ++li) delete *li;
  for (Dead::iterator di(dead.begin()); di != dead.end(); ++di) delete *di;
}

inline double Bookkeeper::ta() {return m_tnext - m_tcurr;}

void Bookkeeper::delta_int() {
  // Print the population state and compute the next print time
  ++m_count;
  m_tcurr = m_tnext;
  m_tnext = m_count * m_twait;
  print_state();
}

void Bookkeeper::delta_ext(double dt, adevs::Bag<Message> const& msgs) {
  m_tcurr += dt;
  adevs::Bag<Message>::const_iterator mi;
  for (mi = msgs.begin(); mi != msgs.end(); ++mi) {
    switch ((*mi).event_type) {
    case DEBUT: debut_ext(*mi); break;
    case DEATH_NAT: death_nat_ext(*mi); break;
    case DEATH_HIV: death_hiv_ext(*mi); break;
    case INFECT: infect_ext(*mi); break;
    case PROGRESS: progress_ext(*mi); break;
    case ART_INIT: treat_init_ext(*mi); break;
    case ART_CHANGE: treat_change_ext(*mi); break;
    case PREP_INIT: prep_init_ext(*mi); break;
    case PREP_CHANGE: prep_change_ext(*mi); break;
    case RESISTANCE: resist_ext(*mi); break;
    default: break;
    }
  }
}

void Bookkeeper::delta_conf(adevs::Bag<Message> const& msgs) {
  delta_ext(0.0, msgs);
  delta_int();
}

inline void Bookkeeper::debut_ext(Message const& m) {
  live.insert(m.person1);
}

inline void Bookkeeper::death_nat_ext(Message const& m) {
  death_ext(m);
}

inline void Bookkeeper::death_hiv_ext(Message const& m) {
  ++m_death_hiv;
  death_ext(m);
}

inline void Bookkeeper::death_ext(Message const& m) {
  if (m.person1->prep() != PREP_NONE) {
    update_py_prep();
    --m_num_prep;
  }

  if (m.person1->treat() != ART_NONE) {
    update_py_art();
    --m_num_art;
  }

  live.erase(m.person1);
  dead.push_front(m.person1);
  ++m_deaths;
  deallocate();
}

inline void Bookkeeper::infect_ext(Message const& m) {
  Person *p(m.person2);
  ++m_incident[p->sex()][p->age_band()];
  m_sum_tdr += (p->virus() != WT);
}

inline void Bookkeeper::progress_ext(Message const& m) {
}

inline void Bookkeeper::treat_init_ext(Message const& m) {
  update_py_art();
  ++m_num_art;
}

inline void Bookkeeper::treat_change_ext(Message const& m) {
  update_py_art();
  m_num_art -= (m.person1->treat() == ART_NONE);
}

inline void Bookkeeper::prep_init_ext(Message const& m) {
  update_py_prep();
  ++m_num_prep;
}

inline void Bookkeeper::prep_change_ext(Message const& m) {
  update_py_prep();
  m_num_prep -= (m.person1->prep() == PREP_NONE);
}

inline void Bookkeeper::resist_ext(Message const& m) {
}

inline void Bookkeeper::update_py_art() {
  m_py_art += m_num_art * (m_tcurr - m_t_art);
  m_t_art = m_tcurr;
}

inline void Bookkeeper::update_py_prep() {
  m_py_prep += (m_tcurr - m_t_prep) * m_num_prep;
  m_t_prep = m_tcurr;
}


void Bookkeeper::print_header() const {
  std::cout << "year"
 	    << " X P1 P2 L1 L2 A1"
 	    << " incidence"
 	    << " sum.HIV sum.tdr sum.HIV.death"
 	    << " DR R1 R2 C1 C2 Q1 Q2"
 	    << " PrEP.X PrEP.P1 PrEP.P2 PrEP.L1 PrEP.L2 PrEP.A1"
 	    << " ART.P1 ART.P2 ART.L1 ART.L2 ART.A1"
 	    << " sum.ART sum.PrEP"
 	    << " prev.m1 prev.m2 prev.m3 prev.m4 prev.m5 prev.m6 prev.m7 prev.m8"
 	    << " prev.w1 prev.w2 prev.w3 prev.w4 prev.w5 prev.w6 prev.w7 prev.w8"
 	    << " inci.m1 inci.m2 inci.m3 inci.m4 inci.m5 inci.m6 inci.m7 inci.m8"
 	    << " inci.w1 inci.w2 inci.w3 inci.w4 inci.w5 inci.w6 inci.w7 inci.w8"
	    << " S.mk0 S.mk1 S.mk2 S.mk3 S.wk0 S.wk1 S.wk2 S.wk3"
	    << " PrEP.mk0 PrEP.mk1 PrEP.mk2 PrEP.mk3 PrEP.wk0 PrEP.wk1 PrEP.wk2 PrEP.wk3"
	    << " mmc.X mmc.Y";

  for (int ki(0); ki < LEVELS; ++ki) {
    for (int ai(0); ai < AGES; ++ai) std::cout << " nwk" << ki << 'b' << ai;
  }
  for (int ki(0); ki < LEVELS; ++ki) {
    for (int ai(0); ai < AGES; ++ai) std::cout << " nmk" << ki << 'b' << ai;
  }
  for (int ki(0); ki < LEVELS; ++ki) {
    for (int ai(0); ai < AGES; ++ai) std::cout << " ywk" << ki << 'b' << ai;
  }
  for (int ki(0); ki < LEVELS; ++ki) {
    for (int ai(0); ai < AGES; ++ai) std::cout << " ymk" << ki << 'b' << ai;
  }

  if (report_partnerships()) {
    std::cout << " couples";
    for (int wa(0); wa < AGES; ++wa) {
      for (int wk(0); wk < LEVELS; ++wk) {
	for (int ma(0); ma < AGES; ++ma) {
	  for (int mk(0); mk < LEVELS; ++mk) {
	    std::cout << " nw" << wa << wk << 'm' << ma << mk;
	  }
	}
      }
    }
  }

  std::cout << '\n';
}

void Bookkeeper::print_state() {
  update_py_art();  // get up-to-date person-years of ART
  update_py_prep(); // get up-to-date person-years of PrEP

  // number of persons in each stage of infection
  size_t X[STAGES];

  // number of persons, stratified by sex, activity, age and disease stae
  size_t N[SEXES][LEVELS][AGES][STAGES];
  size_t I[SEXES][LEVELS][AGES]; // number infected

  // Y[gi][bi][0]: number of gender-gi ageband-bi individuals susceptible
  // Y[gi][bi][1]: number of gender-gi ageband-bi individuals infected
  size_t Y[SEXES][AGES][2];

  // number of prevalent cases by virus. Susceptible persons have
  // virus=VARIANTS. We just count these persons in a dummy variable
  // instead of using conditional logic,
  size_t V[VARIANTS+1]; 

  // Count the number on ART by ART state and stage of infection
  size_t T[STAGES][ART_STATES];

  // used to record numbers of couples. couples[wa][wk][ma][mk] is
  // indexed by woman's age band (wa), woman's activity level (wk)
  // man's age band (ma) and man's activity level (mk)
  size_t ncouples(0);
  size_t couples[AGES][LEVELS][AGES][LEVELS];

  // Count the number of circumcised individuals (the array is
  // stratified by sex and disease stage so to avoid conditional logic
  // in the main accounting loop).
  // mmc[s][h][i]: i=0, uncircumcised; i=1, circumcised
  size_t mmc[SEXES][STAGES][2];

  // Count the numbers on PrEP
  size_t Q[SEXES][LEVELS][STAGES][PREP_STATES];
  size_t Qs[STAGES][PREP_STATES];

  size_t incident(0);
  double incidence;

  std::fill(V, V + VARIANTS + 1, 0);

  for (unsigned int si(0); si < SEXES; ++si) {
    for (unsigned int ai(0); ai < LEVELS; ++ai) {
      for (unsigned int bi(0); bi < AGES; ++bi) {
	for (unsigned int hi(0); hi < STAGES; ++hi) {
	  N[si][ai][bi][hi] = 0;
	}
      }
    }
  }

  if (report_partnerships()) {
    for (unsigned int wa(0); wa < AGES; ++wa) {
      for (unsigned int wk(0); wk < LEVELS; ++wk) {
	for (unsigned int ma(0); ma < AGES; ++ma) {
	  for (unsigned int mk(0); mk < LEVELS; ++mk) {
	    couples[wa][wk][ma][mk] = 0.0;
	  }
	}
      }
    }
  }

  for (unsigned int si(0); si < SEXES; ++si) {
    for (unsigned int hi(0); hi < STAGES; ++hi) {
      mmc[si][hi][0] = mmc[si][hi][1] = 0;
    }
  }

  for (unsigned int hi(0); hi < STAGES; ++hi) {
    for (unsigned int ti(0); ti < ART_STATES; ++ti) {
      T[hi][ti] = 0;
    }
  }

  for (unsigned int si(0); si < SEXES; ++si) {
    for (unsigned int ai(0); ai < LEVELS; ++ai) {
      for (unsigned int hi(0); hi < STAGES; ++hi) {
	for (unsigned int qi(0); qi < PREP_STATES; ++qi) {
	  Q[si][ai][hi][qi] = 0;
	}
      }
    }
  }

  // count the number of persons in key epidemiological states
  for (Live::const_iterator li(live.begin()); li != live.end(); ++li) {
    ++N[(*li)->sex()][(*li)->activity()][(*li)->age_band()][(*li)->stage()];
    ++V[(*li)->virus()];
    ++T[(*li)->stage()][(*li)->treat()];
    ++Q[(*li)->sex()][(*li)->activity()][(*li)->stage()][(*li)->prep()];
    ++mmc[(*li)->sex()][(*li)->stage()][(*li)->circumcised()];
  }

  if (report_partnerships()) {
    int wa, wk, ma, mk;
    for (Live::const_iterator li(live.begin()); li != live.end(); ++li) {
      // iterate through partner list
      // look up partner attributes
      // only looking at women because partnership information stored
      // among men is symmetric
      if ((*li)->sex() == FEMALE && !(*li)->partners.empty()) {
	wa = (*li)->age_band();
	wk = (*li)->activity();
	for (Person::Partners::const_iterator qi((*li)->partners.begin()); qi != (*li)->partners.end(); ++qi) {
	  assert((*qi)->sex() == MALE);
	  ma = (*qi)->age_band();
	  mk = (*qi)->activity();
	  ++couples[wa][wk][ma][mk];
	}
      }
    }

    for (unsigned int wa(0); wa < AGES; ++wa) {
      for (unsigned int wk(0); wk < LEVELS; ++wk) {
	for (unsigned int ma(0); ma < AGES; ++ma) {
	  for (unsigned int mk(0); mk < LEVELS; ++mk) {
	    ncouples += couples[wa][wk][ma][mk];
	  }
	}
      }
    }
  }

  for (unsigned int hi(0); hi < STAGES; ++hi) {
    X[hi] = 0;
    for (unsigned int si(0); si < SEXES; ++si) {
      for (unsigned int ai(0); ai < LEVELS; ++ai) {
	for (unsigned int bi(0); bi < AGES; ++bi) {
	  X[hi] += N[si][ai][bi][hi];
	}
      }
    }
  }

  for (unsigned int si(0); si < SEXES; ++si) {
    for (unsigned int bi(0); bi < AGES; ++bi) {
      Y[si][bi][0] = Y[si][bi][1] = 0;
      for (unsigned int ai(0); ai < LEVELS; ++ai) {
	I[si][ai][bi] = 0;
	Y[si][bi][0] += N[si][ai][bi][SUSCEPTIBLE];
	for (unsigned int hi(ACUTE_WINDOW); hi < STAGES; ++hi) {
	  I[si][ai][bi] += N[si][ai][bi][hi];
	}
	Y[si][bi][1] += I[si][ai][bi];
      }
    }
  }

  // calculate the numbers on or off PrEP in each disease stage
  for (unsigned int hi(0); hi < STAGES; ++hi) {
    for (unsigned int qi(0); qi < PREP_STATES; ++qi) {
      Qs[hi][qi] = 0;
      for (unsigned int si(0); si < SEXES; ++si) {
	for (unsigned int ai(0); ai < LEVELS; ++ai) {
	  Qs[hi][qi] += Q[si][ai][hi][qi];
	}
      }
    }
  }

  // increment the cumulative number of infections since the epidemic
  // began by the number of new infections in the previous time step
  for (unsigned int bi(0); bi < AGES; ++bi) {
    incident += m_incident[MALE][bi] + m_incident[FEMALE][bi];
  }
  m_infected += incident;

  // calculate incidence as a rate per reporting time step, then
  // divide by the time step duration to get annual
  // incidence. X[SUSCEPTIBLE] + m_infected approximates the number
  // susceptible at the beginning of the time step (may not be exact
  // because of recruitment and death during the time step).
  incidence = incident / static_cast<double>(X[SUSCEPTIBLE] + incident);
  incidence /= m_twait;

  std::cout << m_tcurr;

  for (unsigned int hi(0); hi < STAGES; ++hi) {
    std::cout << ' ' << X[hi];
  }
  //  std::cout << ' ' << incidence;
  std::cout << ' ' << incident;
  std::cout << ' ' << m_infected << ' ' << m_sum_tdr << ' ' << m_death_hiv;
  std::cout << ' ' << (V[R1] + V[R2] + V[C1] + V[C2] + V[Q1] + V[Q2])
    	    << ' ' << V[R1] << ' ' << V[R2]
	    << ' ' << V[C1] << ' ' << V[C2]
	    << ' ' << V[Q1] << ' ' << V[Q2];

  for (unsigned int hi(0); hi < STAGES; ++hi) {
    std::cout << ' ' << (Qs[hi][PREP_HIGH] + Qs[hi][PREP_POOR]);
  }

  for (unsigned int hi(ACUTE_WINDOW); hi < STAGES; ++hi) {
    std::cout << ' ' << (X[hi] - T[hi][ART_NONE]);
  }

  // print person-years of ART and person-years of PrEP
  std::cout << ' ' << m_py_art << ' ' << m_py_prep;

  // print age- and sex-stratified prevalence
  for (unsigned int si(MALE); si < SEXES; ++si) {
    for (unsigned int bi(0); bi < AGES; ++bi) {
      std::cout << ' ' << Y[si][bi][1] / static_cast<double>(Y[si][bi][0] + Y[si][bi][1]);
    }
  }
  
  // print age- and sex-stratified incidence
  for (unsigned int si(MALE); si < SEXES; ++si) {
    for (unsigned int bi(0); bi < AGES; ++bi) {
      incidence = m_incident[si][bi] / static_cast<double>(Y[si][bi][0] + m_incident[si][bi]);
      incidence /= m_twait;
      std::cout << ' ' << incidence;
    }
  }

  for (unsigned int si(0); si < SEXES; ++si) {
    for (unsigned int ki(0); ki < LEVELS; ++ki) {
      std::cout << ' ' <<
	Q[si][ki][SUSCEPTIBLE][PREP_NONE]
	+ Q[si][ki][SUSCEPTIBLE][PREP_HIGH]
	+ Q[si][ki][SUSCEPTIBLE][PREP_POOR];
    }
  }

  for (unsigned int si(0); si < SEXES; ++si) {
    for (unsigned int ki(0); ki < LEVELS; ++ki) {
      std::cout << ' ' << (Q[si][ki][SUSCEPTIBLE][PREP_HIGH] + Q[si][ki][SUSCEPTIBLE][PREP_POOR]);
    }
  }

  // fail if any women are circumcised
  for (unsigned int hi(0); hi < STAGES; ++hi) {assert(mmc[FEMALE][hi][1] == 0);}

  size_t mmc_y(0);
  for (unsigned int hi(ACUTE_WINDOW); hi < STAGES; ++hi) mmc_y += mmc[MALE][hi][1];
  std::cout << ' ' << mmc[MALE][SUSCEPTIBLE][1] << ' ' << mmc_y;

  // number of women
  for (unsigned int ki(0); ki < LEVELS; ++ki) {
    for (unsigned int ai(0); ai < AGES; ++ai) {
      std::cout << ' ' << std::accumulate(N[FEMALE][ki][ai], N[FEMALE][ki][ai] + STAGES, 0.0);
    }
  }
  // number of men
  for (unsigned int ki(0); ki < LEVELS; ++ki) {
    for (unsigned int ai(0); ai < AGES; ++ai) {
      std::cout << ' ' << std::accumulate(N[MALE][ki][ai], N[MALE][ki][ai] + STAGES, 0.0);
    }
  }
  // number of infected women
  for (unsigned int ki(0); ki < LEVELS; ++ki) {
    for (unsigned int ai(0); ai < AGES; ++ai) std::cout << ' ' << I[FEMALE][ki][ai];
  }
  // number of infected men
  for (unsigned int ki(0); ki < LEVELS; ++ki) {
    for (unsigned int ai(0); ai < AGES; ++ai) std::cout << ' ' << I[MALE][ki][ai];
  }

  if (report_partnerships()) {
    std::cout << ' ' << ncouples;
    for (unsigned int wa(0); wa < AGES; ++wa) {
      for (unsigned int wk(0); wk < LEVELS; ++wk) {
	for (unsigned int ma(0); ma < AGES; ++ma) {
	  for (unsigned int mk(0); mk < LEVELS; ++mk) {
	    std::cout << ' ' << couples[wa][wk][ma][mk];
	  }
	}
      }
    }
  }

  std::cout << '\n';

  // reset the incidence counter
  for (unsigned int bi(0); bi < AGES; ++bi) {
    m_incident[MALE][bi] = m_incident[FEMALE][bi] = 0;
  }
}

void Bookkeeper::deallocate() {
  if (collect() && (m_deaths >= collect())) {
    // slist::erase has linear time complexity. We use two slist
    // iterators and slist::erase_after (which has constant time
    // complexity) to more efficiently erase elements from the list.
    Dead::iterator di, dj(dead.begin());

    // Since deallocate is called when a person dies, since we are not
    // responsible for the order that modules resolve a person's death
    // (adevs handles that), and since new deaths are at the beginning
    // of Bookkeeper::dead, we skip directly over dead.begin() when
    // checking which Person objects may be deallocated.
    for (di = dead.begin(); di != dead.end();) {
      dj = di; ++dj;
      if (dj != dead.end() && !((*dj)->ref_num())) {
	delete *dj;
	dead.erase_after(di);
      } else {
	++di;
      }
    }

    // Reset the collection counter
    m_deaths = 0;
  }
}

void Bookkeeper::print_pship() const {
  int a, g, k;

  char sex[] = {'m', 'w'};

  std::cerr << "year";
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << " n" << sex[g] << 'b' << a << 'k' << k;
    }
  }
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << " s" << sex[g] << 'b' << a << 'k' << k;
    }
  }
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << " c" << sex[g] << 'b' << a << 'k' << k;
    }
  }
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << " t" << sex[g] << 'b' << a << 'k' << k;
    }
  }
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << " a" << sex[g] << 'b' << a << 'k' << k;
    }
  }
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << " m" << sex[g] << 'b' << a << 'k' << k;
    }
  }
  std::cerr << '\n';

  size_t number[SEXES][AGES][LEVELS]; // population size
  size_t single[SEXES][AGES][LEVELS]; // number of single individuals
  size_t concur[SEXES][AGES][LEVELS]; // number of individuals with concurrent partners 
  size_t couple[SEXES][AGES][LEVELS]; // cumulative partnerships in the past year
  size_t active[SEXES][AGES][LEVELS]; // number of individuals with >0 partners in past year
  size_t legion[SEXES][AGES][LEVELS]; // number of individuals with >1 parnters in past year

  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      std::fill(number[g][a], number[g][a] + LEVELS, 0);
      std::fill(single[g][a], single[g][a] + LEVELS, 0);
      std::fill(concur[g][a], concur[g][a] + LEVELS, 0);
      std::fill(couple[g][a], couple[g][a] + LEVELS, 0);
      std::fill(active[g][a], active[g][a] + LEVELS, 0);
      std::fill(legion[g][a], legion[g][a] + LEVELS, 0);
    }
  }

#ifdef PARTNER_HISTORY
  Person::History::const_iterator hi;
  int np;
  for (Live::const_iterator li(live.begin()); li != live.end(); ++li) {
    g = (*li)->sex();
    a = (*li)->age_band();
    k = (*li)->activity();

    ++number[g][a][k];

    single[g][a][k] += ((*li)->partners.size() == 0);
    concur[g][a][k] += ((*li)->partners.size() >= 2);

    // count number of partners in past year. this is the number of
    // current partners plus the number of breakup events in the past
    // year
    np = (*li)->partners.size();
    for (hi = (*li)->history.begin(); hi != (*li)->history.end(); ++hi) {
      if (m_tcurr - (*hi)->time > 1.0) break;
      np += ((*hi)->type() == BREAKUP);
    }

    couple[g][a][k] += np;
    active[g][a][k] += (np > 0);
    legion[g][a][k] += (np > 1);
  }
#endif // PARTNER_HISTORY

  std::cerr << m_tcurr;
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << ' ' << number[g][a][k];
    }
  }
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << ' ' << single[g][a][k];
    }
  }
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << ' ' << concur[g][a][k];
    }
  }
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << ' ' << couple[g][a][k];
    }
  }
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << ' ' << active[g][a][k];
    }
  }
  for (g = 0; g < SEXES; ++g) {
    for (a = 0; a < AGES; ++a) {
      for (k = 0; k < LEVELS; ++k) std::cerr << ' ' << legion[g][a][k];
    }
  }
  std::cerr << std::endl;
}
