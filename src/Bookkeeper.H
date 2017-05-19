#ifndef BOOKKEEPER_H
#define BOOKKEEPER_H

#include <iostream>
#include <ext/slist>
#include <Person.H>
#include <Message.H>
#include <ModelTypes.H>
#include "adevs.h"

// The Bookkeeper object is responsible for reporting the simulation
// state at user-defined intervals. The Bookkeeper also maintains a
// reference to each living member of the population and is
// responsible for deallocating storage for Person objects. Dead
// individuals are also retained.
//
// Dead individuals are deallocated during simulation using a
// collection policy that periodically interrupts simulation. This
// collection policy attempts to deallocate people whenever at least
// Bookkeeper::collect deaths have occurred since the most recent
// deallocation. Dead individuals are deallocated only if their
// partner list is empty, since otherwise other modules may retain
// references to the person. If Bookkeeper::collect is zero, then dead
// individuals are deallocated only when the Bookkeeper is destroyed.
class Bookkeeper : public adevs::Atomic<Message> {
public:

  // Initialize with the reporting interval (in years) and the
  // deallocation policy parameter. Setting ncollect=0 postpones
  // deallocation of dead individuals until the Bookkeeper is
  // destroyed
  Bookkeeper(double interval = 1.0, size_t ncollect = 0);

  virtual ~Bookkeeper();

  // Establish the initial (living) population. 'Container' must store
  // Person* pointers.
  template<typename Container>  
  void initialize(Container const& people) {
    live.insert(people.begin(), people.end());
    typename Container::const_iterator it;
    for (it = people.begin(); it != people.end(); ++it) {
      if ((*it)->infected()) {
	++m_infected;
	m_sum_tdr += ((*it)->virus() != WT);
      }
    }
    print_header();
  }

  // adevs::Atomic<T> interface
  virtual double ta();
  virtual void delta_int();
  virtual void delta_ext(double dt, adevs::Bag<Message> const& msgs);
  virtual void delta_conf(adevs::Bag<Message> const& msgs);
  virtual void output_func(adevs::Bag<Message>& msgs) {}
  virtual void gc_output(adevs::Bag<Message>& garbage) {}

  // Accessors for the deallocation policy parameter. Set to zero to
  // delay deallocation until the Bookkeeper is destroyed, otherwise
  // the Bookkeeper will attempt to deallocate dead people whenever a
  // n*collect() deaths have occurred for integer values of n.
  inline size_t collect() const {return m_collect;}
  inline void collect(size_t const ncollect) {m_collect = ncollect;}

  // Access flags for partnership reporting. If partnership reporting
  // is enabled, then each state report includes the total number of
  // couples and the number of couples of each possible combination of
  // age bands and sexual activity levels. By default, partnership
  // reporting is disabled.
  //
  // Partnership reporting should be finalized before module
  // initialization, since otherwise the output header may be
  // misaligned with the data.
  inline bool report_partnerships() const {return m_report_pships;}
  inline void report_partnerships(bool const flag) {m_report_pships = flag;}

private:
  typedef std::tr1::unordered_set<Person*> Live;
  typedef __gnu_cxx::slist<Person*> Dead;

  // external event handlers
  void debut_ext(Message const& m);
  void death_nat_ext(Message const& m);
  void death_hiv_ext(Message const& m);
  void death_ext(Message const& m);
  void infect_ext(Message const& m);
  void progress_ext(Message const& m);
  void treat_init_ext(Message const& m);
  void treat_change_ext(Message const& m);
  void prep_init_ext(Message const& m);
  void prep_change_ext(Message const& m);
  void resist_ext(Message const& m);

  // helper functions that maintain state variables tracking
  // cumulative ARV delivery
  void update_py_art();
  void update_py_prep();

  void print_header() const;
  void print_state();

  // print partnership indicators
  void print_pship() const;

  // deallocates dead individuals if possible
  void deallocate();

  double m_tcurr; // bookkeeper internal clock
  double m_tnext; // next report time
  double m_twait; // time between consecutive reports
  size_t m_count; // counts elapsed reporting intervals

  size_t m_incident[SEXES][AGES]; // incidence counters

  size_t m_infected;  // cumulative cases of infection
  size_t m_sum_tdr;   // cumulative cases with transmitted resistant virus
  size_t m_death_hiv; // cumulative HIV deaths

  size_t m_collect; // deallocation policy parameter
  size_t m_deaths;  // counts the number of deaths since the most recent deallocation

  // variables used to track ARV delivery
  size_t m_num_art;  // number receiving ART
  size_t m_num_prep; // number receiving PrEP
  double m_t_art;    // time of last ART delivery update
  double m_t_prep;   // time of last PrEP delivery update
  double m_py_art;   // cumulative person-years of ART delivered
  double m_py_prep;  // cumulative person-years of PrEP delivered

  Live live; // Container for living individuals
  Dead dead; // Container for dead individuals

  bool m_report_pships; // partnership reporting flag
};

#endif // BOOKKEEPER_H
