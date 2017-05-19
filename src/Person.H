#ifndef PERSON_H
#define PERSON_H

#include <cassert>
#include <tr1/unordered_set>
#include <ext/slist>
#include <gsl/gsl_rng.h>
#include <ModelTypes.H>
#include <Record.H>

// Person objects are entities passed around by modules of the model.
// Person attributes dictate sexual behavior. Person objects store
// current sexual partners and event histories
struct Person {
public:

  // Data structure used to store the set of partnerships a person is
  // involved in
  class Partners {
  private:
    // Methods may be implemented inefficiently if the underlying
    // PartnerSet type is changed
    //    typedef std::tr1::unordered_set<Person*> PartnerSet;
    typedef __gnu_cxx::slist<Person*> PartnerSet;
    
  public:
    typedef PartnerSet::iterator iterator;
    typedef PartnerSet::const_iterator const_iterator;

    iterator begin() {return m_partners.begin();}
    iterator end() {return m_partners.end();}
    const_iterator begin() const {return m_partners.begin();}
    const_iterator end() const {return m_partners.end();}

    // Add a person as a partner
    inline void insert(Person* person) {
      m_partners.insert(m_partners.begin(), person);
    }

    // Remove a person as a partner
    inline void erase(Person* person) {
      m_partners.erase(find(person));
    }

    // Return an iterator to person if that person is a partner,
    // otherwise return Partners::end().
    inline iterator find(Person* person) {
      return std::find(begin(), end(), person);
    }

    inline const_iterator find(Person* person) const {
      return std::find(begin(), end(), person);
    }

    // Return the number of partners. This is guaranteed to
    // be O(n), where n is the number of partners.
    inline size_t size() const {return m_partners.size();} 

    // Return true if the partner set is empty, false otherwise
    inline bool empty() const {return m_partners.empty();}

    // Remove all partners
    inline void clear() {return m_partners.clear();}

  private:
    PartnerSet m_partners;
  };

#ifdef HISTORY
  // History wraps a container (singly-linked list) so that if we
  // decide change the underlying type all of the necessary interface
  // changes are concentrated here instead of scattered throughout the
  // codebase
  class History {
  private:
    typedef __gnu_cxx::slist<EventRecord*> RecordSet;
  public:
    typedef RecordSet::iterator iterator;
    typedef RecordSet::const_iterator const_iterator;

    iterator begin() {return m_history.begin();}
    iterator end() {return m_history.end();}
    const_iterator begin() const {return m_history.begin();}
    const_iterator end() const {return m_history.end();}

    // Record an event. Event histories are stored in reverse
    // chronological order to minimize the time and memory required
    // for history tracking
    inline void insert(EventRecord *record) {
      m_history.insert(m_history.begin(), record);
    }

    // Return the number of recorded events. This is guaranteed to be
    // O(n), where n is the number of events.
    inline size_t size() const {return m_history.size();}

    // Return true if there are no recorded events, false otherwise
    inline bool empty() const {return m_history.empty();}

    // Clear the event history
    inline void clear() {return m_history.clear();}
  private:
    RecordSet m_history;
  };
#else
  // If individual event histories are disabled we want to minimize
  // the amount of memory used. In principle we could just declare
  // Person::history inside compile guards. However, instead of
  // forcing modules to wrap history maintenance inside compile guards
  // we just set up a minimal data structure to safely ignore that
  // bookkeeping. The overhead appears to be 1 byte per person
  struct History {
    inline void insert(EventRecord *record) {delete record;}
  };
#endif // HISTORY

  // Default constructor. Person::alive() = true for people
  // constructed this way. Takes time of sexual debut as an input
  // parameter (this may be negative if the person reaches sexual
  // debut before simulation begins). All other attribute values are
  // undefined
  Person(double const tdebut = 0.0);

  // Default destructor
  ~Person();

  // Factory method to construct and initialize a person object. A
  // random number generator must be supplied to initialize some
  // attributes. tdebut, the time of sexual debut, must be supplied by
  // the caller.
  //
  // WARNING: this method assigns the newly constructed person to age
  // band 0. The caller must explicitly set the person's age band if
  // this is not the intended behavior (for example, if build_person
  // is used to construct people who debut before simulation begins).
  //
  // Persons initialized this way are susceptible and
  // ARV-naive. Gender and sexual activity level are initialized
  // randomly according to params.prop_female and
  // params.prop_activity(). See ModelTypes.H for details on these
  static Person* build_person(gsl_rng* rng, const double tdebut);

  inline bool susceptible() const {return stage() == SUSCEPTIBLE;}
  inline bool infected() const {return !susceptible();}

  // Check if the person has detectable HIV infection
  inline bool detectable() const {return stage() > ACUTE_WINDOW;}

  // return time of sexual debut
  inline double debut() const {return m_debut;}

  // Attribute accessors
  inline bool alive() const {return m_alive;}
  inline void alive(bool a) {m_alive = a;}

  inline Sex sex() const {return static_cast<Sex>(m_sex);}
  inline void sex(Sex s) {m_sex = s;}

  inline Activity activity() const {return static_cast<Activity>(m_activity);}
  inline void activity(Activity a) {m_activity = a;}

  inline unsigned int age_band() const {return m_band;}
  inline void age_band(unsigned int a) {m_band = a;}

  inline bool circumcised() const {return m_mmc;}
  inline void circumcised(bool c) {assert((sex() != FEMALE) | (!c)); m_mmc = c;}

  inline Stage stage() const {return static_cast<Stage>(m_stage);}
  inline void stage(Stage s) {m_stage = s;}

  inline Treat treat() const {return static_cast<Treat>(m_treat);}
  inline void treat(Treat t) {m_treat = t;}

  inline PrEP prep() const {return static_cast<PrEP>(m_prep);}
  inline void prep(PrEP p) {m_prep = p;}

  inline Virus virus() const {return static_cast<Virus>(m_virus);}
  inline void virus(Virus v) {m_virus = v;}

  // Records current partners. These are stored locally to a person
  // (as opposed to globally in a set of pairs of people) to simplify
  // finding all partnerships a given individual is involved in. This
  // is important, for example, when a person dies (since all
  // partnerships involving that person must be dissolved) or becomes
  // infected (since his or her partnerships may become
  // serodiscordant, introducing the risk of transmission).
  Partners partners;

  // Person record keeping
  History history;

  // Partnership parameters
  float seek; // partnership base seek rate
  float damp; // concurrency damping term

  // Reference count accessors used to determine when a person can be
  // deallocated. Epidemic modules that modify individuals' states are
  // responsible for calling ref_inc() when they begin tracking an
  // individual and must call ref_dec() when they stop tracking that
  // individual.
  inline unsigned int ref_num() const {return m_count;}
  inline void ref_inc() {++m_count;}
  inline void ref_dec() {assert(m_count > 0); --m_count;}
private:
  double m_debut; // time of sexual debut

  unsigned int m_alive    : 1;
  unsigned int m_sex      : 2; // gender. large enough to hold error value
  unsigned int m_activity : 3; // activity level
  unsigned int m_band     : 4; // age band
  unsigned int m_mmc      : 1; // true if male and circumcised, false otherwise
  unsigned int m_stage    : 4; // Infection status, recorded as infection stage
  unsigned int m_treat    : 3; // Treatment status
  unsigned int m_prep     : 2; // PrEP status
  unsigned int m_virus    : 4; // HIV variant if infected, undefined otherwise
  unsigned int m_count    : 16; // Reference count
};

#endif // PERSON_H
