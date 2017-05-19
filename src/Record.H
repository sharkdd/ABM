#ifndef RECORD_H
#define RECORD_H

#include <ModelTypes.H>

// This is basically enabling run-time type inference on record
// types. There is likely a more idiomatic way of accomplishing this

// TODO: move derived event records to the modules that need them

// Abstract data type used to record events
struct EventRecord {
  EventRecord(double t = infinity) : time(t) {}
  virtual ~EventRecord() {}
  virtual EventType type() const = 0;
  float time;
};

struct DebutRecord : public EventRecord {
  DebutRecord(double t = infinity) : EventRecord(t) {}
  virtual ~DebutRecord() {}
  virtual EventType type() const {return DEBUT;}
};

struct ExitRecord : public EventRecord {
  ExitRecord(double t = infinity) : EventRecord(t) {}
  virtual ~ExitRecord() {}
  virtual EventType type() const {return DEATH_NAT;}
};

struct PartnerRecord : public EventRecord {
  PartnerRecord(double t = infinity, Person *p = NULL) : EventRecord(t), partner(p) {}
  virtual ~PartnerRecord() {}
  virtual EventType type() const {return PARTNER;}
  Person* partner; // newly acquired partner
};

struct BreakupRecord : public EventRecord {
  BreakupRecord(double t = infinity, Person *p = NULL) : EventRecord(t), partner(p) {}
  virtual ~BreakupRecord() {}
  virtual EventType type() const {return BREAKUP;}
  Person* partner; // other partner involved in newly dissolved partnership
};

struct ProgressRecord : public EventRecord {
  ProgressRecord(double t = infinity, Stage s = STAGES) : EventRecord(t), stage(s) {}
  virtual ~ProgressRecord() {}
  virtual EventType type() const {return PROGRESS;}
  Stage stage; // stage of infection resulting from the progression event
};

struct HIVDeathRecord : public EventRecord {
  HIVDeathRecord(double t = infinity) : EventRecord(t) {}
  virtual ~HIVDeathRecord() {}
  virtual EventType type() const {return DEATH_HIV;}
};

struct InfectRecord : public EventRecord {
  InfectRecord(double t = infinity, Person *p = NULL, Virus v = WT) : EventRecord(t), donor(p), virus(v) {}
  virtual ~InfectRecord() {}
  virtual EventType type() const {return INFECT;}
  Person* donor; // infection donor
  Virus virus; // virus present at transmission
};

struct TreatmentRecord : public EventRecord {
  TreatmentRecord(double t = infinity, Treat s = ART_STATES) : EventRecord(t), treat(s) {}
  virtual ~TreatmentRecord() {}
  virtual EventType type() const {return ART_CHANGE;}  
  Treat treat; // new treatment status
};

struct ProphylaxisRecord : public EventRecord {
  ProphylaxisRecord(double t = infinity, PrEP q = PREP_STATES) : EventRecord(t), prep(q) {}
  virtual ~ProphylaxisRecord() {}
  virtual EventType type() const {return PREP_CHANGE;}
  PrEP prep; // new PrEP status
};

struct ResistanceRecord : public EventRecord {
  ResistanceRecord(double t = infinity, Virus v = VARIANTS) : EventRecord(t), virus(v) {}
  virtual ~ResistanceRecord() {}
  virtual EventType type() const {return RESISTANCE;}
  Virus virus; // new virus
};

struct BehaviorRecord : public EventRecord {
  BehaviorRecord(double t = infinity, Activity k = LEVELS) : EventRecord(t), level(k) {}
  virtual ~BehaviorRecord() {}
  virtual EventType type() const {return RISK_CHANGE;}
  Activity level; // new sexual activity level
};

struct CircumcisionRecord : public EventRecord {
  CircumcisionRecord(double t = infinity) : EventRecord(t) {}
  virtual ~CircumcisionRecord() {}
  virtual EventType type() const {return MMC;}
};

#endif // RECORD_H
