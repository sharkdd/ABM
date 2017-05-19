#ifndef MESSAGE_H
#define MESSAGE_H

#include <ModelTypes.H>

struct Message {
  Message() : event_type(EVENTS) {}
  Message(EventType e, Person* p1 = NULL, Person* p2 = NULL) : event_type(e), person1(p1), person2(p2) {}
  ~Message() {}
  EventType event_type;
  Person* person1;
  Person* person2;
};

#endif // MESSAGE_H
