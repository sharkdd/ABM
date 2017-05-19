#ifndef RUNNINGSUM_H
#define RUNNINGSUM_H

#include <iostream>

// RunningSum defines an adaptor for numerically robust running sums
// of a given precision. T must be a floating point type.
template<typename T>
class RunningSum {
public:

  RunningSum(T value = 0.0) : m_value(value), m_error(0.0) {}
  RunningSum(RunningSum<T> const &x) : m_value(x.m_value), m_error(x.m_error) {}

  ~RunningSum() {}

  RunningSum &operator=(RunningSum<T> const &x) {m_value = x.m_value; m_error = x.m_error; return *this;}
  RunningSum &operator=(T const x) {m_value = x; m_error = 0.0; return *this;}

  bool operator==(RunningSum<T> const &x) const {return (m_value == x.m_value) && (m_error == x.m_error);}
  bool operator==(T const x) const {return (m_value == x) && (m_error == 0.0);}

  bool operator<(RunningSum<T> const &x) const {return m_value < x.m_value;}
  bool operator<(T const x) const {return m_value < x;}

  bool operator>(RunningSum<T> const &x) const {return m_value > x.m_value;}
  bool operator>(T const x) const {return m_value > x;}

  RunningSum<T>& operator+=(T const x) {
    const double y(x - m_error);
    const double t(m_value + y);
    m_error = (t - m_value) - y;
    m_value = t;
    return *this;
  }

  RunningSum<T>& operator-=(T const x) {return this->operator+=(-x);}

  operator T() const {return static_cast<T>(m_value);}

  friend std::ostream& operator<<(std::ostream& out, RunningSum<T> const &x) {
    out << x.m_value; return out;
  }

  friend std::istream& operator>>(std::istream& in, RunningSum<T> &x) {
    in >> x.m_value; x.m_error = 0.0; return in;
  }

private:
  T m_value;
  T m_error;
};

#endif // RUNNINGSUM_H
