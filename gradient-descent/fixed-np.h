#ifndef FIXED_H__
#define FIXED_H__

#include <cstdint>
#include <cstring>
#include <cmath>

template <size_t precision>
struct fixed_t
{
public:
  fixed_t() { }
  fixed_t(double x) : val(std::llroundl(x * scale)) { }
  template <size_t P>
  fixed_t(fixed_t<P> other) : val(std::llroundl(static_cast<double>(other) * scale)) { }
  //fixed_t(fixed_t<2*precision> other) : val(static_cast<double>(other) * scale) { }
  //inline operator double() const { return (val&msb) ? -static_cast<double>(-val) / scale : static_cast<double>(val) / scale; }
  inline operator double() const { return static_cast<double>((int64_t)val) / scale; }
  inline operator uint64_t() const { return val; }
  inline fixed_t<precision> operator-() { return -val; }
  inline fixed_t<precision> & operator=(const double x) { val = std::llroundl(x * scale); return *this; }
  inline fixed_t<precision> operator+(const fixed_t<precision> & other) const   { return val + static_cast<uint64_t>(other); }
  inline fixed_t<precision> operator-(const fixed_t<precision> & other) const   { return val - static_cast<uint64_t>(other); }
  inline fixed_t<2*precision> operator*(const fixed_t<precision> & other) const { return val * static_cast<uint64_t>(other); }

inline fixed_t<3*precision> operator*(const fixed_t<2*precision> & other) const { return val * static_cast<uint64_t>(other); }

  inline fixed_t<precision> & operator+=(const fixed_t<precision> & other)      { val += static_cast<uint64_t>(other); return *this; }
  inline fixed_t<precision> & operator-=(const fixed_t<precision> & other)      { val -= static_cast<uint64_t>(other); return *this; }
  inline fixed_t<precision> & operator*=(const fixed_t<precision> & other)      { val *= static_cast<uint64_t>(other) * scale; return *this; }

  inline fixed_t<precision> operator+(const double other) const { return val + other * scale; }
  inline fixed_t<precision> operator-(const double other) const { return val - other * scale; }
  inline fixed_t<precision> operator*(const double other) const { return val * other; }
  inline fixed_t<precision> & operator+=(const double other) { val += other * scale; return *this; }
  inline fixed_t<precision> & operator-=(const double other) { val -= other * scale; return *this; }
  inline fixed_t<precision> & operator*=(const double other) { val *= other; return *this; }

  inline fixed_t<precision> operator/(const double d) { return val * (1.0L / d); }
  inline fixed_t<precision> & operator/=(const double d) { val *= (1.0L / d); return *this; }


/*
  inline fixed_t<precision> operator+(int other) { return *this + static_cast<uint64_t>(other); }
  inline fixed_t<precision> operator-(int other) { return *this - static_cast<uint64_t>(other); }
  inline fixed_t<precision> & operator+=(int other) { return *this += static_cast<uint64_t>(other); }
  inline fixed_t<precision> & operator-=(int other) { return *this -= static_cast<uint64_t>(other); }

  inline fixed_t<precision> operator+(uint64_t other) { return val + other; }
  inline fixed_t<precision> operator-(uint64_t other) { return val - other; }
  inline fixed_t<precision> & operator+=(uint64_t other) { val += other; return val; }
  inline fixed_t<precision> & operator-=(uint64_t	 other) { val -= other; return val; } 
*/

  inline bool operator==(fixed_t<precision> & other)  { return val == static_cast<uint64_t>(other); }



  
  template <size_t P> inline bool operator!=(fixed_t<P> & other)  { return !(*this == other); }
  inline bool operator==(double other) { return (static_cast<uint64_t>(other * scale) & mask) == (val & mask); }
  inline bool operator!=(double other) { return !(*this == other); }
//private:
  fixed_t(uint64_t x) : val(x) {  }
  static constexpr double scale = 1.0L * static_cast<double>(1ULL << precision);
//  static constexpr uint64_t msb = 1ULL << 63;
  static constexpr uint64_t mask = (1ULL << precision) - 1;
  uint64_t val;
};

#endif