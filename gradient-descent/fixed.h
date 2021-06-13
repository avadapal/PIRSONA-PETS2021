/*
PIRSONA: A 1-private 4-party gradient descent algorithm to provide recommendations

Copyright (C) 2021  Adithya Vadapalli, Ryan Henry

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef FIXED_H__
#define FIXED_H__

#include <cstdint>
#include <cstring>
#include <cmath>

template <size_t precision, typename int_type = uint64_t>
struct fixed_t
{
public:
  fixed_t() { }
  fixed_t(double x) : val(std::llroundl(x * scale)) { }
  
   // template <size_t P>
   // fixed_t(fixed_t<P> other) : val(std::llroundl(static_cast<double>(other) * scale)) { }
  
   template <size_t P>
   fixed_t(fixed_t<P> other) : val(std::llroundl((double)(other) * scale)) { }






//   template <size_t P>
//   fixed_t(fixed_t<P> other) : val(std::llroundl(( other.val *  double((scale/ ((double)(1.0L * static_cast<double>(1ULL << P)))) )))) { }



  //template <size_t P>
  //fixed_t(fixed_t<P> other) : val(other.val >> (P  - precision)) { }
  
  //fixed_t(fixed_t<2*precision> other) : val(static_cast<double>(other) * scale) { }
  //inline operator double() const { return static_cast<double>((val &(msb-1)) / scale); }
  
  inline operator double() const { return (val&msb) ? -static_cast<double>(-val) / scale : static_cast<double>(val) / scale; }

  //inline operator double() const { return static_cast<double>((int_type)val) / scale; }
  inline operator int_type() const { return val; }
  inline fixed_t<precision> operator-() { return -val; }
  inline fixed_t<precision> & operator=(const double x) { val = std::llroundl(x * scale); return *this; }
  inline fixed_t<precision> operator+(const fixed_t<precision> & other) const   { return val + static_cast<int_type>(other); }
  inline fixed_t<precision> operator-(const fixed_t<precision> & other) const   { return val - static_cast<int_type>(other); }
  inline fixed_t<2*precision> operator*(const fixed_t<precision> & other) const { return val * static_cast<int_type>(other); }

  inline fixed_t<3*precision> operator*(const fixed_t<2*precision> & other) const { return val * static_cast<int_type>(other); }

  inline fixed_t<precision> & operator+=(const fixed_t<precision> & other)      { val += static_cast<int_type>(other); return *this; }
  inline fixed_t<precision> & operator-=(const fixed_t<precision> & other)      { val -= static_cast<int_type>(other); return *this; }
  inline fixed_t<precision> & operator*=(const fixed_t<precision> & other)      { val *= static_cast<int_type>(other) * scale; return *this; }

  inline fixed_t<precision> operator+(const double other) const { return val + other * scale; }
  inline fixed_t<precision> operator-(const double other) const { return val - other * scale; }
  inline fixed_t<precision> operator*(const double other) const { return val * other; }
  inline fixed_t<precision> & operator+=(const double other) { val += other * scale; return *this; }
  inline fixed_t<precision> & operator-=(const double other) { val -= other * scale; return *this; }
  inline fixed_t<precision> & operator*=(const double other) { val *= other; return *this; }

  inline fixed_t<precision> operator/(const double d) { return val * (1.0L / d); }
  inline fixed_t<precision> & operator/=(const double d) { val *= (1.0L / d); return *this; }




/*
  inline fixed_t<precision> operator+(int other) { return *this + static_cast<int_type>(other); }
  inline fixed_t<precision> operator-(int other) { return *this - static_cast<int_type>(other); }
  inline fixed_t<precision> & operator+=(int other) { return *this += static_cast<int_type>(other); }
  inline fixed_t<precision> & operator-=(int other) { return *this -= static_cast<int_type>(other); }

  inline fixed_t<precision> operator+(int_type other) { return val + other; }
  inline fixed_t<precision> operator-(int_type other) { return val - other; }
  inline fixed_t<precision> & operator+=(int_type other) { val += other; return val; }
  inline fixed_t<precision> & operator-=(int_type	 other) { val -= other; return val; } 
*/

  inline bool operator==(fixed_t<precision> & other)  { return val == static_cast<int_type>(other); }



  
  template <size_t P> inline bool operator!=(fixed_t<P> & other)  { return !(*this == other); }
  inline bool operator==(double other) { return (static_cast<int_type>(other * scale) & mask) == (val & mask); }
  inline bool operator!=(double other) { return !(*this == other); }
//private:
  fixed_t(int_type x) : val(x) {  }
  static constexpr double scale = 1.0L * static_cast<double>(1ULL << precision);
  static constexpr int_type msb = 1ULL << 63;
  static constexpr int_type mask = (1ULL << precision) - 1;
  int_type val;
};

#endif