#include <algorithm>
#include <tuple>
#include <cmath>
#include "dpf.h"

int main(int argc, char ** argv)
{
  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  const size_t nitems = 1024;
  dpf_key<nitems> dpfkey[2];

  __m128i * etout = new __m128i[dpf_key<nitems>::etlength];
  __m128i ** s = new __m128i*[2]; 
    s[0] = new __m128i[nitems];
    s[1] = new __m128i[nitems];
  uint8_t ** t = new uint8_t*[2];
    t[0] = new uint8_t[nitems];
    t[1] = new uint8_t[nitems];
  
  gen(aeskey, 501, dpfkey);
  //evalfull(aeskey, dpfkey[0], etout, s, t);
  ///for (int i = 0; i < nitems/128; ++i) printf("%llx %llx ", etout[i][0], etout[i][1]); printf("\n");
  //evalfull(aeskey, dpfkey[0], etout, s[0], t[0]);
  //for (int i = 0; i < dpf_key<nitems>::etlength; ++i) printf("%llx %llx ", etout[i][0], etout[i][1]); printf("\n\n");
  //evalfull(aeskey, dpfkey[1], etout, s, t);
  //for (int i = 0; i < nitems/128; ++i) printf("%llx %llx ", etout[i][0], etout[i][1]); printf("\n");
  //evalfull(aeskey, dpfkey[1], etout, s[0], t[0]);
  ///for (int i = 0; i < dpf_key<nitems>::etlength; ++i) printf("%llx %llx ", etout[i][0], etout[i][1]); printf("\n");

  

  delete [] s[0]; delete [] s[1]; delete [] s;
  delete [] t[0]; delete [] t[1]; delete [] t;
  delete [] etout;

  return 0;
}