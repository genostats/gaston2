#include <ctype.h> // pour isspace et isdigit

#ifndef _strtoi_
#define _strtoi_

// comme strtof mais renvoie un entier
inline int strtoi(const char * x, char ** end) {
  const char * p = x;
  bool sign = false;
  int val = 0;
  // skip white spaces
  while(isspace(*p)) p++;

  if(*p == '-') {
    sign = true;
    p++;
  } else if(*p == '+') {
    p++;
  }

  while(isdigit(*p)) {
    val *= 10;
    val += (*p - '0');
    p++;
  }

  if(sign) val *= -1;

  if(end != NULL) {
    *end = (char *) p;
  }
  return val;
}

#endif
