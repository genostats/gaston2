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
  
  /* from cppreferences isdigit() :
  "To use these functions safely with plain chars (or signed chars), 
  the argument should first be converted to unsigned char.
  so TODO : maybe cast p into static_cast<unsigned char>(*p)
  */
  while(isdigit(*p)) { // isdigit('\0') == false, so no need to check for NUL
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
