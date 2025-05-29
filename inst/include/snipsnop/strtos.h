#include <ctype.h> // pour isspace
#include <string>

#ifndef _strtos_
#define _strtos_

// comme strtof mais renvoie une std::string. On s'arrête au premier blanc.
inline std::string strtos(const char * x, char ** end) {
  const char * p = x;
  std::string val;
  // skip white spaces
  while(isspace(*p)) p++;

  while(!isspace(*p)) {
    val += *p;
    p++;
  }

  if(end != NULL) {
    *end = (char *) p;
  }
  return val;
}

#endif
