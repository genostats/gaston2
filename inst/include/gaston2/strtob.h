#include <ctype.h> // pour isspace et isdigit
#include <cstring> // pour strcmp
#include <algorithm> // pour std::transform (voir si pas réécrivable)
#include <string> // for std::string
#include <stdexcept> // for std::invalid_argument


#ifndef _strtob_
#define _strtob_

// prends une str de C et tente de lire un bool
// je considère comme bool : 'TRUE'/'FALSE', 'true'/'false', '0'/'1'||digit
inline bool strtob(const char * x, char ** end) {
  const char * p = x;
  bool val = false;
  bool found = false; // to check digits OR strings
  // skip white spaces
  while(isspace(*p)) p++;

  // will consider 0 (false), 1/2/3/4/5/6/7/8/9 (true)
  while(isdigit(*p)) { // isdigit('\0') == false, so no need to check for NUL
    if (!found) {
      // TODO : see if this fixes the "0" == true
    val = (*p != '0');
    found = true;
    }
    p++;
  }

  if (!found){
    char * copy_true = strndup(p, 4); // avoiding pb if trailing whitespace 
    // OR for truefalsefalsetrue : will isolate first
    char * copy_false = strndup(p, 5);
    if (strcmp("true", copy_true) == 0 || strcmp("TRUE", copy_true) == 0) {
      val = true;
      found = true;
      p += 4; // 4 chars in true
    }
    else if (strcmp("false", copy_false) == 0 || strcmp("FALSE", copy_false) == 0 ) {
      val = false;
      found = true;
      p += 5;
    }
  }

  if(end != NULL) {
    *end = (char *) p;
  }

  // if not digit AND not str
  if (!found) throw std::invalid_argument("Failed to convert a char * to a bool.");

  return val;
}

// Prends une std::string et tente de lire un bool
inline bool stob(const std::string& x) {
    
  // skip first whitespaces
    const char * p = &x[0];
    size_t p_count = 0;
    while(isspace(*p)) {
      p++;
      p_count++;
    }
    // lower can start at the real value
    std::string lower = x.substr(p_count, x.length());
    
    // TODO : vérifier que ne casse rien si c'est des int !
    if (lower.rfind("true", 4) == 0 || lower.rfind("TRUE", 4) == 0 || lower.rfind("1", 1) == 0) {
        return true;
    } else if (lower.rfind("false", 5) == 0 || lower.rfind("FALSE", 5) == 0 || lower.rfind("0", 1) == 0) {
        return false;
    } else
        throw std::invalid_argument("Failed to convert this string to a bool: " + lower);
}

#endif