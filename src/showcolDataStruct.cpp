#include <Rcpp.h>

#include <iostream>
#include <stdexcept> // for runtime error

#include "DataStruct.h"
#include "Column.h"
#include "datatype.h"

/* Function to help showing the first values of given (by name) Column in a 
DataStruct. 
Takes the number of characters left on the console line (I want it to fit into a single one)
*/

// [[Rcpp::export]]
std::string showcolDataStruct(Rcpp::XPtr<DataStruct> pDS, std::string colname, int size_left) {
  if (!pDS) throw std::runtime_error("Broken pointer for DataStruct");
  Column & col = pDS->getColumn(colname);

  std::string vals;
  size_t nvals = col.size();
  auto type = col.type();

  for (size_t i = 0; i < nvals; i++) {
    std::string newval = ""; // make sure it IS resetted at every loop
    switch (type) {
      case datatype::INT: {
        newval = std::to_string(col.at<int>(i));
        break;
      }
      case datatype::FLOAT: {
        newval = std::to_string(col.at<float>(i));
        break;
      }
      case datatype::DOUBLE: {
        newval = std::to_string(col.at<double>(i));
        break;
      }
      case datatype::STRING: {
        newval = col.at<std::string>(i);
        break;
      }
      case datatype::BOOL: {
        newval = ((col.at<bool>(i)) ? "TRUE" : "FALSE");
        break;
      }
      default:
        throw std::runtime_error("In show, type is NONE");
    }

    // choosing to always diplay the first value
    // + 1 why ?
    if (i != 0 && (vals.size() + newval.size() + 1) > size_left) {
      vals += "...";
      return vals;
    } else {
    vals += " ";
    vals += newval;
    }

  }
  #if DEBUG_COL
  vals += "Use count:";
  vals += col.handler_use_count();
  #endif
  return vals;
}