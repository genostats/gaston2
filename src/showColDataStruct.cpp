#include <Rcpp.h>

#include <iostream>
#include <stdexcept> // for runtime error

#include "DataStruct.h"
#include "getDataStructPtr.h"
#include "Column.h"
#include "datatype.h"

#include "gastonOptions.h" // definition of na_value

/* Function to help showing the first values of given (by name) Column in a 
DataStruct. 
Takes the number of characters left on the console line (I want it to fit into a single one)
*/

// [[Rcpp::export]]
std::string showColDataStruct(Rcpp::S4 x, std::string colname, int size_left) {
  DataStruct * pDS = getDataStructPtr(x);
  if (!pDS) Rcpp::stop("Broken pointer for DataStruct");
  Column & col = pDS->getColumn(colname);

  std::string vals;
  size_t nvals = col.size();
  auto type = col.type();

  for (size_t i = 0; i < nvals; i++) {
    std::string newval = "";
    switch (type) {
      case datatype::INT: {
        int v = col.at<int>(i);
        if(v == na_value) 
          newval = "NA";
        else
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
    // + 1 for the space
    if (i != 0 && (vals.size() + newval.size() + 1) > size_left) {
      vals += "...";
      return vals;
    } else {
      vals += " ";
      vals += newval;
    }

  }
  return vals;
}
