#include "DataStruct.h"
#include "Column.h"
#include "datatype.h" // for typeToString
#include <iostream>
#include <stdexcept>// for runtime_error exception
#include <Rcpp.h>

// [[Rcpp::export]]
std::string getcolTypeDataStruct(Rcpp::XPtr<DataStruct> pDS, std::string colName) {
    std::string type;
    if (pDS->hasColumn(colName)){
    // copie pas génante parce que uniquement le ptr, pas vect de datas
    Column col = pDS->getColumn(colName);
    type = typeToString(col.type());
    } else {
    throw std::runtime_error("Error when trying to fetch a Column by name, it likely doesn't exist.");
    }
    return type;
}