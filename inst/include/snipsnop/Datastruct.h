#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory> // for shared_ptr
#include "Column.h"

#ifndef _DataStruct_
#define _DataStruct_

struct DataStruct {
    std::vector<Column> cols;
    std::vector<std::string> colNames;

    size_t size() const {
      return cols.size();
    }

    void push_back(Column newcol) { 
      cols.push_back(newcol);
      colNames.push_back("");
    }

    void push_back(Column newcol, std::string name) {
      cols.push_back(newcol);
      colNames.push_back(name);
    }

    // TODO 
    // function setColumn(Column, std::string) qui vérifie si la colonne existe et remplace, et sinon fait un push_back

    Column at(size_t pos) const { 
      return cols.at(pos); 
    }

    Column getColumn(std::string name) const {
      size_t pos = colNames.size(); // default value = out of range
      for(size_t i = 0; i < colNames.size(); i++) {
        if(colNames[i] == name) {
          pos = i;
          break;
        }
      }
      return at(pos);
    }

    DataStruct() {}

    template<typename intVec>
    DataStruct(DataStruct & DS, intVec & keep) {
      for(size_t i = 0; i < DS.size(); i++) {
        push_back(Column(DS.at(i), keep), DS.colNames[i]);
      }
    }
};

#endif
