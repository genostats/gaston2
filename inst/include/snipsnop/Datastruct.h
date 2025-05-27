#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory> // for shared_ptr
#include <istream>
#include "Column.h"

#ifndef _DataStruct_
#define _DataStruct_

class DataStruct {
  private:
    std::vector<Column> cols;
    std::vector<std::string> colNames;

  public:
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

    std::string colName(size_t i) const {
      return colNames[i];
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

    // constructeur vide...
    DataStruct() {}

    // constructeur qui fait une extraction de lignes
    template<typename intVec>
    DataStruct(DataStruct & DS, intVec & keep) {
      for(size_t i = 0; i < DS.size(); i++) {
        push_back(Column(DS.at(i), keep), DS.colNames[i]);
      }
    }

    // constructeur qui lit un fichier sans header
    // il faut donner les types des colonnes et leur nom
    // pas très perfectionné (uniquement space delimited values,
    // pas de prise en compte de quotes, pas bcp de check de format...) 
    // mais ça ira pour lire des fichiers bim / fam
    DataStruct(std::istream & in, std::vector<datatype> & colTypes, std::vector<std::string> & colNames) {
      size_t nb_cols = colTypes.size();
      if(nb_cols != colNames.size()) 
        throw std::runtime_error("colTypes and colNames have different sizes");
      // créations colonnes vides
      for(size_t i = 0; i < nb_cols; i++) 
        push_back(Column(colTypes[i]), colNames[i]);

      // lecture fichier ligne par ligne
      std::string line;
      while(std::getline(in, line)) {
        char * c = (char *) line.c_str();
        for(size_t i = 0; i < nb_cols; i++) {
          c = at(i).push_back_token(c);
          if(*c == 0) break; // fin de ligne
        }
      }
    }

};

#endif
