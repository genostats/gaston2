#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory> // for shared_ptr
#include <istream>
#include "Column.h"

#ifndef _DataStruct_
#define _DataStruct_

#define DEBUG_DS false

class DataStruct {
  private:
    std::vector<Column> cols;
    std::vector<std::string> colNames;

    size_t findColumn(std::string name) const {
      for(size_t i = 0; i < colNames.size(); i++) {
        if(colNames[i] == name) return i;
      }
      return(colNames.size()); // out of range
      // throw std::out_of_range("Column not fund");
    }

  public:
    inline size_t size() const {
      return cols.size();
    }

    inline std::string colName(size_t i) const {
      return colNames.at(i);
    }

    // push_back = ajouter une colonne non nommée
    inline void push_back(const Column & newcol) { 
      cols.push_back(newcol);
      colNames.push_back("");
    }

    // ajouter une colonne nommée
    inline void push_back(const Column & newcol, std::string name) {
if(DEBUG_DS) std::cout << "Entering push_back\n";
      cols.push_back(newcol);
      colNames.push_back(name);
if(DEBUG_DS) std::cout << "Exiting push_back\n";
    }

    // ajouter toutes les colonnes de D
    inline void push_back(DataStruct D) {
      for(size_t i = 0; i < D.size(); i++) {
        push_back(D.at(i), D.colName(i));
      }
    }

    // vérifie si la colonne existe et la remplace, et sinon fait un push_back
    // a priori la colonne est copiée lors du passage de l'argument
    inline void setColumn(const Column & col, std::string name) {
if(DEBUG_DS) std::cout << "Entering setColumn\n";
      size_t pos = findColumn(name);
if(DEBUG_DS) std::cout << "pos = " << pos << "\n";
      if(pos == size()) {
if(DEBUG_DS) std::cout << "pushing\n";
        push_back(col, name);
      } else {
if(DEBUG_DS) std::cout << "replacing\n";
        cols[pos] = col;
      }
if(DEBUG_DS) std::cout << "Exiting setColumn\n";
    }

    // idem avec toutes les colonnes d'uutre data struct (forme simple de merge)
    inline void setColumns(const DataStruct & D) {
if(DEBUG_DS) std::cout << "Entering setColumns\n";
      for(size_t i = 0; i < D.size(); i++) {
if(DEBUG_DS) std::cout << "i = " << i << "\n";
        setColumn(D.at(i), D.colName(i));
      }
if(DEBUG_DS) std::cout << "Exiting setColumns\n";
    }

    inline const Column & at(size_t pos) const { 
      return cols.at(pos); 
    }

    inline const Column & getColumn(std::string name) const {
      return at(findColumn(name));
    }

    inline Column & at(size_t pos) { 
      return cols.at(pos); 
    }

    inline Column & getColumn(std::string name) {
      return at(findColumn(name));
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
if(DEBUG_DS) std::cout << "colonnes vides créées\n";
      // lecture fichier ligne par ligne
      std::string line;
      while(std::getline(in, line)) {
        char * c = (char *) line.c_str();
        for(size_t i = 0; i < nb_cols; i++) {
          c = cols[i].push_back_token(c);
          if(*c == 0) break; // fin de ligne
        }
      }
    }

};

#endif
