#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory> // for shared_ptr
#include <istream>
#include "Column.h"

#include "debug_flags.h"

#ifndef _DataStruct_
#define _DataStruct_

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
    // nombre de colonnes : deux fonctions
    inline size_t size() const {
      return cols.size();
    }

    inline size_t ncol() const {
      return cols.size();
    }

    // nombre de lignes
    inline size_t nrow() const {
      if(size() == 0) return(0);
      return cols[0].size();
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
#if DEBUG_DS
  std::cout << "Entering push_back\n";
#endif
      cols.push_back(newcol);
      colNames.push_back(name);
#if DEBUG_DS
  std::cout << "Exiting push_back\n";
#endif
    }

    // ajouter toutes les colonnes de D
    inline void push_back(DataStruct D) {
      for(size_t i = 0; i < D.size(); i++) {
        push_back(D.at(i), D.colName(i));
      }
    }

    // ---- setColumn ----
    // vérifie si la colonne existe et la remplace, et sinon fait un push_back
    // a priori la colonne est copiée lors du passage de l'argument
    // /!\ si c'est une colonne avec le nom == "", va remplacer la 1ère colonne sans nom !
    inline void setColumn(const Column & col, std::string name) {
#if DEBUG_DS
  std::cout << "Entering setColumn\n";
#endif
      size_t pos = findColumn(name);
#if DEBUG_DS
  std::cout << "pos = " << pos << "\n";
#endif
      if(pos == size()) {
#if DEBUG_DS
  std::cout << "pushing\n";
#endif
        push_back(col, name);
      } else {
#if DEBUG_DS
  std::cout << "replacing\n";
#endif
        cols[pos] = col;
      }
#if DEBUG_DS
  std::cout << "Exiting setColumn\n";
#endif
    }

    // ---- setColumns ----
    // idem avec toutes les colonnes d'un autre data struct (forme simple de merge)
    inline void setColumns(const DataStruct & D) {
#if DEBUG_DS
  std::cout << "Entering setColumns\n";
#endif
      for(size_t i = 0; i < D.size(); i++) {
#if DEBUG_DS
  std::cout << "i = " << i << "\n";
#endif
        setColumn(D.at(i), D.colName(i));
      }
#if DEBUG_DS
  std::cout << "Exiting setColumns\n";
#endif
    }

    // ---- get column at ----
    inline const Column & at(size_t pos) const { 
      return cols.at(pos); 
    }

    // non const version
    inline Column & at(size_t pos) { 
      return cols.at(pos); 
    }

    // ---- get column by name ----
    inline const Column & getColumn(std::string name) const {
      return at(findColumn(name));
    }

    // non const version
    inline Column & getColumn(std::string name) {
      return at(findColumn(name));
    }

    // ---- has column ? ----
    inline bool hasColumn(std::string name) {
      size_t i = findColumn(name);
      return i < size();
    }

    // ---- constructeurs ----
    // constructeur vide...
    DataStruct() {}

    // constructeur qui fait une extraction de lignes
    template<typename intVec>
    DataStruct(DataStruct & DS, intVec & keep) {
      for(size_t i = 0; i < DS.size(); i++) {
        push_back(Column(DS.at(i), keep), DS.colNames[i]);
      }
    }

    // constructeur qui créée des colonnes vides d'un type donné
    DataStruct(std::vector<datatype> & colTypes, std::vector<std::string> & colNames) {
      size_t nb_cols = colTypes.size();
      if(nb_cols != colNames.size()) 
        throw std::runtime_error("colTypes and colNames have different sizes");
      // créations colonnes vides
      for(size_t i = 0; i < nb_cols; i++) 
        push_back(Column(colTypes[i]), colNames[i]);
#if DEBUG_DS
  std::cout << "colonnes vides créées\n";
#endif
    }

    // I am stupid and this is not usefull for now but could be later on
    // // fonction qui enlève les colonnes (ainsi que leurs noms)
    // // spécifiés dans to_remove
    // void remove(DataStruct & D, std::vector<std::string> & to_remove) {
    //   for (const auto& name : to_remove) {
    //     size_t pos = D.findColumn(name);
    //     if (pos != D.size()) {
    //       D.cols.erase(D.cols.begin() + pos);
    //       D.colNames.erase(D.colNames.begin() + pos);
    //     } else {
    //       // commented because it is bound to happen 
    //         //throw std::runtime_error("Failed to remove Column " + name);
    //     }
    //   }
    // }

    // ---- rbind two DataStruct ----
    DataStruct(const DataStruct& dt1, const DataStruct& dt2) {
      size_t nb_cols = dt1.size();
      size_t nb_cols_dt2 = dt2.size();

      for(size_t i = 0; i < nb_cols; i++) {
        std::string name = dt1.colName(i);
        // if name is "", skipping it
        if (name != "") {
          // only keep and append Column that exists in both
          size_t i2 = dt2.findColumn(name);
          if(i2 != nb_cols_dt2) {
            // call constructor of Column that will merge the 2 columns
            push_back(Column(dt1.at(i), dt2.at(i2)), name);
          }
            // if a Column with this name 
            // doesn't exist in dt2, then not adding it
        }
      }
    }

    // ---- read File ----
    //
    // fonction qui lit (append) un fichier sans header dans des colonnes
    // déjà créées
    // pas très perfectionné (uniquement space delimited values,
    // pas de prise en compte de quotes, pas bcp de check de format...) 
    // mais ça ira pour lire des fichiers bim / fam
    void readFile(std::istream & in) {
      size_t nb_cols = size();
      // lecture fichier ligne par ligne
      std::string line;
      while(std::getline(in, line)) {
        char * c = (char *) line.c_str();
        for(size_t i = 0; i < nb_cols; i++) {
          c = cols[i].push_back_token(c); // toute la magie est dans ce membre de Column
          if(*c == 0) break; // fin de ligne
        }
      }
    }

};

#endif
