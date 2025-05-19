#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory> // for shared_ptr
#include "datatype.h"

#ifndef _Column_
#define _Column_

struct Column {
  private:

    datatype type_ = NONE;
    std::shared_ptr<void> handler = nullptr;

    template<typename T>
    void checkType() const {
      if (whichType<T>() != type_) throw std::runtime_error("Trying to access with wrong type !");
    }

 public:

    // un constructeur qui prend comme argument un std vector et initialise handler et type_
    template <typename T>
    Column(std::vector<T> vec) : handler( std::make_shared<std::vector<T>>(vec) ) {
      type_ = whichType<T>();
    }

    // un constructeur qui prend juste un datatype et initialise avec un vecteur vide
    Column(datatype ty) : type_(ty) {
      if(ty == INT) {
        std::vector<int> vec;
        handler = std::make_shared<std::vector<int>>(vec);
      } else if(ty == FLOAT) {
        std::vector<float> vec;
        handler = std::make_shared<std::vector<float>>(vec);
      } else if(ty == DOUBLE) {
        std::vector<double> vec;
        handler = std::make_shared<std::vector<double>>(vec);
      } else if(ty == STRING) {
        std::vector<std::string> vec;
        handler = std::make_shared<std::vector<std::string>>(vec);
      } else {
        throw std::runtime_error("Can't initialize with type NONE");
      }
    }

    // getter pour type_
    datatype type() const {
      return type_;
    }

    // renvoie un pointer vers le std::vector...
    template <typename T>
    std::vector<T> * get() const {
      if (!handler) throw std::runtime_error("Handler was never instanciated");
      checkType<T>();
      return static_cast<std::vector<T>*>(handler.get());
    }

    // une fonction push_back 
    // obligation de faire un check_type... 
    template <typename T>
    void push_back(T x) {
      checkType<T>();
      ((std::vector<T> *) handler.get())->push_back(x);
    }

    // ~Column() {
    //     std::cout << "D° called for a Column of type " << type << " (0=INT,DOUBLE,FLOAT,3=STRING), before destr°, "<< handler.use_count() <<" ref to underliying data\n"; 
    // }

    // un constructeur qui fait une extraction
    // équivaut à vec[keep] en R
    // on template pour pouvoir utiliser n'importe quel type de vecteur d'entiers
    template<typename intVec>
    Column(const Column col, const intVec & keep) {
      datatype ty = col.type();
      if(ty == INT)
        *this = ColumnExtract<int>(col, keep);
      else if(ty == FLOAT)
        *this = ColumnExtract<float>(col, keep);
      else if(ty == DOUBLE)
        *this = ColumnExtract<double>(col, keep);
      else if(ty == STRING)
        *this = ColumnExtract<std::string>(col, keep);
      else
        throw std::runtime_error("Can't extract from type NONE");
    }

    // ceci fait le boulot pour le constructeur ci dessus
    // c'est privé puisqu'il suffit d'appeller le constructeur (qui n'a pas besoin
    // d'etre templaté pour le type de la colonne)
  private:
    template<typename T, typename intVec>
    Column ColumnExtract(const Column col, const intVec & keep) {
      const std::vector<T> * vec = col.get<T>();
      std::vector<T> filtered;
      filtered.reserve(keep.size());
      for(size_t i : keep)
        filtered.push_back(vec->at(i));
      return Column(filtered);
    }
};

#endif
