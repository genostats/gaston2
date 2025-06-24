#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory> // for shared_ptr
#include "datatype.h"
#include <stdlib.h>
#include "strtoi.h"
#include "strtos.h"

#include "debug_flags.h"

#ifndef _Column_
#define _Column_

struct Column {
  private:

    datatype type_ = NONE;
    std::shared_ptr<void> handler = nullptr;

    template<typename T>
    void checkType() const {
      if (whichType<T>() != type_) throw std::runtime_error("Trying to access with wrong type !");
      if(type_ == NONE) throw std::runtime_error("Trying to access column with type NONE"); // est-ce que c'est possible que ça arrive ?
    }

    void checkHandler() const {
      if (!handler) throw std::runtime_error("Handler was not instanciated");
    }

 public:

    // un constructeur qui prend comme argument un std vector et initialise handler et type_
    // le vecteur est copié lors du passage de l'argiment
    template <typename T>
    Column(std::vector<T> vec) : handler( std::make_shared<std::vector<T>>(vec) ) {
      type_ = whichType<T>();
      checkHandler();
    }

    // size of the column
    size_t size() const {
      switch(type_) {
        case INT : { 
           return ((std::vector<int> *) handler.get())->size();
        }
        case FLOAT : {
           return ((std::vector<float> *) handler.get())->size();
        } 
        case DOUBLE : {
           return ((std::vector<double> *) handler.get())->size();
        }
        case STRING : {
           return ((std::vector<std::string> *) handler.get())->size();
        }
        default :
          throw std::runtime_error("In size, type is NONE");
      }
    }


#if DEBUG_COL
    // verbose copy constructor
    Column(const Column & col) : type_(col.type_), handler(col.handler) {
      std::cout << "Copy one column of type " << typeToString(type_) << " (passed by ref)\n";
    }

    // verbose assignement operator    
    Column & operator=(const Column& col) {
      type_ = col.type_;
      handler = col.handler;
      std::cout << "Copy operator on one column of type " << typeToString(type_) << " (passed by ref)\n";
      return *this;
    }

    // verbose move operator, moving TO col
    Column & operator=(Column&& col)
    {
      std::cout << "Using move for column of type " << typeToString(type_) << "\n";
      std::swap(type_, col.type_);
      std::swap(handler, col.handler);
      return *this;
    }
#endif


    // un constructeur qui prend juste un datatype et initialise avec un vecteur vide
    Column(datatype ty) : type_(ty) {
      switch(ty) {
        case INT : { 
            std::vector<int> vec;
            handler = std::make_shared<std::vector<int>>(vec);
            break;
          }
        case FLOAT : {
          std::vector<float> vec;
          handler = std::make_shared<std::vector<float>>(vec);
          break;
        } 
        case DOUBLE : {
          std::vector<double> vec;
          handler = std::make_shared<std::vector<double>>(vec);
          break;
        }
        case STRING : {
          std::vector<std::string> vec;
          handler = std::make_shared<std::vector<std::string>>(vec);
          break;
        }
        default :
          throw std::runtime_error("Can't initialize with type NONE");
      }
      checkHandler();
    }

    // un constructeur qui fait une extraction
    // équivaut à vec[keep] en R
    // on template pour pouvoir utiliser n'importe quel type de vecteur d'entiers
    template<typename intVec>
    Column(const Column col, const intVec & keep) {
      datatype ty = col.type();
      switch(ty) {
        case INT : { 
          *this = ColumnExtract<int>(col, keep);
           break;
        }
        case FLOAT : {
          *this = ColumnExtract<float>(col, keep);
          break;
        } 
        case DOUBLE : {
          *this = ColumnExtract<double>(col, keep);
          break;
        }
        case STRING : {
          *this = ColumnExtract<std::string>(col, keep);
          break;
        }
        default :
          throw std::runtime_error("Can't extract from type NONE");
      }
    }

    // ceci fait le boulot pour le constructeur ci dessus
    // c'est privé puisque l'utilisateur peut appeller le constructeur (qui n'a pas besoin
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

    public : 
    // un constructeur qui retourne une copie
    // d'un append de col2 dans col
    Column(const Column col1, const Column col2) {
      // je pars du principe que les Columns doivent être du même type
      // à voir si possible pour chr ?
      type_ = col1.type();
      if (type_ != col2.type()) {
        throw std::runtime_error("Failing to append columns beacause of mismatched types");
      }

      switch (type_) {
        case INT: {
          *this = ColumnAppend<int>(col1, col2);
          break;
        }
        case FLOAT: {
          *this = ColumnAppend<float>(col1, col2);
          break;
        }
        case DOUBLE: {          
          *this = ColumnAppend<double>(col1, col2);
          break;
        }
        case STRING: {          
          *this = ColumnAppend<std::string>(col1, col2);
          break;
        }
        default:
          throw std::runtime_error("Cannot append columns of type NONE");
      }
    }

    private :
    // /!\ if vec1 has the same ptr as vec2 it will cause big pb :/
    // but this should never happen as they are not handled with ptr
    template <typename T>
    Column ColumnAppend(const Column& col1, const Column& col2) {
      const std::vector<T> * vec1 = col1.get<T>();
      const std::vector<T> * vec2 = col2.get<T>();
      std::vector<T> vec(*vec1);
      vec.insert(vec.end(), vec2->begin(), vec2->end());
      return Column(vec);
    }

  public:
    // getter pour type_
    datatype type() const {
      return type_;
    }

    // renvoie un pointer vers le std::vector...
    template <typename T>
    std::vector<T> * get() const {
      // a priori si handler = null alors type_ = NONE : géré par checkType
      checkType<T>();
      return static_cast<std::vector<T>*>(handler.get());
    }

    // une fonction push_back 
    // on fait un check_type, on cast et on appelle psuh_back
    template <typename T>
    void push_back(T x) {
      checkType<T>();
      ((std::vector<T> *) handler.get())->push_back(x);
    }

    // une fonction at
    // même démarche que push_back
    template <typename T>
    T at(size_t i) const {
      checkType<T>();
      return ((std::vector<T> *) handler.get())->at(i);
    }

    // une fonction push_back qui va convertir une chaine de caractères
    // dans le type cible
    void push_back_convert(const std::string & x) {
      switch(type_) {
        case INT : { 
           ((std::vector<int> *) handler.get())->push_back(std::stoi(x));
           break;
        }
        case FLOAT : {
          ((std::vector<float> *) handler.get())->push_back(std::stof(x));
          break;
        } 
        case DOUBLE : {
          ((std::vector<double> *) handler.get())->push_back(std::stod(x));
          break;
        }
        case STRING : {
          ((std::vector<std::string> *) handler.get())->push_back(x);
          break;
        }
        default :
          throw std::runtime_error("Can't convert string to type NONE");
      }
    }

    // la même qui prend un char *
    void push_back_convert(const char * x) {
      switch(type_) {
        case INT : { 
           ((std::vector<int> *) handler.get())->push_back(atoi(x));
           break;
        }
        case FLOAT : {
          char * end = NULL;
          ((std::vector<float> *) handler.get())->push_back(strtof(x, &end));
          break;
        } 
        case DOUBLE : {
          char * end = NULL;
          ((std::vector<double> *) handler.get())->push_back(strtod(x, &end));
          break;
        }
        case STRING : {
          std::string s(x);
          ((std::vector<std::string> *) handler.get())->push_back(s);
          break;
        }
        default :
          throw std::runtime_error("Can't convert char * to type NONE");
      }
    }


    // et la même qui prend un char * et s'arrête au premier blanc ; renvoie 
    // le pointeur vers le premier caractère non lu. Ceci pour pouvoir lire une 
    // ligne morceau par morceau (utile lors de la lecture de "white delimited values"
    // dans une DataStruct)
    char * push_back_token(const char * x) {
      char * end = NULL;
      switch(type_) {
        case INT : { 
           ((std::vector<int> *) handler.get())->push_back(strtoi(x, &end));
           break;
        }
        case FLOAT : {
          ((std::vector<float> *) handler.get())->push_back(strtof(x, &end));
          break;
        } 
        case DOUBLE : {
          ((std::vector<double> *) handler.get())->push_back(strtod(x, &end));
          break;
        }
        case STRING : {
          ((std::vector<std::string> *) handler.get())->push_back(strtos(x, &end));
          break;
        }
        default :
          throw std::runtime_error("Can't convert char * to type NONE");
      }
      return end;
    }

    // ~Column() {
    //     std::cout << "D° called for a Column of type " << type << " (0=INT,DOUBLE,FLOAT,3=STRING), before destr°, "<< handler.use_count() <<" ref to underliying data\n"; 
    // }

};

#endif
