#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory> // for shared_ptr

#ifndef _Column_
#define _Column_

enum datatype { INT, DOUBLE, FLOAT, STRING };
// TODO : how can I add a custom type for the user ? 

struct Column {
    datatype type_;
    std::shared_ptr<void> handler = nullptr;

  private:
    template <typename T>
    void setType();

  public:
    template <typename T>
    Column(std::vector<T> vec) : handler( std::make_shared<std::vector<T>>(vec) ) {
        setType<T>();
    }

    datatype type() const {
      return type_;
    }

    template <typename T>
    std::vector<T> * get() const;

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
    }

    // ceci fait le boulot pour le constructeur ci dessous
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

// inline should be used for all specializations made in header cos
// "function included in multiple source files must be inline" https://en.cppreference.com/w/cpp/language/inline

//specializations for setType
template <>
inline void Column::setType<int>() {
    type_ = INT;
}

template<>
inline void Column::setType<float>() {
    type_ = FLOAT;
}

template<>
inline void Column::setType<double>() {
    type_ = DOUBLE;
}

template<>
inline void Column::setType<std::string>() {
    type_ = STRING;
}


//specializations for get
template<>
inline std::vector<int> * Column::get() const {
  if (!handler) throw std::runtime_error("Handler was never instanciated");
  if (INT != type_) throw std::runtime_error("Trying to access with wrong type !");
  return static_cast<std::vector<int>*>(handler.get());
}

template<>
inline std::vector<float> * Column::get() const {
    if (!handler) throw std::runtime_error("Handler was never instanciated");
    if (FLOAT != type_) throw std::runtime_error("Trying to access with wrong type !");
    return static_cast<std::vector<float>*>(handler.get());
}

template<>
inline std::vector<double> * Column::get() const {
    if (!handler) throw std::runtime_error("Handler was never instanciated");
    if (DOUBLE != type_) throw std::runtime_error("Trying to access with wrong type !");
    return static_cast<std::vector<double>*>(handler.get());
}

template<>
inline std::vector<std::string> * Column::get() const {
    if (!handler) throw std::runtime_error("Handler was never instanciated");
    if (STRING != type_) throw std::runtime_error("Trying to access with wrong type !");
    return static_cast<std::vector<std::string>*>(handler.get());
}

#endif
