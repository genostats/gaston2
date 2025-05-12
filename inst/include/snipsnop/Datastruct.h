#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory> // for shared_ptr

#ifndef _DataStruct_
#define _DataStruct_

enum datatype { INT, DOUBLE, FLOAT, STRING };
// TODO : how can I add a custom type for the user ? 

struct Column {
    datatype type;
    std::shared_ptr<void> handler = nullptr;

  private:
    template <typename T>
    void setType();

  public:
    template <typename T>
    Column(std::vector<T> vec) : handler( std::make_shared<std::vector<T>>(vec) ) {
        setType<T>();
    }

    template <typename T>
    std::vector<T> * get();

    // ~Column() {
    //     std::cout << "D° called for a Column of type " << type << " (0=INT,DOUBLE,FLOAT,3=STRING), before destr°, "<< handler.use_count() <<" ref to underliying data\n"; 
    // }

};

// inline should be used for all specializations made in header cos
// "function included in multiple source files must be inline" https://en.cppreference.com/w/cpp/language/inline

//specializations for setType
template <>
inline void Column::setType<int>() {
    type = INT;
}

template<>
inline void Column::setType<float>() {
    type = FLOAT;
}

template<>
inline void Column::setType<double>() {
    type = DOUBLE;
}

template<>
inline void Column::setType<std::string>() {
    type = STRING;
}


//specializations for get
template<>
inline std::vector<int> * Column::get() {
  if (!handler) throw std::runtime_error("Handler was never instanciated");
  if (INT != type) throw std::runtime_error("Trying to access with wrong type !");
  return static_cast<std::vector<int>*>(handler.get());
}

template<>
inline std::vector<float> * Column::get() {
    if (!handler) throw std::runtime_error("Handler was never instanciated");
    if (FLOAT != type) throw std::runtime_error("Trying to access with wrong type !");
    return static_cast<std::vector<float>*>(handler.get());
}

template<>
inline std::vector<double> * Column::get() {
    if (!handler) throw std::runtime_error("Handler was never instanciated");
    if (DOUBLE != type) throw std::runtime_error("Trying to access with wrong type !");
    return static_cast<std::vector<double>*>(handler.get());
}

template<>
inline std::vector<std::string> * Column::get() {
    if (!handler) throw std::runtime_error("Handler was never instanciated");
    if (STRING != type) throw std::runtime_error("Trying to access with wrong type !");
    return static_cast<std::vector<std::string>*>(handler.get());
}

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
};

#endif
