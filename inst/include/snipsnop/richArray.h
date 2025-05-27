#include<stdexcept>

// this class allows to declare arrays that can be used with a part of the vector interface
// they have members size(), begin(), end(), at()
// (but no push_back etc)
template<size_t size_, typename T>
class richArray {
  T data[size_];
  public:
  size_t size() { return size_; }
  T & operator[](size_t i) { return data[i]; }
  T & at(size_t i) {
    if(i < size_)
      return data[i];
    else 
      throw std::out_of_range("richArray range error");
  }
  T * begin() { return data; }
  T * end() { return data + size_; }
};

