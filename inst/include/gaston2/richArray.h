#include<stdexcept>

// this class allows to declare arrays that can be used with a part of the vector interface
// they have members size(), begin(), end(), at()
// (but no push_back etc)
template<size_t size_, typename T>
class richArray {
  T data[size_];
  public:
  inline size_t size() { return size_; }
  inline T & operator[](size_t i) { return data[i]; }
  inline T & at(size_t i) {
    if(i < size_)
      return data[i];
    else 
      throw std::out_of_range("richArray range error");
  }
  inline T * begin() { return data; }
  inline T * end() { return data + size_; }
};

