#ifndef MIO_UTILS
#define MIO_UTILS

#include <fstream>

int resizing_file(std::string path, int to_add);
int resizing_stream(std::ofstream & file, int to_add);

#endif // MIO_UTILS