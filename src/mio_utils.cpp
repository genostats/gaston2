/**
 * @file mio_utils.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2023-12-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "mio.hpp"
#include <iostream>
#include <fstream>


int resizing_file(std::string path, int to_add) {
    try {
        std::fstream file(path, std::ios_base::out);
        if (!file ||!file.is_open()) throw std::runtime_error("Failed to open the file.");
        // TODO : think if useful without mio ? Should I not use mio here ? If yes, passed it in args ?
        file.seekp(to_add, std::ios_base::end);
        file.put('\0');
        file.flush();
        file.close(); // redundant I guess, flush done when closing
        return 0;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error when resizing file: " << e.what() << std::endl;
        return -1;
    }
}