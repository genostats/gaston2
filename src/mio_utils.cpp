/**
 * @file mio_utils.cpp
 * @brief file grouping functions not used yet, but useful for later
 * @date 2023-12-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "mio.hpp"
#include <iostream>
#include <fstream>

/**
 * @brief function opening a fstream with path and going  to_add bytes past the end of file, sizing it up
 * 
 * This function uses seekp with an offset from the end of the file, then puts a \0 byte 
 * (thus adding +1 byte in size hence the to_add-1) before flushing and closing. 
 * 
 * @param path 
 * @param to_add 
 * @return int 
 */
int resizing_file(std::string path, int to_add) {
    try {
        std::fstream file(path);
        if (!file ||!file.is_open()) throw std::runtime_error("Failed to open the file.");
        file.seekp(to_add - 1, std::ios_base::end);
        file.put('\0');
        file.flush();
        file.close(); // redundant I guess, flush done when closing
        return 0;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error when resizing file: " << e.what() << std::endl;
        return -1;
    }
}


int resizing_with_mio(std::string path){
    try {
        std::error_code error;
         mio::mmap_sink changing_file = mio::make_mmap_sink(path, 0, mio::map_entire_file, error);
        if (error) 
        { 
            std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
            throw std::runtime_error(errMsg);
        }
        // resizing file with streams without having to close or remap mio
        int to_add = 0; // TODO : to change or to take from input
        resizing_file(path, 0);
        // now writing past previous eof with mio shouldn't cause problems
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}