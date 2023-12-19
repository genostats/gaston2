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
        std::fstream fifi(path, std::ios_base::out);
        if (!fifi ||!fifi.is_open()) throw std::runtime_error("Failed to open the file;");
        // TODO : think if useful without mio ? Should I not use mio here ? If yes, passed it in args ?
        std::error_code error;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error in resizing file: " << e.what() << std::endl;
    }
}


int test_dynamic()
{
    try {
        std::string path = "smol.bin";
        /*
        std::ofstream file(path, std::ios::binary);
        if (!file ||!file.is_open())
        {
            throw std::runtime_error("Failed to open the file;");
        }
        file.seekp(9);
        file.put('\0');
        file.close();
        */

        std::fstream fifi(path, std::ios_base::out);
        if (!fifi ||!fifi.is_open())
        {
            throw std::runtime_error("Failed to open the file;");
        }
        fifi.seekp(9);
        fifi.put('\0');
        fifi.flush();

        std::error_code error;
        // mio::mmap_sink mio::make_mmap_sink<std::string>(const std::string &token, size_t offset, size_t length, std::error_code &error)
        mio::mmap_sink changing_file = mio::make_mmap_sink(path, 0, 10, error);
        if (error) 
        { 
            std::string errMsg = "Error code " + std::to_string(error.value()) + ", Failed to map the file : " + error.message();
            throw std::runtime_error(errMsg);
        }

        //returns the actual number of bytes that were mapped
        std::cout << "Matrix 8 size mapped: " << changing_file.mapped_length() << std::endl;
        // doesn't seem to stop even if I'm after mapped region ?
        for (size_t j = 0; j < 10; j++)
        {
            changing_file.data()[j] = 65 + j%26;
        }

        // changing_file.sync(error);
        // changing_file.unmap();
        // on essaie de le rallonger


        std::cout << fifi.eof() << '\n'; // if 1 then true
        // std::cout << fifi.eof() << '\n'; // if 1 then true
        if (fifi.good()) std::cout << "No problems here\n"; 
        //in fstreams get and put pointers are always the same
        // fifi.seekp(0, std::ios_base::end);
        std::cout << fifi.tellp() << "\n";
        // std::cout << fifi.eof() << '\n'; // if 1 then true
        std::cout << fifi.eof() << '\n'; // if 1 then true
        std::cout << "get = " << fifi.get() << "\n";
        std::cout << fifi.eof() << '\n'; // if 1 then true       
        // prolonger de 10
        // fifi.clear(); // clear error flags before seeking pos past eof !
        
        fifi.seekp(9, std::ios_base::end);
        std::cout << fifi.eof() << '\n'; // if 1 then true
        std::cout << "get = " << fifi.get() << "\n";
        std::cout << fifi.eof() << '\n'; // if 1 then true
        fifi.put('\0');
        fifi.flush();
        

        for (size_t j = 10; j < 20; j++)
        {
            changing_file.data()[j] = 97 + j%26;
        }
        fifi.close();
        /*
        // trying to remap with a different offset, further ? 
        // void map<String>(const String &path, const size_t offset, const size_t length, std::error_code &error)
        changing_file.map(path, 20, 10, error);
        std::cout << "Matrix 8 size mapped: " << changing_file.mapped_length() << std::endl;
        for (size_t j = 0; j < 10; j++)
        {
            changing_file.data()[j] = 41;
        }

        // changing_file.data()[2049] = 65;

        changing_file.sync(error);
        changing_file.unmap();
        //std::ofstream file(path, std::ifstream::binary);
        //file.seekp(2049);
        */
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}