#ifndef BINARYSTORE_H
#define BINARYSTORE_H
#include <string>
#include <stdio.h>
#include <iostream>
#include <algorithm>

class binary_store
{
public:
    binary_store(const std::string file_name,int given_element_size, bool flag);
    void write_element(void *element);
    size_t read_element(void *element);
    void write( void *element, int given_element_size);
    size_t read( void *element, int given_element_size);
    size_t read_element_buffer(void *element,size_t nb_elements);
    void rewind_all();
    void close();
    void open(const std::string file_name,bool write);
private:
    FILE * file_ptr;
    const int element_size;
};


#endif // BINARYSTORE_H
