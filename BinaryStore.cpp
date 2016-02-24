#include "BinaryStore.hpp"

binary_store::binary_store(const std::string file_name,int given_element_size, bool flag):element_size(given_element_size)
{
    open(file_name,flag);
}
void binary_store::open(const std::string file_name,bool flag)
{
    const char * f_name = file_name.c_str();
    file_ptr = fopen(f_name,flag?"wb":"rb");
    if( file_ptr  == NULL )
    {
        std::cerr << "error during fopen" << std::endl;
        exit(1);
    }
}
void binary_store::write_element(void *element)
{

    if (!fwrite(element, element_size, 1, file_ptr))
    {
        std::cerr <<"error: can't fwrite (disk full?)" << std::endl;
        exit(1);
    }

}
size_t binary_store::read_element( void *element)
{
    return fread(element, element_size,1, file_ptr);
}

size_t binary_store::read_element_buffer(void *element,size_t nb_elements)
{
    return fread(element, element_size,nb_elements, file_ptr);
}

void binary_store::write( void *element, int given_element_size)
{
    if (!fwrite(element, given_element_size, 1, file_ptr))
    {
        std::cerr <<"error: can't fwrite (disk full?)" << std::endl;
        exit(1);
    }
}

size_t binary_store::read( void *element, int given_element_size)
{
    return fread(element, given_element_size,1, file_ptr);
}

void binary_store::rewind_all()
{
    rewind(file_ptr);
}

void binary_store::close()
{
        fclose(file_ptr);
}
