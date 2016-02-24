#include "Utility.hpp"
std::string prefix;
time_t rawtime ;
std::string return_file_name(const std::string file_name)
{
    std::string result_file_name;
    if(prefix.length()>0)
        result_file_name=prefix+"."+file_name;
    else
        result_file_name=file_name;
    return   result_file_name;

}
void delete_created_file(const std::string file_name)
{
     
    if( remove( file_name.c_str() ) != 0 )
     std::cerr<<"--- error deleting the file. "<<file_name<<std::endl;
}
