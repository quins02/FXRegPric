#ifndef DOWNLOAD_H
#define DOWNLOAD_H

#include <string>
#include <vector>

size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream);
void Data_Download(std::string Q_loc, std::string API_code);
std::vector <double> TS_Vec(std::string Name);

#endif