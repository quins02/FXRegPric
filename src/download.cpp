#include <stdio.h>
#include <curl/curl.h>
#include <string>
#include <sstream>

#include <vector>
#include <iostream>
#include <fstream>

size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t written = fwrite(ptr, size, nmemb, stream);
    return written;
}

void Data_Download(std::string Q_loc, std::string API_code){
    std::string loc = "https://www.quandl.com/api/v3/datasets/" + Q_loc + ".csv?sort_order=dsc&auth_token=" + API_code;
    CURL *curl;
    FILE *fp;
    CURLcode res;
    //char *url = "https://www.quandl.com/api/v1/datasets/ECB/EURUSD.csv?sort_order=dsc&auth_token=8dM9KB3tW11HYvxXGCBK";
    char *url = &loc[0u];
    char outfilename[FILENAME_MAX] = "tmp/Out.csv";
    curl = curl_easy_init();
    if (curl) {
        fp = fopen(outfilename,"wb");
        curl_easy_setopt(curl, CURLOPT_URL, url);
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
        res = curl_easy_perform(curl);
        /* always cleanup */
        curl_easy_cleanup(curl);
        fclose(fp);
    }
}

std::vector <double> TS_Vec(std::string Name){
    std::vector <double> Vec;
    std::string line;
    std::string Str;
    std::string tmp;
    double V;
    int i =0;   

    std::ifstream data (Name.c_str());

    getline(data, tmp);
    
    while(getline(data,line)){
        std::istringstream Lin(line);
        while(getline(Lin,Str,',')){
            if(i==1){
                std::istringstream tmp(Str);
                tmp >> V;
                Vec.push_back(V);
                //std::cout<<V<<std::endl;
            }
            i++;
        }
        i=0;
    }   

    return Vec;
}