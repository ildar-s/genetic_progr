#include "tpd.h"
#include "common_0.h"

#include <sstream>
#include <fstream>

void Tpd::skip_ws(std::ifstream &file){
    char c;
    while(file.get(c)){
        if(!std::isspace(c)){
            file.unget();
            break;
        }
    }
}


std::vector<std::string> Tpd::fill_(std::ifstream &file, int &nc){
    skip_ws(file);
    std::stringstream ss;

    nc=0;
    char c;
    while(file.get(c)){
        if(c=='\n' || c=='\r')
            break;
        ss.put(c);
        nc++;
    }

    const int N=100;
    char buf[N];
    std::vector<std::string> v;

    const char sep=',';

    if(ss.peek()==sep){
        ss.get(buf,1+1);
    }

    while(ss.get(buf,N,sep)){
        v.push_back(trim(buf));
        ss.get(buf,1+1);
    }

    return v;
}


void Tpd::read_csv(const std::string &filename,
        const bool has_header){

    std::ifstream file(filename,std::ios::in);
    if(!file){
        pr(std::system("ls -l >_tmp_output_of_ls"));
        std::cout << std::ifstream("_tmp_output_of_ls").rdbuf();
        throw xX_0("reading \""+filename+"\" failed");
    }
    

    int nc;
    columns = fill_(file,nc);
    if(!columns.size())
        throw xX_0("\"columns\" is empty");



    std::vector<std::string> first_line;

    if(has_header){
        int nc0;
        first_line = fill_(file,nc0);
        if(nc0>nc)nc=nc0;
    }
    else{
        first_line = columns;
    }


    size_t ncols = first_line.size();
    if(!ncols)
        throw xX_0("\"first_line\" is empty");

    if(ncols!=columns.size()){
        if((int)ncols-(int)columns.size()==1){
            columns.insert(columns.begin(),"empty");
        }
        else{
            throw xX_0("read_csv: wrong num of columns");
        }
    }


    std::vector<double> data_l;

    for(size_t i=0;i<first_line.size();++i){
        double x0 = std::stod(first_line[i]);
        data_l.push_back(x0);
    }
    data.push_back(data_l);


    std::stringstream ss;

    while(file.get(*ss.rdbuf())){

        std::string s = rm_ws(ss.str());
        ss.str(""); 
        ss.clear();
        if(!s.size())
            continue;

        data_l.clear();
        std::size_t i0{0}, idx;

        for(size_t i=0;i<ncols;++i){
            double x0;
            try{
                x0 = std::stod(s.c_str()+i0,&idx);
            }
            catch(...){
                throw xX_0("read_csv: at line \""+s+"\"");
            }
            data_l.push_back(x0);

            i0 += idx+1;

            if(i0>=s.size())
                break;
        }

        data.push_back(data_l);
        skip_ws(file);
    }
}





