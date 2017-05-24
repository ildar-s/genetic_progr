#ifndef COMMON_0_H_
#define COMMON_0_H_

#include <iostream>
#include <random>
#include <cstring>
#include <cassert>
#include <algorithm>



typedef std::vector<int> tvi;
typedef std::vector<size_t> tvsz;
typedef std::vector<double> tvd;
typedef std::vector<tvd> tvvd;

typedef std::pair<tvvd,tvd> txy;

enum class LOG : char {SILENT,INFO,DEBUG};

inline void pr(){
    std::cout << std::endl;
}
template <typename T, typename... Types>
void pr(const T& firstArg, const Types&... args){
    std::cout << firstArg <<" ";
    pr(args...);
}
inline void pr_(){
    std::cout << std::endl;
}
template <typename T, typename... Types>
void pr_(const T& firstArg, const Types&... args){
    std::cout << firstArg;
    pr_(args...);
}
class xX_0{
    public:
        xX_0(const std::string &msg){ pr("ERROR: ",msg);}
};




inline int get_discrete_s(const int beg, const int end, const double lam, 
        std::default_random_engine &dre){
    /* [beg,end] */

    if(end<beg)
        throw xX_0("get_discrete_s: beg="+std::to_string(beg)+" end="+std::to_string(end));

    int n = end+1;
    std::vector<double> pp(n,0);

    for(int i=0;i<end-beg+1;i++){
        double p0 = std::exp(-lam*i);
        pp[beg+i] = p0;
    }
    std::discrete_distribution<int> dd(pp.begin(),pp.end());
    int ret = dd(dre);
    if(ret<beg)
        throw xX_0("get_discrete_s: ret<beg beg="+std::to_string(beg)+" ret="+std::to_string(ret));
    return ret;
}

inline char* trim(char *p){
    while(std::isspace(*p))p++;
    if(*p==0)return p;

    char *e = p+std::strlen(p)-1;
    
    while(std::isspace(*e))e--;
    *(e+1)=0;
    return p;
}



const char ws_chars[] = " \f\n\r\t\v";

inline std::string trim(const std::string &s){
    auto i_beg = s.find_first_not_of(ws_chars);
    if(i_beg==std::string::npos)
        return "";
    auto i_end = s.find_last_not_of(ws_chars);
    return s.substr(i_beg,i_end-i_beg+1);
}

inline std::string rm_ws(const std::string &s){
    static auto isws = [](const char &c){ return !std::isspace(c);    };
    std::string ret;
    std::copy_if(s.begin(),s.end(),std::back_inserter(ret),isws);
    return ret;
}


#endif /* COMMON_0_H_ */
