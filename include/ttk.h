#ifndef TTK_H_
#define TTK_H_

#include "tfunc.h"
#include <string>
#include <map>
#include <memory>


class Ttk{
    public:
        typedef std::vector<std::shared_ptr<Tfunc_base>> Tf;

        Ttk(const std::string &s, const Tf *ff);
        int num_of_ch()const{ return m_ch.size(); }
        std::unique_ptr<Ttk> get_ch(const size_t i)const;

        std::string val;

        int type;

        int func_ind;

        bool isvar;
        int var_ind;
        std::string var_sind;
        double cvalue;

    private:
        std::map<size_t,size_t> m_ch;
        std::string src;
        const Tf *ff;


        /*
         

           if(type==1){ 
                it's function
                func_ind =  index of function in ff
                val = function name not used further because the above func_ind identifies the function
           }
           else{
                it's terminal
                if(isvar){
                    it's variable
                    var_ind = numerical index recovered from []
                    var_sind = text index from []
                    only one of the above inds are defined
                }
                else{
                    it's constant
                    cvalue = it's value
                }
           }



         */




};


#endif //TTK_H_
