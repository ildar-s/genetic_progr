#include "ttk.h"

#include "ttree.h"
#include "brackets.h"

#include <vector>
#include "common_0.h"


const std::string var_symbol{"x"};

Ttk::Ttk(const std::string &s, const Tf *ff): 
    src(s),ff(ff)
{

    const Tf &ff_ = *this->ff;

    val="unknown";

    std::map<size_t,size_t> m = brackets::match_br(src);

    if(m.size())
    if((m.begin()->first==0) && (m.begin()->second==src.size()-1)){
        src.erase(static_cast<size_t>(0),static_cast<size_t>(1));
        src.pop_back();
        m = brackets::match_br(src);
    }

    std::map<size_t,size_t> m_ = brackets::top_level_br(m);


    std::string s_tmp(src);

    for(const auto &t: m_){
        size_t idx = t.first+1;
        size_t n = t.second-t.first-1;
        s_tmp.replace(idx,n,n,'(');
    }


    for(auto it=ff_.begin();it!=ff_.end();++it){

        if((*it)->isprefix())
            continue;

        std::string sep_ = (*it)->get_value_decor();
        std::string::size_type idx = s_tmp.find(sep_);
        if(idx!=std::string::npos){
            m_ch = brackets::split_m(s_tmp,sep_);

            if(static_cast<int>(m_ch.size())!=(*it)->get_num_args()){
                m_ch.clear();
                continue;
            }

            val = sep_;
            func_ind = std::distance(ff_.begin(),it);
            break;
        }
    }



    type=1; // type==1 unless terminal detected below



    if(!m_ch.size()){
        if(m_.size()==1){
            //prefix function


            for(auto it=ff_.begin();it!=ff_.end();++it){

                if(!((*it)->isprefix()))
                    continue;

                std::string sep_ = (*it)->get_value_decor();
                std::string::size_type idx = s_tmp.find(sep_);
                if(idx!=std::string::npos){
                    m_ch = brackets::split_m(s_tmp,sep_);

                    if(m_ch.size()!=1){
                        m_ch.clear();
                        continue;
                    }

                    val = sep_;
                    func_ind = std::distance(ff_.begin(),it);
                    break;
                }
            }
            
            // now checking if a list of args was saved as a single ch.
            if(m_ch.size()){ //prefix f match found in the above loop
                if(ff_[func_ind]->get_num_args()>1){

                    auto t = *m_ch.begin();
                    size_t idx = t.first+1;
                    size_t n = t.second-t.first-1;

                    std::string src_tmp = src.substr(idx,n);
                    Ttk tk_tmp(src_tmp,&Ttree::funcs_l);

                    m_ch.clear();

                    for(const auto &t : tk_tmp.m_ch){
                        m_ch[t.first+idx]=t.second+idx;
                    }
                    
                }
            }
        }
    }



    int ind;
    std::string sind;
    auto is_var = [&ind,&sind](const std::string &s)->bool{
        typedef std::string::size_type ssize;
        ssize idx;

        idx=s.find(var_symbol+"[");
        if(idx==std::string::npos)
            return false;

        idx=s.find_first_not_of(var_symbol+"[",idx);
        if(idx==std::string::npos)
            return false;

        try{
            ind = std::stoi(s.substr(idx));
        }
        catch(...){
            ind = -1;
            ssize idx_=s.find("]",idx);
            if(idx_==std::string::npos)
                return false;

            ssize i0 = idx+var_symbol.size()+1;
            ssize n = idx_-i0;
            sind=s.substr(i0,n);
        }

        return true;
    };



    if(!m_ch.size()){
        type=0; // terminal
        val = src;

        if(is_var(val)){
            isvar = true;
            var_ind = ind;
            if(var_ind<0)
                var_sind = sind;
        }
        else{
            isvar = false;
            try{
                cvalue = std::stod(val);
            }
            catch(...){
                val += "_ERR";
            }
        }
    }


}

std::unique_ptr<Ttk> Ttk::get_ch(const size_t i)const{

    if(i>m_ch.size()-1)
        return nullptr;

    auto it = std::next(m_ch.begin(),i);
    size_t idx = it->first;
    size_t n = it->second-it->first+1;

    std::string s_ = src.substr(idx,n);

    return std::make_unique<Ttk>(s_,ff);
}




