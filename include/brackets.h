#ifndef BRACKETS_H_
#define BRACKETS_H_

#include <stack>

namespace brackets{
    inline std::map<size_t,size_t> match_br(const std::string &s){

        std::map<size_t,size_t> m;
        std::stack<size_t> st;

        for(size_t i=0;i<s.size();++i){
            if(s[i]=='('){
                st.push(i);
            }
            if(s[i]==')'){
                if(!st.size())
                    throw xX_0("reduce_br: wrong closing bracket");
                m[st.top()]=i;
                st.pop();
            }
        }
        if(st.size())
            throw xX_0("reduce_br: wrong opening bracket");

        return m;
    }

    inline std::string reduce_br(const std::string &s){

        std::string s_ = rm_ws(s);

        std::map<size_t,size_t> m=match_br(s_);
        if(!m.size()) return s_;

        auto it0 = m.begin();
        auto it1 = m.begin();
        it1++;

        for(;it1!=m.end();++it0,++it1){
            int open0 = it0->first;
            int close0 = it0->second;
            int open1 = it1->first;
            int close1 = it1->second;
            if((open1-open0==1) && (close1-close0==-1)){
                s_[it0->first]=' ';
                s_[it0->second]=' ';
            }
        }

        s_ = rm_ws(s_);
        return s_;
    }


    inline std::map<size_t,size_t> split_m(const std::string &s,const std::string &sep){

        std::map<size_t,size_t> m;
        if(!s.size())
            return m;

        size_t start = 0;

        while(start<s.size()){

            size_t idx = s.find(sep,start);

            if(idx==std::string::npos){
                m[start]=s.size()-1;
                break;
            }

            if(idx>start){
                m[start]=idx-1;
            }
            start = idx+sep.size();
        }

        return m;
    }


    inline std::map<size_t,size_t> top_level_br(const std::map<size_t,size_t> &br){

        std::map<size_t,size_t> m(br);
        std::map<size_t,size_t> m_;

        while(m.size()){
            auto it = m.begin();
            auto i0 = *it;
            it++;

            std::unordered_set<size_t> idel;

            for(;it!=m.end();++it){
                if(it->first > i0.first && it->second < i0.second)
                    idel.insert(it->first);
            }
            for(const auto i: idel){
                m.erase(i);
            }

            it = m.begin();
            std::copy(it,std::next(it),inserter(m_,m_.begin()));
            m.erase(it);
        }

        return m_;
    }
}


#endif /* BRACKETS_H_ */
