#include "tel.h"
#include "ttree.h"
#include "common_ttrees.h"
#include "common_0.h"



Tel::Tel(const int depth, const std::weak_ptr<Tel> &parent, std::unique_ptr<Ttk> tk_):
    depth(depth),type(tk_->type),parent(parent),tk(std::move(tk_))
{}

Tel::Tel(const int depth,const std::weak_ptr<Tel>& parent,
        Ttree *ptree){

    this->depth = depth;
    this->parent = parent;

    if(
            ptree->urd(dre) > ptree->p.element_function_probability 
            || 
            depth>=ptree->p.max_depth

            ) {

        type=0; // terminal
    }
    else{
        type=1; // function 
    }


    if(type==0){ //terminal
        if( ptree->urd(dre) < ptree->p.terminal_const_probability
               || ptree->p.ndim <= 0 ){
            // const
            value_ind = (*ptree->p_uid_term_c)(dre);
        }
        else{
            //var
            value_ind = (*ptree->p_uid_term_v)(dre);
        }

        term = ptree->terms[value_ind];
    }
    else{ //function
        value_ind = (*ptree->p_uid_func)(dre);
        func = ptree->funcs[value_ind].get();
    }
}
Tel::Tel(const int depth,const std::weak_ptr<Tel>& parent, 
        const Tel& src, const Ttree *ptree){

    this->depth = depth;
    this->parent = parent;
    type = src.type;

    value_ind = src.value_ind;

    if(type==0){ // terminal
        term = ptree->terms[value_ind];
    }
    else{ //function
        func = ptree->funcs[value_ind].get();
    }
}

std::string Tel::get_value_decor() const {
    if(type==0)
        return term.second;
    else
        return func->get_value_decor();
}




