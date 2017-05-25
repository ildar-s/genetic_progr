#include "ttree.h"
#include "common_0.h"

#include <map>
#include <set>
#include <stack>
#include <algorithm>
//#include <thread>


//#include <chrono>
//using namespace std::chrono;


std::vector<std::shared_ptr<Tfunc_base>> Ttree::funcs;
std::vector<std::shared_ptr<Tfunc_base>> Ttree::funcs_l;
bool Ttree::static_up_to_date = false;

//std::uniform_real_distribution<double> Ttree::urd(0,1);
//std::unique_ptr<std::uniform_int_distribution<int>> Ttree::p_uid_term;
//std::unique_ptr<std::uniform_int_distribution<int>> Ttree::p_uid_term_c;
//std::unique_ptr<std::uniform_int_distribution<int>> Ttree::p_uid_term_v;
//std::unique_ptr<std::uniform_int_distribution<int>> Ttree::p_uid_func;
//int Ttree::ndim_st = -1;



double Ttree::J_default(const std::vector<double> &y0,const std::vector<double> &y1){
    assert(y0.size()==y1.size());
    double s=0;
    for(size_t i=0;i<y0.size();++i){
        double z_ = y0.at(i)-y1.at(i);
        s += z_*z_;
    }
    return s;
}

// y_true must be {-1,1}
double Ttree::J_logit(const std::vector<double> &y_true,const std::vector<double> &y_pred){
    assert(y_true.size()==y_pred.size());

    double s=0;
    for(size_t i=0;i<y_true.size();++i){

        double z_ = y_true.at(i) * y_pred.at(i);

        s += std::log(1+std::exp(-z_));
    }
    return s;
}

void Ttree::init_static(){
    funcs.clear();
    funcs.push_back(std::make_shared<Tfunc_plus>());
    funcs.push_back(std::make_shared<Tfunc_minus>());
    //funcs.push_back(std::make_shared<Tfunc_pref_minus>());
    funcs.push_back(std::make_shared<Tfunc_mult>());
    funcs.push_back(std::make_shared<Tfunc_div>());
    funcs.push_back(std::make_shared<Tfunc_log>());
    funcs.push_back(std::make_shared<Tfunc_sin>());
    funcs.push_back(std::make_shared<Tfunc_cos>());
    funcs.push_back(std::make_shared<Tfunc_tanh>());
    funcs.push_back(std::make_shared<Tfunc_pow2>());
    funcs.push_back(std::make_shared<Tfunc_pow3>());
    //funcs.push_back(std::make_shared<Tfunc_min>());
    //funcs.push_back(std::make_shared<Tfunc_max>());

    funcs_l.clear();
    funcs_l.push_back(std::make_shared<Tfunc_list2>());



    static_up_to_date=true;
}

void Ttree::update_ndim(const Ttree_parameters &p){

    if(!static_up_to_date)
        init_static();

    this->p = p;

    consts.clear();
    double c_range = p.consts_max-p.consts_min;
    for(int i=0;i<p.consts_n;++i){
        double v = i*1./(p.consts_n-1)*c_range-p.consts_min;
        consts.push_back(v);
    }

    terms.clear();

    char buf[100];
    for(auto &c : consts){
        sprintf(buf, "%6.2f", c);
        std::string cs(trim(buf));

        if(c<0)
            cs="("+cs+")";
        
        terms.push_back(std::make_pair(&c,cs));
    }

    x.resize(p.ndim);
    for(size_t i=0;i<x.size();++i){
        terms.push_back(std::make_pair(&x[i],"x["+std::to_string(i)+"]"));
    }

    p_uid_term_c.reset(new std::uniform_int_distribution<int>(0,consts.size()-1));
    p_uid_term_v.reset(new std::uniform_int_distribution<int>(consts.size(),consts.size()+p.ndim-1));
    p_uid_term.reset(new std::uniform_int_distribution<int>(0,consts.size()+p.ndim-1));
    p_uid_func.reset(new std::uniform_int_distribution<int>(0,funcs.size()-1));
}

Ttree::Ttree(const Ttree_parameters * const p_):
    err(-1),use_count(0),urd(0,1),J_(J_default)
{
    if(p_!=nullptr){
        update_ndim(*p_);
    }
    else{
        Ttree_parameters p;
        update_ndim(p);
    }
}


Ttree::Ttree(const Ttree_parameters &p,
        const std::function<double(const std::vector<double>&,const std::vector<double>&)> &J):
    Ttree(&p)
{
    if(J!=nullptr){
        if(p.logit)
            this->J_ = J_logit;
        else
            this->J_ = J;
    }

    root = std::make_shared<Tel>(0,std::weak_ptr<Tel>(),this);

    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);

    for(auto& el : els){
        if(el->type==1){ // function

            for(int i=0;i<el->func->get_num_args();++i){
                std::shared_ptr<Tel> ch_  = std::make_shared<Tel>(el->depth+1,el,this);
                el->ch.push_back(ch_);
                els.push_back(ch_);
            }
            
        }
    }
}

Ttree::Ttree(const Ttree &src, const int cpoint):
    Ttree(&src.p)
{
    J_ = src.J_;

    root = std::make_shared<Tel>(0,std::weak_ptr<Tel>(),this);

    const std::shared_ptr<Tel> &root_src = src.root;
    std::shared_ptr<Tel> root_src_ = nullptr;


    if(cpoint<0){
        root_src_ = root_src;
    }
    else {
        std::list<std::shared_ptr<Tel>> els;

        els.push_back(root_src);
        int i00=-1;
        for(const auto& el : els){
            i00++;
            if(el==nullptr)
                continue;

            if(i00==cpoint){
            //if(el->i == cpoint){
                root_src_ = el;
                break;
            }

            std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
        }
    }
    if(root_src_ == nullptr)
        return;

    assign_els(root_src_,root);
}



Ttree::Ttree(const std::string &s, const Ttree_parameters * const p_):
    Ttree(p_)
{

    std::string s_ = reduce_br(s);
    //pr("\ninitial ->",s_);

    root = std::make_shared<Tel>(0,std::weak_ptr<Tel>(),std::make_unique<Ttk>(s_,&funcs));

    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);

    for(auto &el:els){
        if(el->type==1){ // function

            //pr_("  func=\"",el->tk->val,"\" ");

            size_t nch = el->tk->num_of_ch();

            for(size_t i=0;i<nch;++i){
                std::shared_ptr<Tel> ch_  = std::make_shared<Tel>(el->depth+1,el,el->tk->get_ch(i));

                el->ch.push_back(ch_);
                els.push_back(ch_);
            }
        }
        else{
            //pr_("\tterm=\"",el->tk->val,"\" ,is var=",el->tk->isvar);
        }
    }

    if(!els.size())
        throw xX_0("Ttree(s): no elements generated");

    if(p_==nullptr){
        p.max_depth = els.back()->depth;
    }


    // are var_inds represented by strings(e.g. x["column 3"]) ?
    bool is_sind = false;
    for(auto &el:els){
        if(el->type==0){
            if(el->tk->isvar){ //var
                if(el->tk->var_ind<0){
                    is_sind = true;
                    break;
                }
            }
        }
    }

    if(is_sind){
        std::unordered_set<std::string> sinds;

        for(auto &el:els){
            if(el->type==0){
                if(el->tk->isvar){ //var
                    sinds.insert(el->tk->var_sind);
                }
            }
        }
        for(auto &el:els){
            if(el->type==0){
                if(el->tk->isvar){ //var
                    auto pos = sinds.find(el->tk->var_sind);
                    if(pos==sinds.end())
                        throw xX_0("Ttree(): sinds values bug");
                    size_t idx = std::distance(sinds.begin(),pos);
                    el->tk->var_ind = idx;
                }
            }
        }

        if(p_==nullptr){
            p.ndim = sinds.size();
        }
    }
    else{
        int max_ind=-1;
        for(auto &el:els){
            if(el->type==0){
                if(el->tk->isvar){
                    if(el->tk->var_ind>max_ind){
                        max_ind = el->tk->var_ind;
                    }
                }
            }
        }

        if(p_==nullptr){
            p.ndim = max_ind+1;
        }
    }

    if(p.ndim>10000)
        throw xX_0("Ttree(s): ndim too high="+std::to_string(p.ndim));

    if(p_==nullptr){
        update_ndim(p);
    }

    for(auto &el:els){
        if(el->type==1){ // function
            el->value_ind = el->tk->func_ind;
            el->func = funcs[el->value_ind].get();
        }
        else{ //terminal
            if(el->tk->isvar){ //var
                el->value_ind = el->tk->var_ind+consts.size();
            }
            else{//const

                double cval = el->tk->cvalue;

                double min_diff=1e10;
                int idx=-1;

                for(size_t i=0;i<consts.size();++i){
                    double diff = std::abs(consts[i]-cval);
                    if(diff<min_diff){
                        min_diff = diff;
                        idx = i;
                    }
                }
                if(idx<0)
                    xX_0("error parsing const from string");
                el->value_ind = idx;
            }
            el->term = terms[el->value_ind];
        }
    }
}


void Ttree::assign_els(const std::shared_ptr<Tel> &root_src, std::shared_ptr<Tel> &root_dest){

    int depth = root_dest->depth;
    std::weak_ptr<Tel> parent = root_dest->parent;
    root_dest = std::make_shared<Tel>(depth,parent,*root_src,this);

    std::list<std::shared_ptr<Tel>> els_dest;
    els_dest.push_back(root_dest);

    std::list<std::shared_ptr<Tel>> els_src;
    els_src.push_back(root_src);

    auto el_src = els_src.begin();
    auto el_dest = els_dest.begin();

    for(; el_src!=els_src.end(); ++el_src,++el_dest){

        if(*el_src==nullptr) 
           continue;

        int depth = (*el_dest)->depth;

        for(auto& el_ch : (*el_src)->ch){
            std::shared_ptr<Tel> el_ch_dest = nullptr;
            if(el_ch!=nullptr)
                el_ch_dest = std::make_shared<Tel>(depth+1,*el_dest,*el_ch,this);
            (*el_dest)->ch.push_back(el_ch_dest);

            els_src.push_back(el_ch);
            els_dest.push_back(el_ch_dest);
        }
    }
}

typedef std::shared_ptr<Tel> spTel;

void Ttree::replace(const int cpoint, const Ttree &src){

    std::shared_ptr<Tel> *p_el_c = nullptr;
    
    std::list<std::reference_wrapper<spTel>> els;
    els.push_back(root);

    int i00 = -1;
    for(spTel& el : els){
        i00++;
        if(el==nullptr)
            continue;

        //if(el->i==cpoint){
        if(i00==cpoint){
            p_el_c = &el;
            break;
        }

        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
    }

    assign_els(src.root,*p_el_c);
}



double Ttree::eval(const std::vector<double> &x){
    this->x = x;
    
    std::list<std::shared_ptr<Tel>> els;
    if(root!=nullptr)
        els.push_back(root);
    for(auto& el : els){
        if(el==nullptr)
            continue;

        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
        el->type_tmp = el->type;
        if(el->type==0) //terminal
            el->value_tmp = *el->term.first;
    }


    int iii = 0;
    while(true){
        
        if(iii++>1000){
            pr("iii=",iii);
            throw xX_0("mf-error");
        }


        bool any_funcs_left = false;
        for( auto& el : els){
            if(el==nullptr)
                continue;

            if(el->type_tmp==1){
                any_funcs_left = true;

                bool all_terminals = true;
                for(auto& el_ : el->ch){
                    if(el_->type_tmp==1){ //function
                        all_terminals = false;
                        break;
                    }
                }

                if(all_terminals){
                    std::vector<double> args;
                    for(auto& el_ : el->ch){
                        args.push_back(el_->value_tmp);
                    }
                    if(!args.size()){
                       throw xX_0("mf-error");
                    }

                    el->value_tmp = el->func->f(args);
                    el->type_tmp = 0;
                }
            }

        }

        if(!any_funcs_left)
            break;
    }

    double ret_value = 0;
    if(root!=nullptr)
        ret_value = root->value_tmp;

    return ret_value;
}


int Ttree::count() const {
    int n=0;

    std::list<std::shared_ptr<Tel>> els;
    if(root != nullptr)
        els.push_back(root);
    for(const auto& el : els){
        if(el != nullptr)
            std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
        n++;
    }

    return n;
}
void Ttree::list_addr() const {

    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);

    pr("start of listing");
    for(auto& el : els){
        std::cout << el.get() << " ";
        if(el==nullptr)
            continue;
        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
    }
    pr("");
    pr("end of listing");
}

void Ttree::init_step0_w(){

    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);

    for(auto& el : els){
        if(el==nullptr)
            continue;

        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
        el->type_tmp = el->type;
        el->i_left=0;
        el->w=1;
    }

    int iii = 0;
    while(true){
        
        if(iii++>1000){
            pr("iii=",iii);
            throw xX_0("mf-error");
        }

        bool any_funcs_left = false;
        for( auto& el : els){
            if(el->type_tmp==1){
                any_funcs_left = true;

                bool all_terminals = true;
                for(auto& el_ : el->ch){
                    if(el_->type_tmp==1){ //function
                        all_terminals = false;
                        break;
                    }
                }

                if(all_terminals){
                    int w = 1;
                    for(auto& el_ : el->ch){
                        w+=el_->w;
                    }
                    el->w = w;

                    el->type_tmp = 0;
                }
            }

        }

        if(!any_funcs_left)
            break;
    }
    
}

void Ttree::init_step1_i(){
    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);

    int max_depth=0;
    for(auto& el : els){
        if(el==nullptr)
            continue;

        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));

        int depth = el->depth;
        if(depth>max_depth) max_depth = depth;
        if(depth==0) el->i0=0;
    }


    for(int depth=0;depth<=max_depth;++depth){
        for(auto& el : els){
            if(el->depth != depth) continue;

            int i0_par = el->i0;
            int w_right_sib_sum = 0;
            int w_right_ch_sum = 0;  // for i
            int n_ch = el->ch.size();
            
            int iii_ = 0;
            for(auto& el_ : el->ch){
                int i_left;
                if(iii_++<n_ch/2){
                    i_left=0;
                    w_right_ch_sum+=el_->w;
                }
                else{
                    i_left=1;
                }
                el_->i0 = i0_par+w_right_sib_sum + i_left;
                el_->i_left = i_left;
                w_right_sib_sum += el_->w;
            }
            el->i = el->i0 + w_right_ch_sum;
        }
    }
}

void Ttree::init_step2_cells(const int cell_width){

    std::string fmt_s0 = "%-"+std::to_string(cell_width)+"s";
    std::string s_tmp(cell_width,' ');


    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);

    int max_depth=0;
    for(auto& el : els){
        if(el==nullptr)
            continue;

        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));

        int depth = el->depth;
        if(depth>max_depth) max_depth = depth;
    }

    /*
       
 columns              
(d means depth):       d=0     d=1    d=2   
rows:
    i=0              |*-----|x0    |      |
    i=1              |x1    |      |      |
    i=2              |!     | -----|x2    |
    i=3              |*-----|x3    |      |
    i=4              |      | -----|x4    |

    */

    // each el owns a row of cells = std::vector<std::string>
    // vector.size=el->depth+1, i.e. the el's value is always written to the last vectors's slot

    for(auto& el : els){
        int depth = el->depth;
        el->cells.resize(depth+1);
        el->cells[depth].resize(cell_width);

        std::string value_decor = el->get_value_decor();

        int n = std::min(size_t(cell_width),value_decor.size());
        std::snprintf(&(el->cells[depth].front()),n+1,fmt_s0.c_str(),value_decor.c_str());
        
        for(int d=0;d<depth;++d){
            el->cells[d] = s_tmp;
        }

        // by now, cells[row=i,col=0..depth] have width=cell_width

        // see fig. above
        // if it's e.g. x4, then now we need to draw connections to it's parent, x1.
        // symbols: "*" in the fig. is u2514 and u250c, and "-" in the fig. is u2500
        if(depth>0){
            std::string s0;
            if(el->i_left)
                s0 = "  \u2514";
            else
                s0 = "  \u250c";
            for(int i=0;i<cell_width-3;++i){
                s0 += "\u2500";
            }
            el->cells[depth-1] = s0;
        }
    }
    
    // now we need to draw vertical parent connections 
    // for childs sitting further than 1 row away from the parent
    // i.e. symbol "!" in the fig. = u2502 below

    const std::string vcon = "  \u2502";
    const std::string vcon_ = "  \u251c";

    for(int depth=1;depth<=max_depth;++depth){
        for(auto& el : els){
            if(el->depth != depth) continue;

            int par_i = el->parent.lock()->i;
            int i = el->i;

            if(std::abs(par_i-i)>1){
                for(int i_=std::min(i,par_i)+1;i_<std::max(i,par_i);++i_){
                    for(auto& el_ : els){
                        if(el_->i != i_) continue;

                        if(trim(el_->cells[depth-1]).size()<vcon.size()){
                            el_->cells[depth-1] = vcon;
                            el_->cells[depth-1].append(cell_width-3,' ');
                            // vcon.size=5, but only 3 printed chars
                        }
                        else{
                            el_->cells[depth-1].replace(0,vcon_.size(),vcon_,0,vcon_.size());
                        }
                    }
                }
            }
        }
    }


}


void Ttree::viz_init(){
    init_step0_w();
    init_step1_i();
    init_step2_cells(8);
}

void Ttree::viz(int show_ind){
    viz_init();

    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);

    int max_depth=0;
    int max_i=0;
    for(auto& el : els){
        if(el==nullptr)
            continue;

        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));

        int depth = el->depth;
        if(depth>max_depth) max_depth = depth;
        int i = el->i;
        if(i>max_i) max_i = i;
    }

    for(int i=0;i<=max_i;++i){
        for(int depth=0;depth<=max_depth;++depth){
            int i000=0;
            for(auto& el: els){
                if(el->i==i && el->depth==depth){
                    for(auto& s: el->cells) std::cout << s;
                    if(show_ind==1) std::cout << "  ["<< el->i <<"]" ;
                    if(show_ind==2) std::cout << "  ["<< i000 <<"]" ;
                    std::cout << std::endl;
                }
                i000++;
            }
        }
    }
}

bool Ttree::delete_me_from_ch(const std::shared_ptr<Tel> el, bool nullify){

    bool ret = false;

    std::weak_ptr<Tel> par = el->parent;
    if(!par.expired()){
        std::list<std::shared_ptr<Tel>> &ch_ = par.lock()->ch;
        auto pos = std::find(ch_.begin(),ch_.end(),el);
        if(pos==ch_.end()){
            throw xX_0("mf-error");
        }

        if(!nullify)
            ch_.erase(pos);
        else
            *pos = nullptr;
        ret = true;
    }

    return ret;
}

void Ttree::cut(const int cpoint, const bool by_els_ind){

    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);
    for(auto& el : els){
        if(el==nullptr)
            continue;
        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
    }


    for(auto el=els.begin();el!=els.end();++el){

        if(
                (by_els_ind && std::distance(els.begin(),el)==cpoint)
                ||
                (!by_els_ind && (*el)->i==cpoint)
                )
        {
            delete_me_from_ch(*el);
        }
    }
}


int Ttree::get_max_depth_cur() const {
    int max_depth=0;

    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);

    for(auto& el : els){
        if(el==nullptr)
            continue;

        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));

        int depth = el->depth;
        if(depth>max_depth) max_depth = depth;
    }
    return max_depth; 
};
int Ttree::get_depth_by_cp(const int cp) const {

    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);

    int i00=-1;
    int depth = -1;
    for(auto& el : els){
        i00++;
        if(el==nullptr)
            continue;
        if(i00==cp){
            depth = el->depth;
            break;
        }
        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
    }
    return depth;
}

std::unordered_set<int> Ttree::get_cps_by_depth(const int depth_min,const int depth_max) const {

    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);

    std::unordered_set<int> cps;
    int i00=-1;
    for(auto& el : els){
        i00++;
        if(el==nullptr)
            continue;

        int depth0 = el->depth;
        if(
                (depth_max>=0 && depth0>=depth_min && depth0<=depth_max)
                ||
                ( depth_max<0 && depth0>=depth_min)
          )
        {
            cps.insert(i00);
        }
        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
    }
    return cps;
}

void Ttree::tst(){

    std::list<std::shared_ptr<Tel>> els;
    els.push_back(root);
    for(auto& el : els){
        if(el==nullptr)
            continue;
        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
    }


    int cpoint = 3;
    for(auto el=els.begin();el!=els.end();++el){

        if(std::distance(els.begin(),el)==cpoint)     {
            std::weak_ptr<Tel> par = (*el)->parent;
            if(!par.expired()){
                std::list<std::shared_ptr<Tel>> &ch_ = par.lock()->ch;
                auto pos = std::find(ch_.begin(),ch_.end(),*el);
                if(pos==ch_.end()){
                    throw xX_0("mf-error");
                }
                //ch_.erase(pos);
                *pos = nullptr;
            }
        }
    }


    {
        std::list<std::shared_ptr<Tel>> els;
        els.push_back(root);
        for(auto& el : els){
            pr("is nullptr?",el==nullptr);
            if(el != nullptr)
            std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
        }
    }

}

std::vector<double> Ttree::predict_proba(const tvvd &X){
    std::vector<double> y_pred = predict(X);
    for(auto &y : y_pred){
        y = 1/(1+std::exp(-y));
    }
    return y_pred;
}

std::vector<double> Ttree::predict(const tvvd &X){

    std::vector<double> y_pred(X.size());
    for(size_t i=0;i<X.size();++i){
        double y0 = eval(X.at(i));
        y_pred.at(i)=y0;
    }

    return y_pred;
}

void Ttree::eval_v(const txy &xy){

    std::vector<double> y_pred = predict(xy.first);
    err = J_(xy.second,y_pred);
}


tvd Ttree::predict_bin(const tvvd &X,const double thr){

    tvd y_pred = predict_proba(X);

    for(auto &y: y_pred){
        if(y<thr)
            y=-1;
        else
            y=1;
    }

    return y_pred;
}




std::string Ttree::to_string(){


    std::list<std::shared_ptr<Tel>> els;
    if(root!=nullptr)
        els.push_back(root);
    for(auto& el : els){
        if(el==nullptr)
            continue;

        std::copy(el->ch.begin(),el->ch.end(),std::back_inserter(els));
        el->type_tmp = el->type;
        if(el->type==0) //terminal
            el->expr = el->get_value_decor();
    }


    int iii = 0;
    while(true){
        
        if(iii++>1000){
            pr("iii=",iii);
            throw xX_0("mf-error");
        }


        bool any_funcs_left = false;
        for( auto& el : els){
            if(el==nullptr)
                continue;

            if(el->type_tmp==1){
                any_funcs_left = true;

                bool all_terminals = true;
                for(auto& el_ : el->ch){
                    if(el_->type_tmp==1){ //function
                        all_terminals = false;
                        break;
                    }
                }

                if(all_terminals){

                    std::list<std::shared_ptr<Tel>> &chl = el->ch;
                    if(!chl.size())
                        throw xX_0("to_string: no ch for "+el->get_value_decor());
                    std::string sep;

                    el->expr = "";
                    if(el->func->isextra_br()) el->expr += "(";

                    if(el->func->isprefix()){
                        el->expr += el->get_value_decor();
                        sep = ",";
                    }
                    else{
                        sep = el->get_value_decor();
                    }


                    auto last = std::prev(chl.end());
                    el->expr += "(";

                    for(auto it=chl.begin();it!=chl.end();++it){
                        el->expr += (*it)->expr;
                        if(it!=last)
                            el->expr+=sep;
                    }

                    el->expr += ")";
                    if(el->func->isextra_br()) el->expr += ")";


                    el->type_tmp = 0;
                }
            }

        }

        if(!any_funcs_left)
            break;
    }

    std::string ret;
    if(root!=nullptr)
        ret = reduce_br(root->expr);

    return ret;
}





std::map<size_t,size_t> Ttree::match_br(const std::string &s){

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

std::string Ttree::reduce_br(const std::string &s){

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


std::map<size_t,size_t> Ttree::split_m(const std::string &s,const std::string &sep){

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


std::map<size_t,size_t> Ttree::top_level_br(const std::map<size_t,size_t> &br){

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





