#include "ttrees.h"

#include <iomanip>
#include <set>
#include <algorithm>
#include <chrono>
using namespace std::chrono;



std::default_random_engine dre;



Ttrees::Ttrees(const Ttrees_parameters p,const LOG logging,int seed):
    p(p),logging(logging),num_of_co(0),min_err(1e10),urd(0,1),xy(xy_default_initializer)
{
    if(seed>=0)
        dre.seed(seed);

    for(int iii=0;iii<p.ntrees;++iii){
        std::unique_ptr<Ttree> tree0(new Ttree(p.p));
        trees.push_back(std::move(tree0));
    }
};



void Ttrees::fit(const txy &xy, Treport_detail *rd){

    auto tp0 = system_clock::now();

    this->xy = xy;

    min_err = 1e10;
    num_of_co=0;
    errs.reserve(p.ngen_max);


    if(logging>LOG::SILENT) pr("");
    for(int gen=0;gen<p.ngen_max;++gen){

        /* perform cross-over operation with 95% chance */
        if(urd(dre) < 0.95){
            int nc0;
            co(nc0);
            num_of_co += nc0;
            errs.push_back(min_err);
        }
        else{
            /* and mutation with only 5% of chance */
            mut();
        }
        if(min_err<p.err_thr)
            break;


        if(logging>LOG::SILENT && p.ngen_max>1){
            /* decoration of the progress report */
            double pf_ = 100./(p.ngen_max-1);
            int perc = gen*pf_;
            std::cout<<"\r";
            std::cout<<"fit progress: "<<std::setw(3)<<perc<<"%";
            std::cout.flush();
        }
    }
    if(logging>LOG::SILENT) pr("");

    auto dur = duration_cast<milliseconds>(system_clock::now()-tp0);
    double dur_s = dur.count()/1000.0;

    if(rd!=nullptr){
        /* saving some run stats for later reports */
        rd->dur = dur_s;
        rd->num_of_co=num_of_co;
        rd->min_err = min_err;
    }
}

Ttree Ttrees::best_tree(){
    
    double err_ = 1e10;
    int i_=-1;
    for(size_t i=0;i<trees.size();++i){
        if(trees[i]->err<0)
            continue;

        if(trees[i]->err<err_){
            err_=trees[i]->err;
            i_=i;
        }
    }
    if(i_<0)
        throw xX_0("best_tree: empty list");

    return *trees[i_];
}

void Ttrees::print_depth_cnt(const std::string &msg){
    std::map<int,int> depth_cnt;
    for(const auto& t : trees){
        int d = t->get_max_depth_cur();
        depth_cnt[d]++;
    }
    pr(msg);
    for(const auto u : depth_cnt){
        pr(u.first,u.second);
    }
    pr("");
}
void Ttrees::print_use_cnt(const std::string &msg){
    std::map<int,int> use_cnt;
    for(const auto& t : trees){
        use_cnt[t->use_count]++;
    }
    pr(msg);
    for(const auto u : use_cnt){
        pr(u.first,u.second);
    }
    pr("");
}
void Ttrees::print_ind_cnt(const std::string &msg){
    //pr("ind cnt:");
    //for(const auto u : ind_cnt){
    //    pr(u.first,u.second);
    //}
    //pr("");
}



std::vector<size_t> Ttrees::get_sorted_subset(const std::unordered_set<size_t> &excl){

    int ntrees = trees.size();
    if(!ntrees)
        throw xX_0("get_sorted_subset");

    std::uniform_int_distribution<size_t> uid(0,ntrees-1);

    std::unordered_set<size_t> i_ss0;

    int N = std::max(1,int(p.ratio_ss*ntrees));

    for(int i=0;i<N;++i){
        size_t rv = uid(dre);
        if(excl.find(rv)==excl.end()){
            i_ss0.insert(rv);
            Ttree &t = *trees[rv];
            if(t.err<0)
                t.eval_v(xy);
        }
    }
    std::vector<size_t> ii(i_ss0.begin(),i_ss0.end());

    auto less = [&](size_t i0, size_t i1){return trees[i0]->err < trees[i1]->err;};
    std::sort(ii.begin(),ii.end(),less);

    return ii;
}
Ttree Ttrees::co_(const Ttree &tree0, const Ttree &tree1){

    int cp1 = get_discrete_s(1,tree1.count()-1,p.lam,dre);
    __debug__cnt_cp1[cp1]++;

    int depth1 = tree1.get_depth_by_cp(cp1);
    int depth_left = p.max_depth_ex-depth1;
    depth_left = std::max(depth_left,1);

    int max_depth_cur0 = tree0.get_max_depth_cur();
    int d0 = max_depth_cur0-depth_left;
    std::unordered_set<int> cps0 = tree0.get_cps_by_depth(d0+1);

    int cp0_ind = get_discrete_s(0,cps0.size()-1,p.lam,dre);
    auto it = cps0.begin();
    std::advance(it,cp0_ind);
    int cp0 = *it;


    Ttree tree0_(tree0,cp0);
    Ttree tree1_(tree1);
    tree1_.replace(cp1,tree0_);

    //pr("cp1=",cp1,"depth1=",depth1,"depth_left=",depth_left,"d0=",d0,"cp0=",cp0);

    //std::cout << cp0 <<","<< cp1 << "   ";
    //std::cout << tree1_.get_max_depth_cur() << "   ";
    //pr("new tree count=",cp0+(tree1count-n1),tree1_);

    return tree1_;
}
int Ttrees::co_and_add_best(const size_t i0, const size_t i1, const bool erase_worst){

    Ttree &t0 = *trees[i0];
    Ttree &t1 = *trees[i1];

    int i_worst = -1;
    if(t1.count()<2)
        return i_worst;

    ind_cnt[i0]++;
    ind_cnt[i1]++;

    std::unique_ptr<Ttree> tree2(new Ttree(co_(t0,t1)));


    t0.use_count++;
    t1.use_count++;

    Ttree &t2 = *tree2;
    t2.eval_v(xy);
    t2.use_count++;

    if(t0.err>=t1.err){
        if(t0.err>=t2.err){//t0 highest (worst)
            i_worst = i0;
            trees.push_back(std::move(tree2));
        }
        else{ /*t2 highest*/ }
    }
    else{ //t1 > t0
        if(t1.err>=t2.err){ //t1 highest
            i_worst = i1;
            trees.push_back(std::move(tree2));
        }
        else{ /*t2 highest */ }
    }
    double min_err0 = std::min(std::min(t0.err,t1.err),t2.err);
    if(min_err0<min_err)
        min_err = min_err0;

    if(erase_worst)
        if(i_worst>=0)
            trees.erase(std::next(trees.begin(),i_worst));

    return i_worst;
}
/* cross-over operation  */ 
void Ttrees::co(int &num_of_co){

    std::unordered_set<size_t> excl;
    std::vector<size_t> ii = get_sorted_subset(excl);
    if(ii.size()<1)
        throw xX_0("co: get_sorted_subset resulted in empty list, try increasing \"ratio_ss\"");

    size_t i0 = ii[0];

    //excl.insert(i0);
    ii = get_sorted_subset(excl);
    if(ii.size()<1)
        throw xX_0("co: get_sorted_subset resulted in empty list, try increasing \"ratio_ss\"");

    size_t i1 = ii[0];

    int i_worst = co_and_add_best(i0,i1,true);
    if(i_worst>=0)
        num_of_co = 1;
    else
        num_of_co = 0;
}


Ttree Ttrees::mut_(const Ttree &tree0){

    int cp0 = get_discrete_s(0,tree0.count()-1,p.lam,dre);
    int depth0 = tree0.get_depth_by_cp(cp0);

    int depth_left = p.max_depth_ex-depth0;
    depth_left = std::max(depth_left,0);

    Ttree_parameters p_(p.p);
    p_.max_depth = depth_left;
    Ttree tree1_(p_); // create new random tree

    Ttree tree0_(tree0); // copy tree0
    tree0_.replace(cp0,tree1_);

    return tree0_;
}
void Ttrees::mut(){

    size_t ntrees = trees.size();
    if(!ntrees)
        throw xX_0("mut");

    std::uniform_int_distribution<size_t> uid(0,ntrees-1);
    size_t i0 = uid(dre);
    Ttree &t0 = *trees[i0];
    if(t0.err<0)
        t0.eval_v(xy);

    std::unique_ptr<Ttree> tree1(new Ttree(mut_(t0)));
    Ttree &t1 = *tree1;
    t1.eval_v(xy);

    if(t0.err>t1.err){
        trees.erase(std::next(trees.begin(),i0));
        trees.push_back(std::move(tree1));

        if(t1.err<min_err)
            min_err = t1.err;
    }
}






