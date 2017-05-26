#ifndef TTREES_H_
#define TTREES_H_

#include "ttree.h"
#include "common_ttrees.h"
#include "common_0.h"

#include <map>
#include <fstream>
#include <iomanip>

struct Ttrees_parameters{
    Ttrees_parameters():
        ntrees(1000),
        max_depth_ex(5),
        ngen_max(30000),
        err_thr(150),
        lam(0.1),
        ratio_ss(0.01)
    {}

    int ntrees;
    int max_depth_ex;
    int ngen_max;
    double err_thr;
    double lam;
    double ratio_ss;
    Ttree_parameters p;

    void print_all(){
        pr("ntrees =",ntrees);
        pr("max_depth_ex =",max_depth_ex);
        pr("ngen_max =",ngen_max);
        pr("err_thr =",err_thr);
        pr("lam =",lam);
        pr("ratio_ss =",ratio_ss);
        p.print_all();
    }
    void save(const bool is_bin,const std::string &filename){
        std::ofstream(filename,std::ios::out|std::ios::binary)<< 
            is_bin << std::endl <<
            ntrees << std::endl <<
            max_depth_ex << std::endl <<
            ngen_max << std::endl <<
            err_thr << std::endl <<
            lam << std::endl <<
            ratio_ss << std::endl <<
            p.ndim << std::endl <<
            p.max_depth << std::endl <<
            p.element_function_probability << std::endl <<
            p.terminal_const_probability << std::endl <<
            p.logit << std::endl <<
            p.consts_min << std::endl <<
            p.consts_max << std::endl <<
            p.consts_n << std::endl <<
            p.loss_type;
    }
    void load(const std::string &filename, bool &is_bin){
        std::ifstream(filename,std::ios::in)>> 
            is_bin >>
            ntrees >> 
            max_depth_ex >> 
            ngen_max >> 
            err_thr >> 
            lam >> 
            ratio_ss >> 
            p.ndim >> 
            p.max_depth >> 
            p.element_function_probability >> 
            p.terminal_const_probability >> 
            p.logit >>
            p.consts_min >>
            p.consts_max >>
            p.consts_n >>
            p.loss_type;
    }

};

struct Treport_detail{
    double dur;
    double min_err;
    int num_of_co;
};

class Ttrees{
    public:
        Ttrees(const Ttrees_parameters p,const LOG logging=LOG::INFO,int seed=-1);

        void fit(const txy &xy,Treport_detail *rd=nullptr);

        Ttree best_tree();

        void print_depth_cnt(const std::string &msg);
        void print_use_cnt(const std::string &msg);
        void print_ind_cnt(const std::string &msg);

        std::vector<double> errs;
        std::map<int,int> __debug__cnt_cp1;

    private:
        std::vector<size_t> get_sorted_subset(const std::unordered_set<size_t> &excl);
        Ttree co_(const Ttree &tree0, const Ttree &tree1);
        int co_and_add_best(const size_t i0, const size_t i1, const bool erase_worst);
        void co(int &num_of_co);
        void co_full_elite_sel();
        void co_full_rand_sel();
        void co_rand_sel();
        Ttree mut_(const Ttree &tree0);
        void mut();


        Ttrees_parameters p;
        LOG logging;

        std::vector<std::unique_ptr<Ttree>> trees;
        std::map<size_t,int> ind_cnt;
        int num_of_co;
        double min_err;
        std::uniform_real_distribution<double> urd;

        std::reference_wrapper<const txy> xy;
        txy xy_default_initializer;
};


#endif /* TTREES_H_ */

