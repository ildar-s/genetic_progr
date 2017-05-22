#ifndef TTREES_H_
#define TTREES_H_

#include "ttree.h"
#include "common_ttrees.h"
#include "common_0.h"

#include <map>

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
};

struct Treport_detail{
    double dur;
    double min_err;
    int num_of_co;
};

class Ttrees{
    public:
        Ttrees(const Ttrees_parameters p,int seed=-1);

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

        std::vector<std::unique_ptr<Ttree>> trees;
        std::map<size_t,int> ind_cnt;
        int num_of_co;
        double min_err;
        std::uniform_real_distribution<double> urd;

        std::reference_wrapper<const txy> xy;
        txy xy_default_initializer;
};


#endif /* TTREES_H_ */

