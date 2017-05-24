#ifndef TTREE_H_
#define TTREE_H_

#include "tel.h"
#include "tfunc.h"
#include "common_0.h"

#include <unordered_set>
#include <map>
#include <random>
#include <memory>



struct Ttree_parameters{
    Ttree_parameters():
        ndim(0),
        max_depth(2),
        element_function_probability(0.85),
        terminal_const_probability(0.2),
        logit(false)
    {}
    Ttree_parameters(const Ttree_parameters &src):
        ndim(src.ndim),
        max_depth(src.max_depth),
        element_function_probability(src.element_function_probability),
        terminal_const_probability(src.terminal_const_probability),
        logit(src.logit)
    {}
    void print_all(){
        pr("ndim =",ndim);
        pr("max_depth =",max_depth);
        pr("element_function_probability =",element_function_probability);
        pr("terminal_const_probability =",terminal_const_probability);
        pr("logit =",logit);
    }

    int ndim;
    int max_depth;
    double element_function_probability;
    double terminal_const_probability;
    bool logit;
};

class Ttree{
    public:
        friend class Tel;
        Ttree(const int ndim);
        Ttree(const Ttree_parameters &p,const std::function<double(const std::vector<double>&, 
                const std::vector<double>&)> &J=nullptr);
        Ttree(const Ttree &src, const int cpoint=-1);
        Ttree(const std::string &s, const Ttree_parameters * const p_=nullptr);

        void update_ndim(const int ndim);

        void assign_els(const std::shared_ptr<Tel> &root_src, std::shared_ptr<Tel> &root_dest);
        void replace(const int cpoint, const Ttree &src);
        bool delete_me_from_ch(const std::shared_ptr<Tel> el, bool nullify=false);

        double eval(const std::vector<double> &x); // x.size() = num_of_d
        void viz(int show_ind=0);
        void cut(const int cpoint, const bool by_els_ind=false);
        std::shared_ptr<Tel> copy(const int cpoint);
        void tst();
        int count() const;
        void list_addr() const;
        int get_max_depth() const { return p.max_depth; };
        int get_max_depth_cur() const;
        int get_depth_by_cp(const int cp) const;
        std::unordered_set<int> get_cps_by_depth(const int depth_min,const int depth_max=-1) const;

        static void init_static();

        void eval_v(const txy &xy);
        tvd predict(const tvvd &X);
        tvd predict_proba(const tvvd &X);
        tvd predict_bin(const tvvd &X,const double thr=0.5);

        std::string to_string();

        double err;
        int use_count;

        static std::vector<double> consts;

        std::uniform_real_distribution<double> urd;
        std::unique_ptr<std::uniform_int_distribution<int>> p_uid_term;
        std::unique_ptr<std::uniform_int_distribution<int>> p_uid_term_c;
        std::unique_ptr<std::uniform_int_distribution<int>> p_uid_term_v;
        std::unique_ptr<std::uniform_int_distribution<int>> p_uid_func;


        static double J_default(const std::vector<double> &y0,const std::vector<double> &y1);
        static double J_logit(const std::vector<double> &y_true,const std::vector<double> &y_pred);

        static std::map<size_t,size_t> match_br(const std::string &s);
        static std::string reduce_br(const std::string &s);
        static std::map<size_t,size_t> split_m(const std::string &s,const std::string &sep);
        static std::map<size_t,size_t> top_level_br(const std::map<size_t,size_t> &br);

        static std::vector<std::shared_ptr<Tfunc_base>> funcs;
        static std::vector<std::shared_ptr<Tfunc_base>> funcs_l;

        static bool static_up_to_date;

    private:
        void init_step0_w();
        void init_step1_i();
        void init_step2_cells(const int cell_width=8);
        void viz_init();

        std::shared_ptr<Tel> root;

        Ttree_parameters p;

        std::vector<double> x;
        std::vector<std::pair<double*,std::string>> terms;

        //static int ndim_st;
        std::function<double(const std::vector<double>&, const std::vector<double>&)> J_;
};

#endif /* TTREE_H_ */
