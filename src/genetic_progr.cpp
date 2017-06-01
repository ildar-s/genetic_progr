#include "tfunc.h"
#include "metrics.h"
#include "ttree.h"
#include "ttrees.h"
#include "common_0.h"
#include "tpd.h"



#include <set>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stack>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std::chrono;




system_clock::time_point make_time_point(int day,int mon,
        int year,int hour,int min,int sec=0){

    struct std::tm t;
    t.tm_mday = day;
    t.tm_mon = mon-1;
    t.tm_year = year-1900;
    t.tm_hour = hour;
    t.tm_min = min;
    t.tm_sec = sec;
    t.tm_isdst = -1;
    std::time_t tt = std::mktime(&t);
    if(tt==-1)
        throw xX_0("make_time_point: make_time failed");
    return system_clock::from_time_t(tt);
}




/* must be csv file, sep=',', first column always ignored  */
txy read_xy(bool &is_bin, const std::string &filename, const bool has_y,
        const bool has_header=true, const LOG logging=LOG::INFO,
        std::vector<std::string> *columns=nullptr){

    Tpd pd;
    txy xy;
    try{
        pd.read_csv(filename,has_header);

        if(pd.data.size()){
            if((has_y && pd.data[0].size()<3)||(!has_y && pd.data[0].size()<2))
                throw xX_0("bad size of data array while reading <"+filename+">");
        }
    }
    catch(const std::exception &e){
        throw xX_0("read_csv failed");
    }


    if(columns!=nullptr)
        *columns = pd.columns;


    tvvd &X = xy.first;
    X.resize(pd.data.size());


    if(has_y){
        tvd &y = xy.second;
        y.resize(pd.data.size());

        for(size_t i=0;i<X.size();++i){
            y[i] = pd.data[i].back();
            X[i].assign(pd.data[i].begin()+1,pd.data[i].end()-1);
        }



        std::set<double> y_unique(y.begin(),y.end());
        if(y_unique.size()==2){
            is_bin=true;

            if(logging==LOG::DEBUG){
                std::cout<<"--> labels in input data: ";
                for(const auto &t: y_unique)std::cout<<t<<" ";
            }

            if(y_unique!=std::set<double>{-1,1}){

                double y0 = *y_unique.begin();

                for(auto &t : y){
                    if(t==y0){
                        t=-1;
                    }
                    else{
                        t=1;
                    }
                }
            }

            if(logging==LOG::DEBUG){
                std::cout<<";  finally using labels: ";
                for(const auto &t: std::set<double>(y.begin(),y.end())) std::cout<<t<<" ";
                pr("");
            }

        }
        else{
            is_bin=false;
        }
    }
    else{
        for(size_t i=0;i<X.size();++i){
            X[i].assign(pd.data[i].begin()+1,pd.data[i].end());
        }
    }

    return xy;
}

void to_csv(const txy &xy, const std::string &filename,
        const std::vector<std::string> * const columns=nullptr){

    const tvvd &X = xy.first;
    const tvd &y = xy.second;

    auto it_x = X.begin();
    auto it_y = y.begin();

    std::ofstream file(filename,std::ios::out|std::ios::binary);

    if(columns!=nullptr){
        std::ostream_iterator<std::string> it{file,","};
        std::copy(columns->begin(),columns->end()-1,it);
        file<<columns->back()<<std::endl;
    }

    std::ostream_iterator<double> it{file,","};

    int idx = 0;
    for(;it_x!=X.end();++it_x,++it_y){
        const tvd &x = *it_x;

        *it++ = idx++;
        std::copy(x.begin(),x.end(),it);
        file << *it_y << std::endl;
    }
}

tvsz bootstrap(const size_t size){

    std::uniform_int_distribution<size_t> uid(0,size-1);
    std::vector<size_t> inds;

    for(size_t i=0;i<size;++i){
        size_t ind = uid(dre);
        inds.push_back(ind);
    }

    return inds;
}

typedef std::tuple<txy,std::set<size_t>,std::set<size_t>> txyss;

/* extracting bootstrapped subset of objects */
txyss get_bs_xy(const txy &xy){
    const tvvd &X = xy.first;
    const tvd &y = xy.second;

    if(X.size()!=y.size())
        throw xX_0("X.size != y.size");

    tvsz inds = bootstrap(y.size());

    txyss xyss;
    txy &xy_ = std::get<0>(xyss);

    for(const auto i: inds){
        xy_.first.push_back(X[i]);
        xy_.second.push_back(y[i]);
    }

    std::set<size_t> &ss_yes = std::get<1>(xyss);
    ss_yes.insert(inds.begin(),inds.end());

    tvsz v_all(inds.size());
    std::iota(v_all.begin(),v_all.end(),0);

    std::set<size_t> &ss_no = std::get<2>(xyss);
    std::set_difference(v_all.begin(),v_all.end(),ss_yes.begin(),ss_yes.end(),
            std::inserter(ss_no,ss_no.begin()));


    return xyss;
}


tvd predict_c(std::vector<Ttree> &ts, const tvvd &X, const bool is_bin){

    tvd y_pred(X.size());

    for(auto &t: ts){
        tvd y_pred0;
        if(is_bin)
            y_pred0 = t.predict_bin(X);
        else
            y_pred0 = t.predict(X);

        auto it_y0 = y_pred0.begin();
        auto it_y = y_pred.begin();
        for(;it_y!=y_pred.end();++it_y0,++it_y){
            *it_y += *it_y0;
        }
    }

    if(is_bin){
        for(auto &y:y_pred){
            if(y<0)
                y=-1;
            else
                y=1;
        }
    }
    else{
        const size_t n = ts.size();
        for(auto &y:y_pred) y /= n;
    }

    return y_pred;
}


void fit(const txy &xy_train, Ttrees_parameters &p, const bool is_bin,
        const bool use_bs,
        const int n_iter, const std::string f_file_name, 
        const LOG logging=LOG::INFO){

    Treport_detail rd;

    p.p.ndim = xy_train.first[0].size();

    std::vector<Ttree> best_trees;

    /* prediction will be accumulated here */
    tvd bx(xy_train.first.size());

    std::ofstream(f_file_name,std::ios::out|std::ios::binary);
    p.save(is_bin,f_file_name+".p");

    for(int iter=0;iter<n_iter;++iter){
        if(logging > LOG::SILENT && n_iter>1){
            pr("\niteration no",iter);
        }



        /* init new forest */
        Ttrees trs(p,logging);


        /* tuple for holding bootstrapped xy data and oob indices */
        txyss xyss;

        /* fitting on either bs subsample or full xy_train */
        if(use_bs){
            xyss = get_bs_xy(xy_train);
            txy &xy_train_bs = std::get<0>(xyss);

            trs.fit(xy_train_bs,&rd);
        }
        else{
            trs.fit(xy_train,&rd);
        }


        if(logging > LOG::SILENT){
            pr("* min err =",rd.min_err);
            pr("fitting time =",rd.dur);
            pr("number of cross-overs =",rd.num_of_co);
        }

        /* collecting best trees */
        Ttree t_best = trs.best_tree();
        best_trees.push_back(t_best);

        std::ofstream(f_file_name,std::ios::app|std::ios::out|std::ios::binary)<<
            t_best.to_string()<<std::endl;


        /*  Steps below are for monitoring of the current quality */


        std::set<size_t> oob_inds;

        if(logging > LOG::SILENT){
            if(use_bs){
                oob_inds = std::get<2>(xyss);
            }
            else{
                for(size_t i=0;i<bx.size();++i) oob_inds.insert(i);
            }
        }



        /* 0) Gathering the oob subsample */
        txy xy_test0;
        if(logging > LOG::SILENT){
            for(const auto i: oob_inds){
                xy_test0.first.push_back(xy_train.first[i]);
                xy_test0.second.push_back(xy_train.second[i]);
            }
        }


        /* 1) Quality of the last best_tree: */
        tvd y_pred;
        if(logging > LOG::SILENT){
            if(is_bin){
                y_pred = t_best.predict_bin(xy_test0.first);
                double acc = metrics::accuracy(xy_test0.second,y_pred);
                pr("acc_oob(last best_tree) =",acc);
            }
            else{
                y_pred = t_best.predict(xy_test0.first);
                double err = metrics::J(p.p.loss_type,xy_test0.second,y_pred);
                pr("err(last best_tree) =",err);
            }
        }


        /* 2) Quality of the ensemble of best_tree-s so far */
        /*    each best_tree only predicts within its oob subsample, but over time, 
         *    the entire sample range [0..bx.size) becomes populated */
        if(logging > LOG::SILENT && n_iter>1){
            /* accumulating prediction: 
             * current best_tree predicts within its oob range (or full range is no bs used) */
            auto it_ind = oob_inds.begin();
            auto it_y = y_pred.begin();
            for(;it_ind!=oob_inds.end();++it_ind,++it_y){
                size_t i = *it_ind;
                bx[i] += *it_y;
            }

            /* quality of the ensemble so far  */
            if(is_bin){
                /* binary view of the cumulative prediction so far */
                tvd bx_b = bx;
                for(auto &t:bx_b) t = t<0 ? -1 : 1;

                double acc = metrics::accuracy(xy_train.second,bx_b);
                pr("acc_oob(ensemble) =",acc);
            }
            else{
                /*  bx  ->  bx average*/
                tvd bx_n = bx;
                for(auto &t : bx_n) t /= (iter+1);

                double err = metrics::J(p.p.loss_type,xy_test0.second,bx_n);
                pr("err(last best_tree) =",err);
            }
        }


    }
}

void regressor(Ttrees_parameters &p, const std::string &d_filename, 
        const bool csv_has_header, const std::string &f_filename,
        const int n_iter, const bool use_bs, const LOG logging=LOG::INFO){
    
    txy xy;
    bool is_bin;
    try{
        xy = read_xy(is_bin,d_filename,true,csv_has_header,logging);
        if(logging == LOG::DEBUG){
            if(is_bin)
                pr("--> binary classification");
            else
                pr("--> regression");
        }
    }
    catch(const std::exception &e){
        return;
    }

    fit(xy, p, is_bin, use_bs, n_iter, f_filename,logging);
}

void predict_test(Ttrees_parameters &p, const std::string &d_filename, 
        const bool csv_has_header, const std::string &i_filename,
        const std::string &p_filename,const LOG logging=LOG::INFO){

    bool is_bin;
    p.load(i_filename+".p",is_bin);
    if(logging == LOG::DEBUG){
        if(is_bin)
            pr("--> binary classification");
        else
            pr("--> regression");
    }

    std::ifstream file(i_filename,std::ios::in);

    std::vector<std::string> t_strings(
            std::istream_iterator<std::string>{file},
             std::istream_iterator<std::string>{});


    std::vector<Ttree> ts;
    for(const auto &s : t_strings){
        Ttree tree_(s,&p.p);
        ts.push_back(tree_);
    }

    bool tmp_;
    txy xy = read_xy(tmp_,d_filename,false,csv_has_header,logging);
    tvvd &X_test = xy.first;

    tvd y_pred = predict_c(ts,X_test,is_bin);

    if(p_filename.size()){
        std::ofstream file(p_filename,std::ios::out|std::ios::binary);
        std::ostream_iterator<double> it{file,"\n"};
        std::copy(y_pred.begin(),y_pred.end(),it);
    }
}


void example_tree_viz(Ttrees_parameters &p){
    Ttree tree(p.p);
    tree.viz();
    std::string s = tree.to_string();
    pr("string:",s);
}


class base_targ_key{
    public:
        virtual void add_options(po::options_description *desc)=0;
        virtual void fill_p(po::variables_map *vm)=0;
};

template<typename T>
class targ_key: public base_targ_key {
    public:
        targ_key(const std::string &arg_name, po::value_semantic *vs,
                const std::string &arg_text,T *ptr):
            arg_name(arg_name),vs(vs),arg_text(arg_text),ptr(ptr)
            {};

        void add_options(po::options_description *desc){
            desc->add_options()
                (arg_name.c_str(),vs,arg_text.c_str())
                ;
        };

        void fill_p(po::variables_map *vm){
            if(vm->count(arg_name)){
                *ptr = (*vm)[arg_name].as<T>();
            }
        }

        std::string arg_name;
        po::value_semantic *vs;
        std::string arg_text;
        T *ptr;
};




int main(int argc, char** argv){
    
    Ttrees_parameters p;

    std::vector<std::unique_ptr<base_targ_key>> arg_keys;


    arg_keys.push_back(
        std::make_unique<targ_key<int>>(
            "ntrees",po::value<int>()->default_value(1000),"number of trees" ,&p.ntrees)
            );
    arg_keys.push_back(
        std::make_unique<targ_key<int>>(
            "max_depth",po::value<int>()->default_value(2),
            "max depth of the tree, in the initial forest",&p.p.max_depth )
            );
    arg_keys.push_back(
        std::make_unique<targ_key<int>>("max_depth_ex",po::value<int>()->default_value(5),
            "max depth of the tree, after transformations",&p.max_depth_ex)
            );
    arg_keys.push_back(
        std::make_unique<targ_key<int>>("ndim",po::value<int>()->default_value(0),
         "dimensionality of the parameter space. Only need to set if no data is provided, "\
         "e.g. for random tree vizualization", &p.p.ndim)
            );
    arg_keys.push_back(
        std::make_unique<targ_key<int>>("ngen_max",po::value<int>()->default_value(30000),
            "max generations",&p.ngen_max)
            );
    arg_keys.push_back(
        std::make_unique<targ_key<double>>("err_thr",po::value<double>()->default_value(150),
            "error threshold used in stop criteria",&p.err_thr)
            );
    arg_keys.push_back(
        std::make_unique<targ_key<double>>("lambda",po::value<double>()->default_value(0.1, "0.1"),
         "is a parameter controlling depth of cross-over points: "\
         "lower lambdas mean higher chance of picking lower-depth "\
         "points thus compensating for the natural abundance of high-depth "\
         "nodes in trees (e.g. leafs)",&p.lam)
            );
    arg_keys.push_back(
        std::make_unique<targ_key<double>>("ratio_ss",po::value<double>()->default_value(0.01, "0.01"),
         "ratio of trees considered for cross-over to the total number of trees: "\
         "the selection among only a subset of the forest gives better chance to "\
         "the not-so-good trees to be selected too.",&p.ratio_ss)
            );
    arg_keys.push_back(
        std::make_unique<targ_key<double>>("func_prob",po::value<double>()->default_value(0.85, "0.85"),
         "probability of creating a function node (as opposed to creating a terminal node) "\
         "in forest initialization",&p.p.element_function_probability)
            );
    arg_keys.push_back(
        std::make_unique<targ_key<double>>("term_const_prob",po::value<double>()->default_value(0.2, "0.2"),
         "probability of creating a constant terminal node (as opposed to creating a variable terminal node) "\
         "in forest initialization",&p.p.terminal_const_probability)
            );

    arg_keys.push_back(
        std::make_unique<targ_key<double>>("consts_min",po::value<double>()->default_value(-10, "-10.0"),
         "minimum value of terminal constants",&p.p.consts_min)
            );
    arg_keys.push_back(
        std::make_unique<targ_key<double>>("consts_max",po::value<double>()->default_value(10, "10.0"),
         "maximum value of terminal constants",&p.p.consts_max)
            );
    arg_keys.push_back(
        std::make_unique<targ_key<int>>("consts_n",po::value<int>()->default_value(20, "20"),
         "number of different terminal constants",&p.p.consts_n)
            );

    std::string arg_loss_function;
    arg_keys.push_back(std::make_unique<targ_key<std::string>>("loss_function",
            po::value<std::string>()->default_value("squared"),
            "loss function type, supported types: squared and logistic",&arg_loss_function)
            );

    po::options_description desc_gp("Parameters of GP algorithm");

    for(auto &k : arg_keys){
        k->add_options(&desc_gp);
    }

    
    po::options_description desc_io("I/O control");

    desc_io.add_options() ("cfg",po::value<std::string>(),
            "GP configuration file, "\
            "if not present, will attempt to read \"default.cfg\", "\
            "if not present either, will use default values");
    desc_io.add_options() ("data,d",po::value<std::string>(),
            "CSV file with data, first column always ignored, "\
            "if it's training data then last column - labels (target)");
    desc_io.add_options() ("final_regressor,f",po::value<std::string>(),
            "final regressor file");
    desc_io.add_options() ("initial_regressor,i",po::value<std::string>(),
            "initial regressor file");
    desc_io.add_options() ("predictions,p",po::value<std::string>(),
            "file to output predictions to");
    desc_io.add_options() ("csv_has_header",po::value<std::string>()->default_value("yes"),
            "yes or no");


    po::options_description desc_wf("Workflow control");

    desc_wf.add_options() ("logging",po::value<std::string>()->default_value("info"),
            "logging level: \"silent\", \"info\", or \"debug\"");
    desc_wf.add_options() ("n_iter",po::value<int>()->default_value(1),
            "number of fitting procedures and resulting decision trees");
    desc_wf.add_options() ("bs", po::bool_switch()->default_value(false), 
            "use bootstrapped subsample for each fitting procedure");
    desc_wf.add_options() ("testonly,t", po::bool_switch()->default_value(false), 
            "ignore label information and just test");
    desc_wf.add_options() ("random_seed",po::value<int>(),"random seed");

    desc_wf.add_options() ("random_tree_demo", po::bool_switch()->default_value(false), 
            "vizualize a random tree and quit");


    po::options_description desc("All options");
    desc.add_options() ("help","help message");

    desc.add(desc_gp).add(desc_io).add(desc_wf);

    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc,argv,desc),vm);
    }
    catch(const std::exception &e){
        pr("couldn't parse command line agruments",e.what());
        return 1;
    }
    po::notify(vm);

    if(vm.count("help")){
        std::cout << desc << "\n";
        return 1;
    }

    std::stack<std::string> cfg_filenames;
    cfg_filenames.push("default.cfg");

    if(vm.count("cfg")){
        cfg_filenames.push(vm["cfg"].as<std::string>());
    }
    try{
        po::store(po::parse_config_file<char>(cfg_filenames.top().c_str(),desc),vm);
    }
    catch(const std::exception &e){
        pr("couldn't read any GP config file, will use defaults");
    }


    for(auto &k : arg_keys){
        k->fill_p(&vm);
    }
    if( (arg_loss_function=="logistic")||(arg_loss_function=="log")||
            (arg_loss_function=="logit")  ){
        p.p.loss_type = 1;
    }
    else{
        p.p.loss_type = 0;
    }


    std::string d_filename,f_filename,i_filename,p_filename;
    if(vm.count("data")){
        d_filename = vm["data"].as<std::string>();
    }
    if(vm.count("final_regressor")){
        f_filename = vm["final_regressor"].as<std::string>();
    }
    if(vm.count("initial_regressor")){
        i_filename = vm["initial_regressor"].as<std::string>();
    }
    if(vm.count("predictions")){
        p_filename = vm["predictions"].as<std::string>();
    }

    bool csv_has_header = vm["csv_has_header"].as<std::string>()=="yes";


    int n_iter = vm["n_iter"].as<int>();
    bool use_bs = vm["bs"].as<bool>();


    std::map<std::string,LOG> lmap{
        {"silent",LOG::SILENT},{"info",LOG::INFO},{"debug",LOG::DEBUG} };

    LOG logging = lmap[vm["logging"].as<std::string>()];


    int i_seed;
    if(vm.count("random_seed")){
        i_seed = vm["random_seed"].as<int>();
    }
    else{
        auto tp_ = make_time_point(01,01,2017,10,0);
        i_seed = duration_cast<seconds>(system_clock::now()-tp_).count();
    }
    dre.seed(i_seed);


    try{
        if(vm["random_tree_demo"].as<bool>()){
            example_tree_viz(p);
            return 0;
        }
        if(vm.count("final_regressor")){
            regressor(p,d_filename,csv_has_header,f_filename,n_iter,use_bs,logging);
        }
        if(vm["testonly"].as<bool>() && (i_filename!="")){
            predict_test(p,d_filename,csv_has_header,i_filename,p_filename,logging);
        }
    }
    catch(...){
        return 1;
    }

    if(logging == LOG::DEBUG){
        pr("\nGP parameters used:");
        p.print_all();
    }


    return 0;
}



