#ifndef TFUNC_H_
#define TFUNC_H_


#include <vector>
#include <string>
#include <cmath>
#include <functional>


class Tfunc_base{
    protected:
        Tfunc_base(const int num_of_args, const bool isprefix_,
                const std::string &decor_value,
                const std::function<double(const std::vector<double>&)> &l_f,
                const bool isextra_br_=false
                ): 
            num_of_args(num_of_args),isprefix_(isprefix_),decor_value(decor_value),l_f(l_f),
            isextra_br_(isextra_br_)
        {
            //pr("protected constr called",this);
        }

    public:
        Tfunc_base(){
            //pr("base constr called",this);
        }
        virtual ~Tfunc_base() noexcept {};
        double f(const std::vector<double> &x) const {
            return l_f(x);
        }

        std::string get_value_decor() const {
            return decor_value;
        };

        int get_num_args() const {
            return num_of_args;
        }
        bool isprefix()const {
            return isprefix_;
        }
        bool isextra_br()const {
            return isextra_br_;
        }


        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }

        static constexpr const char* class_name = "Tfunc_base";


    private:
        int num_of_args;
        bool isprefix_;
        std::string decor_value;
        std::function<double(const std::vector<double>&)> l_f;
        bool isextra_br_;
};

class Tfunc_plus : public Tfunc_base {
    public:
        Tfunc_plus():
            Tfunc_base(2,false,"+",
                  [](const std::vector<double> &x) -> double { return x[0]+x[1]; }
                    ) 
            {   /*pr("Tfunc_plus constr called",this);*/  }

        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_minus : public Tfunc_base {
    public:
        Tfunc_minus():
            Tfunc_base(2,false,"-",
                  [](const std::vector<double> &x) -> double { return x[0]-x[1]; }
                    ) 
            { /* pr("Tfunc_minus constr called",this); */  }
        
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_pref_minus : public Tfunc_base {
    public:
        Tfunc_pref_minus():
            Tfunc_base(1,true,"-",
                  [](const std::vector<double> &x) -> double { return -x[0]; },
                  true
                    ) 
            {}
        
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_mult : public Tfunc_base {
    public:
        Tfunc_mult():
            Tfunc_base(2,false,"*",
                  [](const std::vector<double> &x) -> double { return x[0]*x[1]; }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_div : public Tfunc_base {
    public:
        Tfunc_div():
            Tfunc_base(2,false,"/",
                  [](const std::vector<double> &x) -> double { 
                  //return std::abs(x[1])<1e-7 ? 1.0 :  x[0]/x[1]; }
                  return std::abs(x[1])<1e-4 ? 1e4 :  x[0]/x[1]; }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_log : public Tfunc_base { // ln
    public:
        Tfunc_log():
            Tfunc_base(1,true,"log",
                  [](const std::vector<double> &x) -> double { 
                  return x[0]<1e-2 ? -4.6 :  std::log(x[0]); }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_sin : public Tfunc_base {
    public:
        Tfunc_sin():
            Tfunc_base(1,true,"sin",
                  [](const std::vector<double> &x) -> double { return std::sin(x[0]); }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_cos : public Tfunc_base {
    public:
        Tfunc_cos():
            Tfunc_base(1,true,"cos",
                  [](const std::vector<double> &x) -> double { return std::cos(x[0]); }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};

class Tfunc_tanh : public Tfunc_base {
    public:
        Tfunc_tanh():
            Tfunc_base(1,true,"tanh",
                  [](const std::vector<double> &x) -> double { return std::tanh(x[0]); }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};

class Tfunc_pow2 : public Tfunc_base {
    public:
        Tfunc_pow2():
            Tfunc_base(1,true,"pow2",
                  [](const std::vector<double> &x) -> double { return x[0]*x[0]; }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_pow3 : public Tfunc_base {
    public:
        Tfunc_pow3():
            Tfunc_base(1,true,"pow3",
                  [](const std::vector<double> &x) -> double { return x[0]*x[0]*x[0]; }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_min : public Tfunc_base {
    public:
        Tfunc_min():
            Tfunc_base(2,true,"min",
                  [](const std::vector<double> &x) -> double { return std::min(x[0],x[1]); }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_max : public Tfunc_base {
    public:
        Tfunc_max():
            Tfunc_base(2,true,"max",
                  [](const std::vector<double> &x) -> double { return std::max(x[0],x[1]); }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};
class Tfunc_list2 : public Tfunc_base {
    public:
        Tfunc_list2():
            Tfunc_base(2,false,",",
                  [](const std::vector<double> &x) -> double { return 0; }
                    ) {}
        virtual auto new_this() -> decltype(this) {
            return new std::remove_reference<decltype(*this)>::type;
        }
};

#endif /* TFUNC_H_ */

