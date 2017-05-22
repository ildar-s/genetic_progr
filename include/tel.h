#ifndef TEL_H_
#define TEL_H_

#include "tfunc.h"
#include "ttk.h"

#include <list>
#include <memory>


class Ttree;

class Tel{
    public:
        friend class Ttree;
        Tel(const int depth,const std::weak_ptr<Tel>& parent,
                Ttree *ptree);
        Tel(const int depth,const std::weak_ptr<Tel>& parent, 
                const Tel& src, const Ttree *ptree);
        Tel(const int depth,const std::weak_ptr<Tel> &parent,
                std::unique_ptr<Ttk> tk_);

        std::string get_value_decor() const;

        ~Tel(){}



    private:
        int depth;
        int type;
        int value_ind;

        Tfunc_base* func;
        std::pair<double*,std::string> term;

        int type_tmp;
        double value_tmp;
        int w,i_left,i0,i;
        std::vector<std::string> cells;
        std::string expr;

        std::weak_ptr<Tel> parent;
        std::list<std::shared_ptr<Tel>> ch;
        std::unique_ptr<Ttk> tk;
        double cvalue; // TMP!!!
};



#endif /* TEL_H_ */
