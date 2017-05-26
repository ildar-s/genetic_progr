#ifndef METRICS_H_
#define METRICS_H_

#include "common_0.h"

namespace metrics{

    inline double accuracy(const tvd &y_true, const tvd &y_pred){
        int s = 0;
        size_t n = y_true.size();
        for(size_t i=0;i<n;++i){
            s += static_cast<int>(y_true[i]==y_pred[i]);
        }

        return s*1.0/n;
    }

    inline double J_squared(const tvd &y0,const tvd &y1){

        double s=0;
        for(size_t i=0;i<y0.size();++i){
            double z_ = y0.at(i)-y1.at(i);
            s += z_*z_;
        }
        return s;
    }

    // y_true must be {-1,1}
    inline double J_logit(const tvd &y_true, const tvd &y_pred){

        double s=0;
        for(size_t i=0;i<y_true.size();++i){

            double z_ = y_true.at(i) * y_pred.at(i);

            s += std::log(1+std::exp(-z_));
        }
        return s;
    }

    inline double J(const int loss_type, const tvd &y_true, const tvd &y_pred){
        if(loss_type==0){
            return J_squared(y_true,y_pred);
        }
        else{
            return J_logit(y_true,y_pred);
        }
    }


}
#endif /* METRICS_H_ */

