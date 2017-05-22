#ifndef TPD_H_
#define TPD_H_

#include <string>
#include <vector>

/*
 * Simiplistic class for reading csv files
 *
 * */

class Tpd{
    public:
        std::vector<std::vector<double>> data;
        std::vector<std::string> columns;

        void read_csv(const std::string &filename,
                const bool has_header);

    private:
        inline static void skip_ws(std::ifstream &file);
        inline static std::vector<std::string> fill_(std::ifstream &file,int &nc);

};


#endif /* TPD_H_ */
