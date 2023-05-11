/*
 * NashGame.cpp
 *
 *  Created on: Nov 28, 2022
 *      Author: Ahmet Tikna
 */

#include "NashGame.h"

NashGame::NashGame() {
        // TODO Auto-generated constructor stub

}

NashGame::NashGame(PayoffVector row_player_payoff, PayoffVector col_player_payoff){

        const Eigen::MatrixXd payoff_matrix_1 = this->convert_vvd_to_matrix(row_player_payoff);
        const Eigen::MatrixXd payoff_matrix_2 = this->convert_vvd_to_matrix(col_player_payoff);

        this->RowPlayer    = payoff_matrix_1;
        this->ColumnPlayer = payoff_matrix_2;

        std::string path_leader = "payoff_leader.csv";
        std::ofstream payoff_leader(path_leader);

        std::string path_follower = "payoff_follower.csv";
        std::ofstream payoff_follower(path_follower);

        Eigen::IOFormat CleanFmt(4, 0, ", ", "\n");

        payoff_leader << payoff_matrix_1.format(CleanFmt);
        payoff_follower << payoff_matrix_2.format(CleanFmt);

        std::cout << " First (Row) Player " << std::endl;
        std::cout << payoff_matrix_1 << std::endl;

        std::cout << " Second (Column) Player " << std::endl;
        std::cout << payoff_matrix_2 << std::endl;

}


NashGame::~NashGame() {
        // TODO Auto-generated destructor stub
}


Eigen::MatrixXd NashGame::convert_vvd_to_matrix(std::vector<std::vector<double>> vvd){

    std::size_t n_rows = vvd.size();
    std::size_t n_cols = vvd.at(0).size();

    Eigen::MatrixXd result(n_rows, n_cols);
    result.row(0) = Eigen::VectorXd::Map(&vvd[0][0], n_cols);

    for (std::size_t i = 1; i < n_rows; i++) {
        result.row(i) = Eigen::VectorXd::Map(&vvd[i][0], n_cols);
    }

    return result;
}




Eigen::MatrixXd NashGame::make_tableau(Eigen::MatrixXd &M) {
    //==================//
    int m = M.rows();
    int n = M.cols();
    //==================//

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m,m);
    Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(m, 1);
    Eigen::MatrixXd C(m, m+n+1);

    C<< M,I, ones;

    return C;

}




Eigen::Index NashGame::find_pivot_row(const Eigen::MatrixXd& tableau, int column_index) {
    Eigen::Index id;
    Eigen::VectorXd ratios = tableau.col(column_index).array() / tableau.col(tableau.cols() - 1).array();
    ratios.array().maxCoeff(&id);
    return id;
}



std::vector<int> NashGame::non_basic_variables(const Eigen::MatrixXd &tableau) {

    Eigen::MatrixXd columns = tableau.block(0, 0, tableau.rows(), tableau.cols() - 1).transpose();
    std::vector<int> non_basic_vars;

    for (int i = 0; i < columns.rows(); i++) {
        Eigen::VectorXd vec_row = columns.row(i).array();
        Eigen::SparseMatrix<double> sparse_vec(vec_row.sparseView());
        if (sparse_vec.nonZeros() != 1) {
            non_basic_vars.push_back(i);
        }
    }
    return non_basic_vars;
}



std::vector<int> NashGame::pivot_tableau(Eigen::MatrixXd &tableau, int column_index) {

        std::vector<int> original_labels = non_basic_variables(tableau);
    const int pivot_row_index = find_pivot_row(tableau, column_index);
    double pivot_element = tableau(pivot_row_index, column_index);

    for (int i = 0; i < tableau.rows(); i++) {
        if (i != pivot_row_index) {
            tableau.row(i) = tableau.row(i) * pivot_element - tableau.row(pivot_row_index) * tableau(i, column_index);
        }
    }

    std::vector<int> non_basic_vars = this->non_basic_variables(tableau);


    for(auto& org_var: original_labels){
        std::vector<int>::iterator position = std::find(non_basic_vars.begin(), non_basic_vars.end(), org_var);
        if (position != non_basic_vars.end())
                non_basic_vars.erase(position);
    }

    return non_basic_vars;
}



Eigen::MatrixXd NashGame::shift_tableau(Eigen::MatrixXd tableau, int num_rows, int num_cols){
    Eigen::MatrixXd shifted_tableau(num_rows, num_cols);

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols-1; j++) {
            shifted_tableau(i, j) = tableau(i, j % num_cols);
        }
        shifted_tableau(i, num_cols-1) = 1.0;
    }

    return shifted_tableau;
}



Eigen::VectorXd NashGame::tableau_to_strategy(Eigen::MatrixXd tableau, std::vector<int> basic_labels, const int strategy_labels){
        std::vector<double> vertex;
    for (int column=0; column < strategy_labels; column++) {
        if (std::find(basic_labels.begin(), basic_labels.end(),column)!=basic_labels.end()) {
            for (int i = 0; i < tableau.rows(); i++) {
                if (tableau(i, column) != 0) {
                    vertex.push_back(tableau(i, tableau.cols()-1) / tableau(i, column));
                }
            }
        } else {
            vertex.push_back(0);
        }
    }


    double *ptr_vertex = &vertex[0];
    Eigen::Map<Eigen::VectorXd> strategy(ptr_vertex, vertex.size());

    double sum = strategy.sum();

    std::cerr << "	=	=	=	=	=	=	=	=	=	=	=	" << std::endl;
    for(auto& nn_bs_var: basic_labels)
        std::cerr << "nn_bsc_var  : " << nn_bs_var << std::endl;

    std::cerr << "strategy_labels  : " << strategy_labels << std::endl;
    std::cerr << strategy / sum << std::endl;
    std::cerr << "	=	=	=	=	=	=	=	=	=	=	=	" << std::endl;

    return strategy / sum;
}




// 0
void NashGame::combination(int n, int r, std::vector<std::vector<double>> &powerset){

        std::vector<double> v(n);
        std::fill(v.begin(), v.begin() + r, true);

        do {
                std::vector<double> elements;
                for (int i = 0; i < n; ++i) {
                   if (v[i]) {
                           elements.push_back(i);
                   }
                }
                if(elements.size() != 0)																//	(excluding the empty set)
                        powerset.push_back(elements);
        } while (std::prev_permutation(v.begin(), v.end()));

}


// 1
std::vector<std::vector<double>> NashGame::powerset(int n){

        std::vector<std::vector<double>> powerset;

        #pragma omp parallel for													// parallelize
        for(int i=0; i <= n; i++)
                this->combination(n, i, powerset);

        return powerset;

}




// 2
bool NashGame::solve_indifference(const PayoffMatrix& A, Eigen::VectorXd &prob, const std::vector<double> &rows, const std::vector<double> &columns){

    const int m = A.rows();
    const int n = A.cols();
    const int rows_m = rows.size();
    Eigen::MatrixXd M(m, n);
    std::vector<double> rows_rotated(rows);

    if (m != 0) {
                std::rotate(rows_rotated.begin(), rows_rotated.begin() + (rows_rotated.size()- 1), rows_rotated.end());				//	Vector is shifted forward by one

                for (int i = 0 ; i < m-1; ++i) {
                        M.row(i) = A.row(rows[i]) - A.row(rows_rotated[i]);
                }
    }

    Eigen::MatrixXd Z(1, n);
    Z.row(0).setOnes();

    for (int j : columns)
        Z(0, j) = 0;
    M.row(m-1) = Z;

    //	Set the last row to ones
    Z.row(0).setOnes();
    Eigen::MatrixXd M_new(M.rows()+1, M.cols());
    M_new << M, Z;

    Eigen::VectorXd b(m+1);
    b.setZero();
    b(m) = 1;

    Eigen::FullPivLU<Eigen::MatrixXd> lu(M_new);
    prob = lu.solve(b);

    if((prob.array() >= 0).all()){
        return true;
    }

    return false;
}




// 3
SupportPairs NashGame::potential_support_pairs(bool non_degenerate){

        int p1_num_strategies = this->RowPlayer.rows();
    int p2_num_strategies = this->RowPlayer.cols();
    auto p1_supports = this->powerset(p1_num_strategies);

    std::vector<std::pair<std::vector<double>, std::vector<double>>> result;
        #pragma omp parallel for													// parallelize
    for (const auto& support1 : p1_supports) {
        if (support1.empty()) continue;
        auto p2_supports = powerset(p2_num_strategies);
        for (const auto& support2 : p2_supports) {
            if (support2.empty()) {
                continue;
            }
            result.emplace_back(support1, support2);
        }
    }
    return result;
}



bool NashGame::obey_support(Eigen::VectorXd strategy_prob, std::vector<double> support_vec){

        double *ptr_support = &support_vec[0];
    Eigen::Map<Eigen::VectorXd> support(ptr_support, support_vec.size());

        if (strategy_prob.size() == 0) {
        return false;
    }
    for (int i = 0; i < strategy_prob.size(); ++i) {
        if ((support.array() == i).any() && strategy_prob(i) <= this->tol) {
            return false;
        }
        if ((support.array() != i).all() && strategy_prob(i) > this->tol) {
            return false;
        }
    }
    return true;
}




ProbabilityVectors NashGame::indifference_strategies(){

        ProbabilityVectors probVector;

        int tolerance = std::min(this->tol, 0.);

        SupportPairs pairs = this->potential_support_pairs(false);

        for(auto& pair: pairs){

            Eigen::VectorXd prob1, prob2;

            bool res1 = this->solve_indifference(this->ColumnPlayer.transpose(), prob1, pair.second, pair.first);
            bool res2 = this->solve_indifference(this->RowPlayer, prob2, pair.first, pair.second);

             if(res1 & res2){
                 bool os_1 = this->obey_support(prob1, pair.first);
                 bool os_2 = this->obey_support(prob2, pair.second);

                 if(os_1 & os_2){
                         probVector.push_back({prob1, prob2, pair.first, pair.second});
                 }
             }
        }

        return probVector;
}



bool NashGame::is_Nash(ProbabilityVector pv){

        std::vector<double> row_support_indices, column_support_indices;

        for (const auto& double_element : std::get<2>(pv)) {
                row_support_indices.push_back(double_element);
        }

        for (const auto& double_element : std::get<3>(pv)) {
                column_support_indices.push_back(double_element);
        }

        auto u = std::get<1>(pv);
        Eigen::VectorXd row_payoffs = this->RowPlayer * u;

        auto v = std::get<0>(pv);
        Eigen::VectorXd column_payoffs = this->ColumnPlayer.transpose() * v;

        double* ptr_row_ids = &row_support_indices[0];
        Eigen::Map<Eigen::ArrayXd> row_sups(ptr_row_ids, row_support_indices.size());

        double* ptr_col_ids = &column_support_indices[0];
        Eigen::Map<Eigen::ArrayXd> col_sups(ptr_col_ids, column_support_indices.size());


        Eigen::VectorXd row_support_payoffs = row_payoffs(row_sups);
        Eigen::VectorXd column_support_payoffs = column_payoffs(col_sups);

        return (
                row_payoffs.maxCoeff() == row_support_payoffs.maxCoeff()
                && column_payoffs.maxCoeff() == column_support_payoffs.maxCoeff()
        );
}



void NashGame::support_enumeration(){

        int num_eq = 0;																											//	The number of equilibria in the matrix
        std::vector<std::pair<Eigen::VectorXd,Eigen::VectorXd>> nash_equilibrias;

        ProbabilityVectors probability_vectors = this->indifference_strategies();

        for(auto& prob_elmnt : probability_vectors){
                bool is_nash = this->is_Nash(prob_elmnt);
                //bool is_nash = true;
                if(is_nash){
                        nash_equilibrias.push_back({std::get<0>(prob_elmnt),std::get<1>(prob_elmnt)});
                        num_eq++;
                }
        }

        if(num_eq%2 == 0){
                std::cerr << "==========================================" << std::endl;
                std::cerr << "====== WARNING : Degenerate Game =========" << std::endl;
                std::cerr << "==========================================" << std::endl;
        }else{
                for(auto& eq : nash_equilibrias){
                        std::cerr << std::endl;
                        std::cerr << "============== NASH EQUILIBRIA ===========" << std::endl;
                        std::cerr << "==========================================" << std::endl;
                        std::cerr << "Row player    (1) : [" << eq.first.transpose()   << "]" << std::endl;
                        std::cerr << "Column player (2) : [" << eq.second.transpose() << "]" << std::endl;
                        std::cerr << "==========================================" << std::endl;
                }
        }


}






