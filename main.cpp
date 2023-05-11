#include <iostream>
#include "NashGame.h"


int main() {

    int row, col;

    std::cerr << "What is the number of rows in the matrix?" << std::endl;
    std::cin  >> row;
    std::cerr << "What is the number of columns in the matrix?" << std::endl;
    std::cin  >> col;

    //==================================================================================================//
    std::cerr << "Enter values for the payoff matrix of the row player (player 1): " << std::endl;
    std::vector<std::vector<double>> row_payoff( row , std::vector<double> (col, 0));
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
        	std::cerr << "row_payoff[" << i << "][" << j << "]" << std::endl;
        	std::cin >> row_payoff[i][j];
        }
    }
    //==================================================================================================//

    //==================================================================================================//
    std::cerr << "Enter values for the payoff matrix of the column player (player 2): " << std::endl;
    std::vector<std::vector<double>> col_payoff( row , std::vector<double> (col, 0));
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
        	std::cerr << "col_payoff[" << i << "][" << j << "]" << std::endl;
        	std::cin >> col_payoff[i][j];
        }
    }
    //==================================================================================================//

	NashGame QV;

	const Eigen::MatrixXd payoff_matrix_1 = QV.convert_vvd_to_matrix(row_payoff);
	const Eigen::MatrixXd payoff_matrix_2 = QV.convert_vvd_to_matrix(col_payoff);

	NashGame NG(payoff_matrix_1, payoff_matrix_2);

    Eigen::MatrixXd Q(3,4);
    Q << 3,4,1,4,
         9,4,2,2,
         12,8,1,2;


    NG.support_enumeration();


    return 0;

}















