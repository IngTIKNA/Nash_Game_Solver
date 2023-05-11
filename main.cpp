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

    /*
        std::vector<std::vector<double>> row_payoff;
    std::vector<std::vector<double>> col_payoff;

    std::vector<double> row_row_1{ 1, 1, -1 };
    std::vector<double> row_row_2{ 2,-1,  0 };
    row_payoff.push_back(row_row_1);
    row_payoff.push_back(row_row_2);

    std::vector<double> col_row_1{ 0.5, -1, 0.5 };
    std::vector<double> col_row_2{ -1, 3, 2 };
    col_payoff.push_back(col_row_1);
    col_payoff.push_back(col_row_2);
        */

        NashGame NG(row_payoff, col_payoff);

        NG.support_enumeration();




    return 0;

}















