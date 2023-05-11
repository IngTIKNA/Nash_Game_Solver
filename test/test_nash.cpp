/*
 * test_nash.cpp
 *
 *  Created on: April 15, 2023
 *      Author: Ahmet Tikna
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "../NashGame.h"


//__________________________________________________________________________________________________________________//
//                                                SCENARIO 01                                                       //
//__________________________________________________________________________________________________________________//

SCENARIO(" A bimatrix(3x3) game (Non-degenerate) Nash-Game ")
{

    GIVEN("3x3 bimatrix game - Find the powerset") {

        //==========================================================//
        //                     Payoff Matrices                      //
        //==========================================================//
        std::vector<std::vector<double>> row_player_payoff;
        std::vector<std::vector<double>> col_player_payoff;

        std::vector<double> row_row_1{ 3.1, 3.1, 3.1 };
        std::vector<double> row_row_2{   3,   3,   3 };
        std::vector<double> row_row_3{-2.9, 2.9, 2.9 };
        row_player_payoff.push_back(row_row_1);
        row_player_payoff.push_back(row_row_2);
        row_player_payoff.push_back(row_row_3);

        std::vector<double> col_row_1{   3,   4,   1 };
        std::vector<double> col_row_2{   3,   2,   1 };
        std::vector<double> col_row_3{   3,   2,   1 };
        col_player_payoff.push_back(col_row_1);
        col_player_payoff.push_back(col_row_2);
        col_player_payoff.push_back(col_row_3);
        //==========================================================//


        NashGame NG(row_player_payoff, col_player_payoff);

        //==============================================================================//
        //                             Generate power sets                              //
        //==============================================================================//
        int p1_num_strategies = row_player_payoff.size();                               //	row size
        int p2_num_strategies = col_player_payoff.size();                               //	row size

        std::vector<std::vector<double>> powerset1 = NG.powerset(p1_num_strategies);
        std::vector<std::vector<double>> powerset2 = NG.powerset(p2_num_strategies);
        //==============================================================================//


        //============================= TEST CASE 1 ================================//
        //==================================================================//		//
        WHEN("Number of elements in the powerset(excluding the empty set)") //      //
        {
        REQUIRE(powerset1.size() == 7);
        REQUIRE(powerset2.size() == 7);
        }                                                                   //      //
        //==================================================================//      //
        //==========================================================================//

        Eigen::VectorXd prob;
        std::vector<double> rows_1;
        std::vector<double> columns_1;

        rows_1.push_back(0);
        columns_1.push_back(0);

        bool res = NG.solve_indifference(NG.RowPlayer, prob, rows_1, columns_1);

        //============================= TEST CASE 2 ================================//
        //==================================================================//      //
        WHEN(" Verify prob vector"){

        REQUIRE(prob.size() == 3);
        REQUIRE(prob(0,0) == 1);
        REQUIRE(prob(1,0) == 0);
        REQUIRE(prob(2,0) == 0);

        }
        //==================================================================//      //
        //==========================================================================//



    }

}


