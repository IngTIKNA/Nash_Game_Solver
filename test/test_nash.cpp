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
//												SCENARIO 01															//
//__________________________________________________________________________________________________________________//

SCENARIO(" A bimatrix(3x3) game (Non-degenerate) Nash-Game ")
{
 
    GIVEN("3x3 bimatrix game - Find the powerset") {

        Eigen::MatrixXd first_payoff(3,3);
        first_payoff << 3.1,3.1,3.1,
                        3,3,3,
                        -2.9,2.9,2.9;

        Eigen::MatrixXd second_payoff(3,3);
        second_payoff<< 3,4,1,
                        3,2,1,
                        3,2,1;

        NashGame NG(first_payoff, second_payoff);

        int p1_num_strategies = first_payoff.rows();
        std::vector<std::vector<double>> powerset = NG.powerset(p1_num_strategies);

        std::cerr << "Powerset size : " <<  powerset.size() << std::endl;

        //for(auto& set: powerset){
        //    for(auto& elm: set)
        //        std::cerr <<  elm  << " - " << std::endl;    
        //    std::cerr << std::endl;
        //}

        //============================= TEST CASE 1 ================================//
		//==================================================================//		//
		WHEN("Number of elements in the powerset(excluding the empty set)"){//	    // 
																			//		//
		REQUIRE(powerset.size() == 7);     									//		//	assert
																			//		//
		}																	//		//
		//==================================================================//		//
		//==========================================================================//



        Eigen::VectorXd prob;        
        std::vector<double> rows_1;
        std::vector<double> columns_1;

        rows_1.push_back(0);
        columns_1.push_back(0);

        bool res = NG.solve_indifference(first_payoff, prob, rows_1, columns_1);

        std::cerr << "prob : " << std::endl;
        std::cerr << prob(0,0) << std::endl;
        std::cerr << prob(1,0) << std::endl;
        std::cerr << prob(2,0) << std::endl;

        //============================= TEST CASE 2 ================================//
		//==================================================================//		//
		WHEN(" Verify prob vector"){                                        //	    // 
																			//		//
		REQUIRE(prob.size() == 3);     									    //		//	assert
        REQUIRE(prob(0,0) == 1);     									    //		//	assert
        REQUIRE(prob(1,0) == 0);     									    //		//	assert
        REQUIRE(prob(2,0) == 0);     									    //		//	assert
																			//		//
		}																	//		//
		//==================================================================//		//
		//==========================================================================//










    }
        
}


