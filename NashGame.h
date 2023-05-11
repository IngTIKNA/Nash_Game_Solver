/*
 * NashGame.h
 *
 *  Created on: Nov 28, 2022
 *      Author: Ahmet Tikna
 */

#ifndef NASHGAME_H_
#define NASHGAME_H_

#include <iostream>
#include <bits/stdc++.h>
#include <vector>
#include <tuple>
#include <algorithm>
#include <string>
#include <fstream>
#include <exception>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>


typedef std::vector<std::vector<double>> PayoffVector;
typedef Eigen::MatrixXd PayoffMatrix;
typedef std::vector<std::pair<std::vector<double>, std::vector<double>>> SupportPairs;
typedef std::tuple<Eigen::VectorXd, Eigen::VectorXd, std::vector<double>, std::vector<double> > ProbabilityVector;
typedef std::vector<ProbabilityVector> ProbabilityVectors;


class NashGame {
public:
	/* First Player **/
	PayoffMatrix RowPlayer;
	/* Second Player **/
	PayoffMatrix ColumnPlayer;

	/* Default Constructor **/
	NashGame();

	/* Parameterized Constructor **/
	NashGame(PayoffVector row_player_payoff, PayoffVector col_player_payoff);

	/* Destructor **/
	virtual ~NashGame();

	Eigen::MatrixXd convert_vvd_to_matrix(std::vector<std::vector<double>> vvd);

	Eigen::MatrixXd make_tableau(Eigen::MatrixXd &M);

	Eigen::Index find_pivot_row(const Eigen::MatrixXd& tableau, int column_index);

	std::vector<int> non_basic_variables(const Eigen::MatrixXd &tableau);

	std::vector<int> pivot_tableau(Eigen::MatrixXd &tableau, int column_index);

	Eigen::MatrixXd shift_tableau(Eigen::MatrixXd tableau, int num_rows, int num_cols);

	Eigen::VectorXd tableau_to_strategy(Eigen::MatrixXd tableau, std::vector<int> basic_labels, const int strategy_labels);

	/* This function yields a set of combinations of n objects taken r **/
	void combination(int n, int r, std::vector<std::vector<double>> &powerset);

	/* This function finds the all possible combinations **/
	std::vector<std::vector<double>> powerset(int n);

	/* This function computes the indifference equations **/
	bool solve_indifference(const PayoffMatrix& A, Eigen::VectorXd &prob, const std::vector<double> &rows, const std::vector<double> &columns);

	/* This function returns the support pairs **/
	SupportPairs potential_support_pairs(bool non_degenerate = false);

	/* This function determines if the strategy has the given support **/
	bool obey_support(Eigen::VectorXd strategy_prob, std::vector<double> support_vec);

	/* This function produces probability vectors**/
	ProbabilityVectors indifference_strategies();

	/* This function evaluates a given probability vector to determine if it is Nash equilibrium **/
	bool is_Nash(ProbabilityVector pv);

	/* Run support enumeration method **/
	void support_enumeration();


private:
	/* Tolerance*/
	double tol = 1e-16;

};

#endif /* NASHGAME_H_ */

