# Nash_Game_Solver

In noncooperative game theory, a game is defined by a set of strategies for each player and a payoff for each strategy.
There is a general solution technique for such games, known as Nash equilibrium, in which every strategy is a best response to the fixed strategies of other players. "Best response" case holds that all equilibrium strategies must lead to the maximal. There are algorithms for computing Nash equilibria between two players given either in strategic form or in extended form.

This solver uses the support enumeration method for computing the Nash Equilibrium in a 2-player bimatrix game.


## Installation

Nash Game Solver requires [Eigen](https://github.com/libigl/eigen) to run.

Install Nash Game Solver dependencies.

```sh
chmod +x install.sh
./install.sh

```

## References

- [Algorthmic Game Theory](https://www.cs.cmu.edu/~sandholm/cs15-892F13/algorithmic-game-theory.pdf) - Noam Nisan
- [Nashpy](https://github.com/drvinceknight/Nashpy)