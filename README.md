# Genetic-Algorithm---One-Point-Crossover

Genetic algorithm aimed at identifying the maximum of a positive function within a specific domain, with the predetermined function included in the code: 

double fitness(double x){
        return pow(x, 2)*this->parameters[0] + x*this->parameters[1] + this->parameters[2];
    } 

## Input data:

- Population Size: This represents the initial number of chromosomes to be used in our algorithm.

- Function Domain: This is the range of values in which the function will be evaluated to find the maximum.

- Parameters for the Maximization Function: The coefficients of the second-degree polynomial that define the objective function.

- Precision: The number of significant digits used to discretize the domain interval.

- Crossover Probability: The probability that two chromosomes will undergo genetic crossover in each iteration.

- Mutation Probability: The probability that a chromosome will undergo mutation in each iteration.

- Number of Algorithm Steps: The total number of iterations of the genetic algorithm.


## Output:

A descriptive text file, which meticulously documents the operations performed in the first stage of the algorithm. This file includes information such as the initial population, selection probabilities for each chromosome, cumulative probabilities defining selection intervals, highlighting the selection process based on generating a random number and binary search to determine the corresponding interval, chromosomes participating in recombination, randomly generated breakpoints, and the chromosomes resulting from recombination, as well as the population after recombination and the population after mutations.

Furthermore, for subsequent generations, the file will record only the maximum and average values of performance achieved. This detailed and transparent process helps us understand and evaluate the evolution of the genetic algorithm in the search for the maximum of the given function.

#### Example [output](https://github.com/BogdanPopel/Genetic-Algorithm---One-Point-Crossover/files/12715563/out.txt), generated using this [input](https://github.com/BogdanPopel/Genetic-Algorithm---One-Point-Crossover/files/12715570/data.txt).
 


