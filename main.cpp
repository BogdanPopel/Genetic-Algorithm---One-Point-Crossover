#include <iostream>
#include <bits/stdc++.h>
#include <random>
using namespace std;

ifstream in("data.txt");
ofstream out("out.txt");

class Chromosome{
private:
    int chromosome_length;
    pair<int, int> domain;
    int precision;
    double selection_probability;
    vector<bool> genes;
public:
    Chromosome(pair<int, int> domain, int precision){
        this->domain = domain;
        this->precision = precision;
        this->chromosome_length =(int) round(abs(log2( (domain.second - domain.first)*pow(10, precision) )));
        //this->genes.resize(chromosome_length);
        //generam genele daca nu sunt date
        for(int i = 0; i<this->chromosome_length;i++)
        {
            this->genes.push_back(rand()%2);
        }
    }

    vector<bool>& get_genome()
    {
        return this->genes;
    }
    int get_len(){return this->chromosome_length;}
    void set_selection_probability(double x){this->selection_probability = x;}
    double get_selection_probability(){return this->selection_probability;}
    void set_genome(vector<bool> v){this->genes.clear(); this->genes.assign(v.begin(), v.end());}
    double get_value(){
        double to_add =  (this->domain.second - this->domain.first) / (pow(2,this->chromosome_length)-1) ;
        double to_return = 0;
        for(int i = 0; i < this->chromosome_length; i++)
        {
            if(this->genes[i] == 1)
                to_return += to_add;

            to_add *= 2;
        }
        return to_return + this->domain.first;
    }

    string genome_toString(){
        string genome;
        for(auto g:genes)
        {
            if(g)
                genome += "1";
            else
                genome +="0";
        }
        return genome;
    }
};

class Mutation{
private:
    int population_size;
    int precision;
    int steps;
    double crossover_probability;
    double mutation_probability;
    double performance;
    pair<int,int> domain;
    int parameters[3];
    vector<Chromosome> population;
    vector<double> selection_interval; //q
public:
    Mutation(int n, int d_first, int d_second, int param1, int param2, int param3, int precision, double cp, double mp, int s )
    {
        this->population_size = n;
        this->domain.first = d_first;
        this->domain.second = d_second;
        this->parameters[0] = param1;
        this->parameters[1] = param2;
        this->parameters[2] = param3;
        this->precision = precision;
        this->crossover_probability = cp;
        this->mutation_probability = mp;
        this->steps = s;
        for(int i = 0; i<this->population_size; i++)
            this->population.emplace_back(Chromosome(domain, precision));
        double t = 0;
        for(int i = 0; i<this->population_size; i++){
            t += fitness(this->population[i].get_value());
        }
        this->performance = t;
        double aux = 0;
        set_selection_pb();
        for(int i = 0; i < this->population_size; i++){
            aux+=this->population[i].get_selection_probability();
            this->selection_interval.push_back(aux);
        }
        this->selection_interval.emplace(this->selection_interval.begin(), 0);
        this->selection_interval.emplace(this->selection_interval.end(), 1);
        sort(this->selection_interval.begin(), this->selection_interval.end());

    }

    int get_steps(){return this->steps;}
    double get_performance(){return this->performance;}

    double fitness(double x){
        return pow(x, 2)*this->parameters[0] + x*this->parameters[1] + this->parameters[2];
    }

    void set_selection_pb(){
        for(int i = 0; i < this->population_size; i++){
            double probability = fitness(this->population[i].get_value()) / get_performance();
            this->population[i].set_selection_probability(probability);
        }
    }
    //selectia proportionala, metoda ruletei
    vector<Chromosome> get_selected_population(){
        vector<Chromosome> new_one;
        for(int i = 0; i < this->population_size-1; i++)
        {
            double u = random_u();
            int chromosome = search_chr_on_selection_interval(u);
            new_one.emplace_back(this->population[chromosome]);
            out<<"u= "<<u<<"; now select chromosome "<<chromosome+1<<'\n';
        }
        return new_one;
    }
    double random_u(){
        random_device rd;
        default_random_engine eng(rd());
        uniform_real_distribution<double> distr(0, 1);
        eng.seed(rand()%1001);
        return distr(eng);
    }
    int search_chr_on_selection_interval(double u){
        int left = 1;
        int right = this->population_size-1;
        while(left<=right){
            int m = (left+right)/2;
            if(this->selection_interval[m+1] >= u && this->selection_interval[m-1]<=u)
                return m;
            else
                if(this->selection_interval[m] > u)
                    right = m-1;
                else
                    left = m+1;
        }
        return left-1;
    }
    Chromosome get_best_chr(){
        Chromosome best = population[0];
        double max = -99999;
        for(auto& c : population){
            double fit = fitness(c.get_value());
            if(fit > max)
            {
                max = fit;
                best = c;
            }
        }
        return best;
    }
    vector<Chromosome> cross_chromosomes(Chromosome& c1, Chromosome& c2, int cut)
    {
        vector<Chromosome> result;
        vector<bool> c1_genome = c1.get_genome();
        vector<bool> c2_genome = c2.get_genome();
        vector<bool> c1_crossed_genome = c1_genome;
        vector<bool> c2_crossed_genome = c2_genome;
        for(int i = 0; i<cut; i++){
            c1_crossed_genome[i] = c2_genome[i];
            c2_crossed_genome[i] = c1_genome[i];
        }
        c1.set_genome(c1_crossed_genome);
        c2.set_genome(c2_crossed_genome);
        result.emplace_back(c1);
        result.emplace_back(c2);
        return result;
    }

    void Start_Mutation(){
        out<<"Initial population:\n";
        for(int i = 0; i<this->population_size; i++)
        {
            out<<i+1<<" -> "<<this->population[i].genome_toString() <<" x= "<<this->population[i].get_value() << " f(x)= "<< fitness(population[i].get_value()) << '\n';
        }
        //print sel prob

        out<<"Selection probability\n";
        for(int i = 0; i < this->population_size; i++)
        {
            out<<i+1<<" probability: "<<this->population[i].get_selection_probability()<<'\n';
        }
        //print sel interval
        out<<"Selection interval:\n";
        for(int i = 0; i< this->selection_interval.size(); i++)
        {
            out<<this->selection_interval[i]<<" - ";
        }
        out<<'\n';
        Chromosome best = get_best_chr();
        this->population = move(get_selected_population());
        this->population.emplace_back(best);
        out<<"After selection:\n";
        for(int i = 0; i<this->population_size; i++)
        {
            out<<i+1<<" -> "<<this->population[i].genome_toString() <<" x= "<<this->population[i].get_value() << "f(x)= "<< fitness(population[i].get_value()) << '\n';
        }
        out<<"Crossover population:\n";
        vector<int> crossover_population_id;
        for(int i = 0; i<this->population_size; i++){
            double u = random_u();
            out<<i+1<<" - "<<this->population[i].genome_toString()<<" u="<<u;
            if(u < this->crossover_probability){
                crossover_population_id.emplace_back(i);
                out<<" =>selected for crossover\n";
            }
            else
                out<<'\n';
        }

        for(int i = 1; i < crossover_population_id.size(); i+=2)
        {
            out<<"Crossing now:\n";
            out<< crossover_population_id[i-1] + 1<<" - "<<population[crossover_population_id[i-1]].genome_toString()<<" and ";
            out<<crossover_population_id[i] + 1<<" - "<<population[crossover_population_id[i]].genome_toString()<<" ";
            int cut = rand()%this->population[i-1].get_len();
            out<<"at point: "<< cut<<'\n';
            vector<Chromosome> crossed_chromosomes;
            crossed_chromosomes = cross_chromosomes(this->population[crossover_population_id[i-1]], this->population[crossover_population_id[i]], cut);
            this->population[crossover_population_id[i-1]] = crossed_chromosomes[0];
            this->population[crossover_population_id[i]] = crossed_chromosomes[1];
            out<<"After Cross: "<<crossover_population_id[i-1] + 1<<" - "<<population[crossover_population_id[i-1]].genome_toString()<<" and ";
            out<<crossover_population_id[i] + 1<<" - "<<population[crossover_population_id[i]].genome_toString()<<'\n';
        }
        out<<"Population after crossover:\n";
        for(int i = 0; i<this->population_size; i++)
        {
            out<<i+1<<" -> "<<this->population[i].genome_toString() <<" x= "<<this->population[i].get_value() << " f(x)= "<< fitness(population[i].get_value()) << '\n';
        }

        vector<int> mutated_population_id;
        out<<"\nMutation probability: "<<this->mutation_probability;
        out<<'\n';
        out<<"Modified Chromosomes:\n";
        for(int i = 0; i<this->population_size; i++)
        {
            double u = random_u();
            if(u < this->mutation_probability)
            {
                mutated_population_id.push_back(i);
                out<<i<<'\n';
            }
        }
        //mutatie rara
        if(mutated_population_id.size())
            {
                for(int i = 0; i<mutated_population_id.size(); i++)
                {
                    int poz_rand = rand() % this->population[mutated_population_id[i]].get_len();
                    vector<bool> genome = this->population[mutated_population_id[i]].get_genome();
                    genome[poz_rand] = !genome[poz_rand];
                    this->population[mutated_population_id[i]].set_genome(genome);
                }

            }
        out<<"Population after mutation:\n";
        for(int i = 0; i<this->population_size; i++)
        {
            out<<i+1<<" -> "<<this->population[i].genome_toString() <<" x= "<<this->population[i].get_value() << " f(x)= "<< fitness(population[i].get_value()) << '\n';
        }

        out<<"Max_performance:\n";
        Chromosome c = get_best_chr();
        out<<fitness(c.get_value())<<'\n';
        double total = 0;
        for(int i = 0; i<this->population.size(); i++)
        {
            total += fitness(this->population[i].get_value());
        }
        out<<"Average performance:\n";
        out<<total/this->population_size;
    }

};

int main()
{
    int n, d_first, d_second, param1, param2, param3, precision, s;
    double cp, mp;
    in>>n>>d_first>>d_second>>param1>>param2>>param3>>precision>>cp>>mp>>s;
    Mutation *mutation = new Mutation(n, d_first, d_second, param1, param2, param3, precision, cp, mp, s);
    mutation->Start_Mutation();
    return 0;
}
