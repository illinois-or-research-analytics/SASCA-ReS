# SASCA(-ReS): Scalable Agent-based Simulator for Citation Analysis (with Recency-emphasized Sampling)
SASCA, or Scalable Agent-based Simulator for Citation Analysis, as the name suggests, is a scalable agent-based modeling simulator that can begin with a small seed network and simulate an exponential network growth to reach sizes of 100 million nodes and more. Currently, SASCA is implemented in modern C++ and can easily be parallelized across hundreds of cores.

## Dependencies
- C++ >= 20
- OpenMP >= 4.0
- cmake >= 3.23
- Eigen3 (can be locally installed via [setup.sh](setup.sh))
- PCG (can be locally installed via [setup.sh](setup.sh))

## One time setup
Run [setup.sh](setup.sh) to locally install [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and [PCG](https://www.pcg-random.org/). Alternatively, just ensure that both Eigen and PCG libraries are discoverable by cmake.

## How to build
SASCA is a standard cmake project. [easy_build_and_compile.sh](easy_build_and_compile.sh) is provided for user convenience.

## How to run
The general command is given below.
```console
abm --config <Configuration file>
```

An example configuration file, shown below, is also provided in the docs folder [Example configuration file](docs/example.ini)
Brief comments for each flag is listed inside the example configuration file. Detailed explanations for each of the flags are at the end of this section. All flags are required to be present in the configuration file. Unused flags can simply be left empty but the flag itself must still be present.


```
[Environment]
nodelist=<FILE> ; csv with (node_id, publication_year) on each line
edgelist=<FILE> ; csv with (source,target) on each line
recency_table=<FILE> ; csv derived from a real-world network
out_degree_bag=<FILE> ; csv derived from a real-world network
growth_rate=<DOUBLE> ; floating point value e.g., 0.03 for 3%
num_cycles=<INT> ; integer value e.g., 30 for a 30-year simulation
recency_bins=<STRING> ; string with comma separated values for each bin
planted_agents=<(optional) FILE> ; csv with planted agent in the format (year,fitness_lag_duration,fitness_peak_value,fitness_peak_duration,count) on each line
start_from_checkpoint=<BOOL> ; boolean value e.g., true or false for whether to start from a checkpoint or not

[Agent]
alpha=<FLOAT> ; floating point value specifying the alpha for neighborhood
use_alpha=<BOOL> ; boolean value e.g., true or false for whether to use alpha or not
same_year_citations=<DOUBLE> ; floating point value e.g., 0.12 for 12%
fully_random_citations=<DOUBLE> ; floating point value e.g., 0.05 for 5%
preferential_weight=<DOUBLE> ; floating point value e.g., 0.33
fitness_weight=<DOUBLE> ; floating point value e.g., 0.33
fitness_value_min=<INT> ; minimum fitness value (inclusive)
fitness_value_max=<INT> ; maximum fitness value (inclusive)
neighborhood_sample=<INT> ; maximum number of nodes to sample from each neighborhood


[General]
output_file=<FILE> ; output csv edgelist
auxiliary_information_file=<FILE> ; output auxiliary information file
log_file=<FILE> ; output log file
num_processors=<INT> ; integer valued maximum parallelism allowed
log_level=<INT> ; 0, 1, and 2 for silent, info, and verbose
```

For example, a simulation with no planted agents can be done using a configuration file structured as below.
```
...
recency_bins=1,5,10
planted_agents=
start_from_checkpoint=false
...
```

In some cases, such as when `use_alpha` is false, the value for `alpha` must be empty.

```
...
[Agent]
alpha=
use_alpha=false
...
```
### Individual flags
#### Environment flags
- `nodelist`: csv file with (node\_id, publication\_year) on each line. Both headers and their values are required. This is the starting seed network nodelist. When starting from a checkpoint, so when `start_from_checkpoint` is true, then `nodelist` supplied should be the auxiliary information file from the previous run with no modifications.
- `edgelist`: csv file with (source,target) on each line. This is the starting seed network edgelist. When starting from a checkpoint, so when `start_from_checkpoint` is true, then `edgelist` supplied should be the output file from the previous run with no modifications. Note that the header line for this CSV does not actually matter as long as the first column represents the source node and the second columns repreresns the target node. For example, the header line could be source,target or citing,cited, or u,v or any valid csv row.
- `recency_table` csv file derived from a real-world network. This file is usually created by taking a real-world citation network and getting the count Y which equals the number of edges that represent a citation across x years. This is encoded as a row in the csv such as x,y for each x. Note that for any simulation, a row for x needs to exist for all x from 1 to last year in the simulation - oldest publication year in the seed node set. Note that the header line for this CSV is required but the string values of the header line is not enforced. The only requirement is that starting from the second line of the csv that there are two comma separated integer values for each line.
- `out_degree_bag`: csv file derived from a real-world network. This file is usually created by taking real-world citation network and getting the out-degree of each node in the network. This is encoded as a row in the csv such as node\_id,out\_degree for each node in the network. Note that the node\_id column is unused but is useful for keeping track of which real-world publication was used. The headers in this file are not enforced except that two columns are required and that the second column must be integers.
- `growth_rate`: floating point value, e.g., 0.03 for 3%, which serves as the exponent for the exponential growth formula used in determining how many new agents should spawn per cycle of simulation.
- `num_cycles`: integer value e.g., 30 for a 30-year simulation
- `recency_bins`: This string is a user supplied comma separated string such that each comma separated value in the string is the start of the bin boundary. For example, a string "1,2,5,10,20" represents a binning for recency such that the first bin contains all publications that are less than 2 year old, meaning 1 <= current year - publication year < 2. It follows that the second bin now are publications that are at least 2 years old but at most 4 years old (2 <= current year - publication year < 5). The last bin is implied to go on until infinity so in this example, the last bin contains all publications that are at least 20 years old (20 <= current\_year - publication year < infinity).
- `planted_agents`: an optional csv file with planted agents. The required headers are (year,pa\_weight,fit\_weight,out\_degree,alpha,fitness\_lag\_duration,fitness\_peak\_value,fitness\_peak\_duration,count). Although the header line with all of the columns are required, only some values are required. For each line in the csv, values for year and count are required. Empty strings are interpreted as omittetd values. For these omitted values, the simulation default will be chosen. For example, in an experiment with random agents, a row string with "1,,,10,,,,,1" will be taken as planting a single agent in year 1 with out-degree 10. Preferential attachment weight, fitness weight, alpha, and fitness values will be whatever the simulation assigns this agent based on the input parameters. Note that each planted agent is chosen randomly from the given year and will override any assigned values. Planting too many agents will distort the distribution of affected values.
- `start_from_checkpoint` boolean value e.g., true or false for whether to start from a checkpoint or not. Check notes about `nodelist` and `edgelist` if this flag is set to true.

#### Agent flags
- `alpha` floating point value specifying the alpha for neighborhood. This value for alpha determines the proportion of citations that are made to the 1-hop nodes of the generator node relative to the total number of citations that the agent will make to the generator node's neighborhood. It can be left to be -1 for a random model or a constant value such as 0.5.
- `use_alpha` boolean value e.g., true or false for whether to use alpha or not. When alpha i
- `same_year_citations`: floating point value e.g., 0.12 for 12%. This value determines the proportion of new agents in a given year that cite another agent from the same year.
- `fully_random_citations`: floating point value e.g., 0.05 for 5%. This value determines the proportion of assigned out-degree that goes towards citing any node from the graph in a fully random manner.
- `preferential_weight`: floating point value e.g., 0.33. This value can be left -1 for agents with random weights or a constant value for static agents.
- `fitness_weight`: floating point value e.g., 0.33. This value can be left -1 for agents with random weights or a constant value for static agents.
- `fitness_value_min`: minimum possible integer fitness value for each agent.
- `fitness_value_max`: maximum possible integer fitness value for each agent.
- `neighborhood_sample`: maximum number of nodes to sample from each neighborhood. The neighborhoods of the genartor will be sampled such that the final list of nodes is at most this number.


#### General flags
- `output_file`: output csv edgelist file.
- `auxiliary_information_file`: output auxiliary information file a.k.a. the nodelist. This file contains the information for each node such as assigned out-degree, assigned alpha values, assigned fitness values, final in-degree, and so on.
- `log_file`: output log file. This log file contains information about the run itself such as the runtime. The configuration file flags will be parsed and the parsed results will be written to this file before the simulation begins. Additionally, an extra file is created in the same directory as the log file containing the wall-time of individual stages inside the simulator.
- `num_processors`: integer valued maximum parallelism allowed. Used for setting the maximum number of workers.
- `log_level`: enum flag where 0, 1, and 2 correspond to silent, info, and verbose, respectively.
