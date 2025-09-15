#include <iostream>

#include "argparse.hpp"
#include "INIReader.h"
#include "abm.h"


int main(int argc, char* argv[]) {
    argparse::ArgumentParser main_program("abm");

    main_program.add_description("Agent Based Modelling");
    main_program.add_argument("--config")
        .required()
        .help(".ini Configuration fileuration file. Example provided in docs/");

    try {
        main_program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << main_program;
        return 1;
    }
    std::string config_file = main_program.get<std::string>("--config");
    INIReader reader(config_file);
    int error_code = reader.ParseError();
    if (error_code != 0) {
        std::cout << "Invalid config file. Error in reading " + config_file + "\n";
        if (error_code == -1) {
            std::cout << "File open error" << std::endl;
        } else if (error_code == -2) {
            std::cout << "Out of memory error" << std::endl;
        } else {
            std::cout << "Error line in config file is " + std::to_string(error_code)  << "\n";
        }
        return 1;
    }

    /* main_program.add_argument("--edgelist") */
    /*     .required() */
    /*     .help("Input edgelist (source, target)"); */
    /* main_program.add_argument("--nodelist") */
    /*     .required() */
    /*     .help("Input nodelist (node, year)"); */
    /* main_program.add_argument("--out-degree-bag") */
    /*     .required() */
    /*     .help("Input out-degree bag (, out-degree)"); */
    /* main_program.add_argument("--recency-probabilities") */
    /*     .required() */
    /*     .help("Input recency bag (year, probability)"); */
    /* main_program.add_argument("--planted-nodes") */
    /*     .default_value("") */
    /*     .help("Planted nodes file (year, fitness lag duration, fitness peak value, fitess peak duration, count)"); */
    /* main_program.add_argument("--preferential-weight") */
    /*     .default_value(double(-1)) */
    /*     .help("Preferential attachment weight") */
    /*     .scan<'g', double>(); */
    /* main_program.add_argument("--fitness-weight") */
    /*     .default_value(double(-1)) */
    /*     .help("Fitness weight") */
    /*     .scan<'g', double>(); */
    /* main_program.add_argument("--minimum-preferential-weight") */
    /*     .default_value(double(-1)) */
    /*     .help("Minimum preferential attachment weight") */
    /*     .scan<'g', double>(); */
    /* main_program.add_argument("--minimum-fitness-weight") */
    /*     .default_value(double(-1)) */
    /*     .help("Minimum fitness weight") */
    /*     .scan<'g', double>(); */
    /* main_program.add_argument("--fully-random-citations") */
    /*     .default_value(double(0.05)) */
    /*     .help("Constant percentage for radom citations") */
    /*     .scan<'g', double>(); */
    /* main_program.add_argument("--growth-rate") */
    /*     .required() */
    /*     .help("Growth rate") */
    /*     .scan<'g', double>(); */
    /* main_program.add_argument("--num-cycles") */
    /*     .required() */
    /*     .help("Number of years") */
    /*     .scan<'d', int>(); */
    /* main_program.add_argument("--same-year-proportion") */
    /*     .required() */
    /*     .help("Growth rate") */
    /*     .scan<'g', double>(); */
    /* main_program.add_argument("--neighborhood-sample") */
    /*     .default_value(int(-1)) */
    /*     .help("Number to sample from the 1 and 2 hop neighborhoods") */
    /*     .scan<'d', int>(); */
    /* main_program.add_argument("--output-file") */
    /*     .required() */
    /*     .help("Output clustering file"); */
    /* main_program.add_argument("--auxiliary-information-file") */
    /*     .required() */
    /*     .help("Auxillary information file"); */
    /* main_program.add_argument("--log-file") */
    /*     .required() */
    /*     .help("Output log file"); */
    /* main_program.add_argument("--num-processors") */
    /*     .default_value(int(1)) */
    /*     .help("Number of processors") */
    /*     .scan<'d', int>(); */
    /* main_program.add_argument("--log-level") */
    /*     .default_value(int(1)) */
    /*     .help("Log level where 0 = silent, 1 = info, 2 = verbose") */
    /*     .scan<'d', int>(); */

/* [Environment] */
/* nodelist=<FILE> ; csv with (node_id, publication_year) on each line */
/* edgelist=<FILE> ; csv with (source,target) on each line */
/* recency_table=<FILE> ; csv derived from a real-world network */
/* out_degree_bag=<FILE> ; csv derived from a real-world network */
/* growth_rate=<DOUBLE> ; floating point value e.g., 0.03 for 3% */
/* num_cycles=<INT> ; integer value e.g., 30 for a 30-year simulation */

/* [Agent] */
/* alpha=<FLOAT> ; floating point value specifying the alpha for neighborhood */
/* recency_bins=<FILE> ; csv file for recency bins */
/* same_year=<DOUBLE> ; floating point value e.g., 0.12 for 12% */
/* fully_random_citations=<DOUBLE> ; floating point value e.g., 0.05 for 5% */
/* preferential_weight=<DOUBLE> ; floating point value e.g., 0.33 */
/* fitness_weight=<DOUBLE> ; floating point value e.g., 0.33 */


/* [General] */
/* output_file=<FILE> ; output csv edgelist */
/* auxiliary_information_file=<FILE> ; output auxiliary information file */
/* log_file=<FILE> ; output log file */
/* num_processors=<INT> ; integer valued maximum parallelism allowed */
/* log_level=<INT> ; 0, 1, and 2 for silent, info, and verbose */
/* reader.Get("user", "name", "UNKNOWN") << ", email=" */
    std::string edgelist = reader.Get("Environment", "edgelist", "");
    std::string nodelist = reader.Get("Environment", "nodelist", "");
    std::string out_degree_bag = reader.Get("Environment", "out_degree_bag", "");
    std::string recency_table = reader.Get("Environment", "recency_table", "");
    std::string planted_nodes = reader.Get("Environment", "planted_nodes", "");
    double growth_rate = reader.GetReal("Environment", "growth_rate", -42);
    int num_cycles = reader.GetInteger("Environment", "num_cycles", -42);
    double fully_random_citations = reader.GetReal("Agent", "fully_random_citations", -42);
    double preferential_weight = reader.GetReal("Agent", "preferential_weight", -42);
    double fitness_weight = reader.GetReal("Agent", "fitness_weight", -42);
    int fitness_value_min = reader.GetInteger("Agent", "fitness_value_min", -42);
    int fitness_value_max = reader.GetInteger("Agent", "fitness_value_max", -42);
    double minimum_preferential_weight = reader.GetReal("Agent", "minimum_preferential_weight", -42);
    double minimum_fitness_weight = reader.GetReal("Agent", "minimum_fitness_weight", -42);
    double same_year_citations = reader.GetReal("Agent", "same_year_citations", -42);
    int neighborhood_sample = reader.GetInteger("Agent", "neighborhood_sample", -42);
    double alpha = reader.GetReal("Agent", "alpha", -42);
    double minimum_alpha = reader.GetReal("Agent", "minimum_alpha", -42);
    std::string use_alpha_string = reader.Get("Agent", "use_alpha", "");
    bool use_alpha = false;
    if (use_alpha_string == "true") {
        use_alpha = true;
    } else if (use_alpha_string == "false") {
        use_alpha = false;
    } else {
        std::cerr << "Valid values for required flag use_alpha are 'true' or 'false'" << std::endl;
        return 1;
    }
    std::string start_from_checkpoint_string = reader.Get("Environment", "start_from_checkpoint", "");
    bool start_from_checkpoint = false;
    if (start_from_checkpoint_string == "true") {
        start_from_checkpoint = true;
    } else if (start_from_checkpoint_string == "false") {
        start_from_checkpoint = false;
    } else {
        std::cerr << "Valid values for required flag start_from_checkpoint are 'true' or 'false'" << std::endl;
        std::cerr << start_from_checkpoint_string + " is invalid" << std::endl;
        return 1;
    }
    std::string recency_bins = reader.Get("Environment", "recency_bins", "");
    std::string output_file = reader.Get("General", "output_file", "");
    std::string auxiliary_information_file = reader.Get("General", "auxiliary_information_file", "");
    std::string log_file = reader.Get("General", "log_file", "");
    int num_processors = reader.GetInteger("General", "num_processors", -42);
    int log_level = reader.GetInteger("General", "log_level", -41) - 1;
    ABM* abm = new ABM(edgelist, nodelist, out_degree_bag, recency_table, recency_bins, alpha, minimum_alpha, use_alpha, start_from_checkpoint, planted_nodes, fully_random_citations, preferential_weight, fitness_weight, fitness_value_min, fitness_value_max, minimum_preferential_weight, minimum_fitness_weight, growth_rate, num_cycles, same_year_citations, neighborhood_sample, output_file, auxiliary_information_file, log_file, num_processors, log_level);
    abm->main();
    delete abm;
}
