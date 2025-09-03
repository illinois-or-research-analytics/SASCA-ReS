#include "abm.h"
#include <iomanip>
#pragma omp declare reduction(merge_int_pair_vecs : std::vector<std::pair<int, int>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(merge_str_int_pair_vecs : std::vector<std::pair<std::string, int>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(merge_int_vecs : std::vector<int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

#pragma omp declare reduction(custom_merge_vec_int : std::vector<std::pair<int, int>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) initializer(omp_priv = decltype(omp_orig){})

int ABM::WriteToLogFile(std::string message, Log message_type) {
    if(this->log_level >= message_type) {
        std::chrono::time_point<std::chrono::steady_clock> now = std::chrono::steady_clock::now();
        std::string log_message_prefix;
        if(message_type == Log::info) {
            log_message_prefix = "[INFO]";
        } else if(message_type == Log::debug) {
            log_message_prefix = "[DEBUG]";
        } else if(message_type == Log::error) {
            log_message_prefix = "[ERROR]";
        }
        auto days_elapsed = std::chrono::duration_cast<std::chrono::days>(now - this->start_time);
        auto hours_elapsed = std::chrono::duration_cast<std::chrono::hours>(now - this->start_time - days_elapsed);
        auto minutes_elapsed = std::chrono::duration_cast<std::chrono::minutes>(now - this->start_time - days_elapsed - hours_elapsed);
        auto seconds_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - this->start_time - days_elapsed - hours_elapsed - minutes_elapsed);
        auto total_seconds_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - this->start_time);
        log_message_prefix += "[";
        log_message_prefix += std::to_string(days_elapsed.count());
        log_message_prefix += "-";
        log_message_prefix += std::to_string(hours_elapsed.count());
        log_message_prefix += ":";
        log_message_prefix += std::to_string(minutes_elapsed.count());
        log_message_prefix += ":";
        log_message_prefix += std::to_string(seconds_elapsed.count());
        log_message_prefix += "]";

        log_message_prefix += "(t=";
        log_message_prefix += std::to_string(total_seconds_elapsed.count());
        log_message_prefix += "s)";
        this->log_file_handle << log_message_prefix << " " << message << '\n';

        if(this->num_calls_to_log_write % 1 == 0) {
            std::flush(this->log_file_handle);
        }
        this->num_calls_to_log_write ++;
    }
    return 0;
}

void ABM::ReadPlantedNodes() {
    char delimiter = ',';
    std::ifstream planted_nodes_stream(this->planted_nodes);
    std::string line;
    int line_no = 1;
    while(std::getline(planted_nodes_stream, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        std::string year = current_line[0];
        std::string fitness_lag_duration = current_line[1];
        std::string fitness_peak_value = current_line[2];
        std::string fitness_peak_duration = current_line[3];
        std::string count = current_line[4];
        this->planted_nodes_map[std::stoi(year)][line_no]["fitness_lag_duration"] = std::stoi(fitness_lag_duration);
        this->planted_nodes_map[std::stoi(year)][line_no]["fitness_peak_value"] = std::stoi(fitness_peak_value);
        this->planted_nodes_map[std::stoi(year)][line_no]["fitness_peak_duration"] = std::stoi(fitness_peak_duration);
        this->planted_nodes_map[std::stoi(year)][line_no]["count"] = std::stoi(count);
        line_no += 1;
        /* std::cout << "adding " << current_line[1] << " to  bag " << std::endl; */
    }
}

void ABM::ReadOutDegreeBag() {
    char delimiter = ',';
    std::ifstream out_degree_bag_stream(this->out_degree_bag);
    // MARK: error
    std::string line;
    int line_no = 0;
    while(std::getline(out_degree_bag_stream, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        if(line_no != 0) {
            this->out_degree_bag_vec.push_back(std::stoi(current_line[1]));
        }
        line_no ++;
        /* std::cout << "adding " << current_line[1] << " to  bag " << std::endl; */
    }
}


std::unordered_map<int, double> ABM::GetBinnedRecencyProbabilities(Graph* graph, int current_year) {
    std::unordered_map<int, double> binned_recency_probabilities;
    double binned_recency_sum = 0;
    for(const auto& [year_diff, count] : this->recency_counts_map) {
        int current_bin_index = this->GetBinIndex(year_diff);
        binned_recency_probabilities[current_bin_index] += count;
        binned_recency_sum += count;
    }
    for(const auto& recency_pair : binned_recency_probabilities) {
        binned_recency_probabilities[recency_pair.first] /= binned_recency_sum;
    }
    return binned_recency_probabilities;
}


void ABM::ReadRecencyProbabilities() {
    char delimiter = ',';
    std::ifstream recency_probabilities_stream(this->recency_table);
    std::string line;
    int line_no = 0;
    while(std::getline(recency_probabilities_stream, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        if(line_no != 0) {
            int integer_year_diff = std::stoi(current_line[0]);
            /* double probability = std::stod(current_line[1]); */
            int count = std::stoi(current_line[1]);
            if(integer_year_diff > 0) {
                this->recency_counts_map[integer_year_diff] = count;
            }
        }
        line_no ++;
    }
}

std::unordered_map<int, int> ABM::BuildContinuousNodeMapping(Graph* graph) {
    int next_node_id = 0;
    std::unordered_map<int, int> continuous_node_mapping;
    for(auto const& node : graph->GetNodeSet()) {
        continuous_node_mapping[node] = next_node_id;
        next_node_id ++;
    }
    return continuous_node_mapping;
}

std::unordered_map<int, int> ABM::ReverseMapping(std::unordered_map<int, int> mapping) {
    std::unordered_map<int, int> reverse_mapping;
    for(auto const& [key,val] : mapping) {
        reverse_mapping[val] = key;
    }
    return reverse_mapping;
}


void ABM::FillInDegreeArr(Graph* graph, const std::unordered_map<int, int>& continuous_node_mapping, int* in_degree_arr) {
    for(auto const& node: graph->GetNodeSet()) {
        // MARK: debug
        if (!continuous_node_mapping.contains(node)) {
            std::cerr << "continuous mapping doesn't contain the node " + std::to_string(node) << std::endl;
            exit(1);
        }
        int continuous_id = continuous_node_mapping.at(node);
        in_degree_arr[continuous_id] = graph->GetInDegree(node);
    }
}

void ABM::InitializeFitness(Graph* graph) {
    this->AssignPeakFitnessValues(graph, graph->GetNodeSet());
    this->AssignFitnessLagDuration(graph, graph->GetNodeSet());
    this->AssignFitnessPeakDuration(graph, graph->GetNodeSet());
}

void ABM::FillFitnessArr(Graph* graph, const std::unordered_map<int, int>& continuous_node_mapping, int current_year, int* fitness_arr) {
    for(auto const& node : graph->GetNodeSet()) {
        int fitness_peak_value = graph->GetIntAttribute("fitness_peak_value", node);
        /* int fitness_lag_duration = graph->GetIntAttribute("fitness_lag_duration", node); */
        /* int fitness_peak_duration = graph->GetIntAttribute("fitness_peak_duration", node); */
        /* int published_year = graph->GetIntAttribute("year", node); */
        int continuous_index = continuous_node_mapping.at(node);
        /* if (published_year + fitness_lag_duration > current_year) { */
        /*     fitness_arr[continuous_index] = 1; */
        /* } else if (published_year + fitness_lag_duration + fitness_peak_duration >= current_year) { */
        fitness_arr[continuous_index] = fitness_peak_value;
        /* } else { */
        /*     double decayed_fitness_value = fitness_peak_value / pow(current_year - published_year - fitness_lag_duration - fitness_peak_duration + 1, this->fitness_decay_alpha); */
        /*     fitness_arr[continuous_index] = decayed_fitness_value; */
        /* } */
    }
}

int ABM::GetMaxYear(Graph* graph) {
    int max_year = -1;
    bool is_first = true;
    for(auto const& node : graph->GetNodeSet()) {
        int current_node_year = graph->GetIntAttribute("year", node);
        if (is_first) {
            max_year = current_node_year;
            is_first = false;
        }
        if (current_node_year > max_year) {
            max_year = current_node_year;
        }
    }
    return max_year;
}

int ABM::GetMaxNode(Graph* graph) {
    int max_node = -1;
    bool is_first = true;
    for(auto const& node : graph->GetNodeSet()) {
        if (is_first) {
            max_node = node;
            is_first = false;
        }
        if (node > max_node) {
            max_node = node;
        }
    }
    return max_node;
}

int ABM::GetFinalGraphSize(Graph* graph) {
    int current_graph_size = graph->GetNodeSet().size();
    for(int i = 0; i < this->num_cycles; i ++) {
        int num_new_nodes = std::ceil(current_graph_size * this->growth_rate);
        current_graph_size += num_new_nodes;
    }
    return current_graph_size;
}

void ABM::PopulateWeightArrs(double* pa_weight_arr, double* fit_weight_arr, int len) {
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    if(this->preferential_weight != -1 &&  this->fitness_weight != -1) {
        for(int i = 0; i < len; i ++) {
            double pa_uniform = this->preferential_weight;
            double fit_uniform = this->fitness_weight;
            double sum = pa_uniform + fit_uniform;
            pa_weight_arr[i] = (double)pa_uniform / sum;
            fit_weight_arr[i] = (double)fit_uniform / sum;
        }
    } else {
        // MARK:is this why?
        for(int i = 0; i < len; i ++) {
            std::uniform_real_distribution<double> first_weights_uniform_distribution{0, 1};
            double first_uniform = first_weights_uniform_distribution(generator);
            std::vector<double> current_weight_array{first_uniform, 1 - first_uniform};
            std::ranges::shuffle(current_weight_array, generator);
            std::ranges::shuffle(current_weight_array, generator);
            std::ranges::shuffle(current_weight_array, generator);
            for(size_t j = 0; j < current_weight_array.size(); j ++) {
                current_weight_array[j] = std::round(current_weight_array[j] * 1000.0) / 1000.0;
            }
            pa_weight_arr[i] = current_weight_array[0];
            fit_weight_arr[i] = current_weight_array[1];
        }
    }
}

void ABM::PopulateAlphaArr(double* alpha_arr, int len) {
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);

    if(this->alpha < 0) {
        for(int i = 0; i < len; i ++) {
            double alpha_uniform = this->alpha_uniform_distribution(generator);
            alpha_uniform = std::round(alpha_uniform * 1000.0) / 1000.0;
            alpha_arr[i] = alpha_uniform;
        }
    } else if(this->minimum_alpha > 0) {
        for(int i = 0; i < len; i ++) {
            std::uniform_real_distribution<double> minimum_alpha_uniform_distribution{minimum_alpha, 1};
            double alpha_uniform = minimum_alpha_uniform_distribution(generator);
            alpha_arr[i] = alpha_uniform;
        }
    } else {
        for(int i = 0; i < len; i ++) {
            alpha_arr[i] = this->alpha;
        }
    }
}

void ABM::PopulateOutDegreeArr(int* out_degree_arr, int len) {
    std::uniform_int_distribution<int> outdegree_index_uniform_distribution{0, (int)(this->out_degree_bag_vec.size() - 1)};
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    for(int i = 0; i < len; i ++) {
        int index_uniform = outdegree_index_uniform_distribution(generator);
        out_degree_arr[i] = this->out_degree_bag_vec[index_uniform];
    }
}

void ABM::UpdateGraphAttributesWeights(Graph* graph, int next_node_id, double* pa_weight_arr, double* fit_weight_arr, int len) {
    for(int i = 0; i < len; i ++) {
        int current_node_id = next_node_id + i;
        graph->SetDoubleAttribute("preferential_attachment_weight", current_node_id, pa_weight_arr[i]);
        graph->SetDoubleAttribute("fitness_weight", current_node_id, fit_weight_arr[i]);
    }
}

void ABM::UpdateGraphAttributesOutDegrees(Graph* graph, int next_node_id, int* out_degree_arr, int len) {
    for(int i = 0; i < len; i ++) {
        int current_node_id = next_node_id + i;
        graph->SetIntAttribute("assigned_out_degree", current_node_id, out_degree_arr[i]);
    }
}

std::vector<int> ABM::GetGraphAttributesGeneratorNodes(Graph* graph, int new_node) const {
    std::vector<int> generator_nodes;
    std::string generator_node_string = graph->GetStringAttribute("generator_node_string", new_node);
    std::stringstream ss(generator_node_string);
    std::string current_value;
    while(std::getline(ss, current_value, ';')) {
        generator_nodes.push_back(std::stoi(current_value));
    }
    return generator_nodes;
}

void ABM::UpdateGraphAttributesGeneratorNodes(Graph* graph, int new_node, const std::vector<int>& generator_nodes) {
    std::string generator_node_string;
    generator_node_string += std::to_string(generator_nodes.at(0));
    for(size_t i = 1; i < generator_nodes.size(); i ++) {
        generator_node_string += ";";
        generator_node_string += std::to_string(generator_nodes.at(i));
    }
    graph->SetStringAttribute("generator_node_string", new_node, generator_node_string);
}

void ABM::CalculateTanhScores(std::unordered_map<int, double>& cached_results, int* src_arr, double* dst_arr, int len) {
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < len; i ++) {
        double current_dst = -1;
        if (src_arr[i] < 10000) {
            current_dst = cached_results[src_arr[i]];
        } else {
            current_dst = this->peak_constant * std::tanh((pow(src_arr[i], 3)/this->delay_constant)*(1/this->peak_constant));
            /* current_dst = (src_arr[i] * this->gamma) + 1; */
        }
        dst_arr[i] = current_dst;
        sum += current_dst;
    }
    #pragma omp parallel for
    for(int i = 0; i < len; i ++) {
        dst_arr[i] /= sum;
    }
}

void ABM::CalculateExpScores(std::unordered_map<int, double>& cached_results, int* src_arr, double* dst_arr, int len) {
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for(int i = 0; i < len; i ++) {
        double current_dst = -1;
        /* if (src_arr[i] < 1000) { */
        /*     current_dst = cached_results[src_arr[i]]; */
        /* } else { */
        current_dst = std::max(pow(src_arr[i], this->gamma), 1.0) + 1;
            /* current_dst = (src_arr[i] * this->gamma) + 1; */
        /* } */
        dst_arr[i] = current_dst;
        sum += current_dst;
    }
    #pragma omp parallel for
    for(int i = 0; i < len; i ++) {
        dst_arr[i] /= sum;
    }
}

void ABM::FillSameYearSourceNodes(std::set<int>& same_year_source_nodes, int current_year_new_nodes) {
    size_t num_same_year_source_nodes = (size_t)std::floor(current_year_new_nodes * this->same_year_citations);
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    std::uniform_int_distribution<int> int_uniform_distribution(0, current_year_new_nodes - 1);
    while(same_year_source_nodes.size() != num_same_year_source_nodes) {
        int current_source = int_uniform_distribution(generator);
        if (same_year_source_nodes.count(current_source) == 0) {
            same_year_source_nodes.insert(current_source);
        }
    }
}

int ABM::MakeSameYearCitations(const std::set<int>& same_year_source_nodes, int num_new_nodes, const std::unordered_map<int, int>& reverse_continuous_node_mapping, int* citations, int current_graph_size) {
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    std::uniform_int_distribution<int> int_uniform_distribution(0, num_new_nodes - 1);
    int current_citation = int_uniform_distribution(generator);
    while(same_year_source_nodes.contains(current_citation)) {
        current_citation = int_uniform_distribution(generator);
    }
    citations[0] = reverse_continuous_node_mapping.at(current_graph_size + current_citation);
    // MARK: debug
    if (citations[0] == 0) {
        std::cerr << "same year target chosen as node id 0" << std::endl;
        std::cerr << "reverse continuous node mapping of " + std::to_string(current_graph_size) + " + " + std::to_string(current_citation) << std::endl;
        exit(1);
    }
    return 1;
}

int ABM::MakeUniformRandomCitationsFromGraph(Graph* graph, const std::unordered_map<int, int>& reverse_continuous_node_mapping, std::vector<int>& generator_nodes, int* citations, int num_cited_so_far, int num_citations) {
    if (num_citations <= 0) {
        return 0;
    }
    int actual_num_cited = num_citations;
    std::set<int> selected;
    for(int i = 0; i < num_cited_so_far; i ++) {
        selected.insert(citations[i]);
    }
    for(size_t i = 0; i < generator_nodes.size(); i ++) {
        selected.insert(generator_nodes.at(i));
    }
    if ((int)graph->GetNodeSet().size() - (int)selected.size() <= num_citations) {
        actual_num_cited = (int)graph->GetNodeSet().size() - (int)selected.size();
        for(int i = 0; i < actual_num_cited; i ++) {
            int current_cited_node = reverse_continuous_node_mapping.at(i);
            // MARK: debug
            if (current_cited_node == 0) {
                std::cerr << "reverse continuous_node_mapping has node id 0" << std::endl;
                exit(1);
            }
            if (!selected.contains(current_cited_node)) {
                citations[num_cited_so_far + i] = current_cited_node;
                selected.insert(current_cited_node);
            }
        }
    } else {
        pcg_extras::seed_seq_from<std::random_device> rand_dev;
        pcg32 generator(rand_dev);
        std::uniform_int_distribution<int> int_uniform_distribution(0, (int)(graph->GetNodeSet().size() - 1));
        int current_citation_index = 0;
        /* while(selected.size() != num_cited_so_far + generator_nodes.size() + actual_num_cited) { */
        while(current_citation_index < actual_num_cited) {
            int current_citation = int_uniform_distribution(generator);
            int current_cited_node = reverse_continuous_node_mapping.at(current_citation);
            // MARK: debug
            if (current_cited_node == 0) {
                std::cerr << "reverse continuous_node_mapping has node id 0" << std::endl;
                exit(1);
            }
            if (!selected.contains(current_cited_node)) {
                citations[num_cited_so_far + current_citation_index] = current_cited_node;
                selected.insert(current_cited_node);
                current_citation_index ++;
            }
        }
    }
    return actual_num_cited;
}

int ABM::MakeUniformRandomCitations(Graph* graph, const std::unordered_map<int, int>& continuous_node_mapping, int current_year, const std::vector<int>& candidate_nodes, int* citations, int current_graph_size, int num_citations) {
    if (num_citations <= 0) {
        return 0;
    }
    if (candidate_nodes.size() <= 0) {
        return 0;
    }

    int actual_num_cited = num_citations;
    if (candidate_nodes.size() <= (size_t)num_citations) {
        actual_num_cited = candidate_nodes.size();
        for(int i = 0; i < actual_num_cited; i ++) {
            citations[i] = candidate_nodes.at(i);
            // MARK: debug
            if (candidate_nodes.at(i) == 0) {
                std::cerr << "candidate nodes in last bin contains node id 0" << std::endl;
                exit(1);
            }
        }
    } else {
        pcg_extras::seed_seq_from<std::random_device> rand_dev;
        pcg32 generator(rand_dev);
        std::uniform_int_distribution<int> int_uniform_distribution(0, (int)(candidate_nodes.size() - 1));
        std::set<int> selected;
        int current_citation_index = 0;
        while((int)selected.size() != actual_num_cited) {
            int current_selected_index = int_uniform_distribution(generator);
            int current_node = candidate_nodes.at(current_selected_index);
            // MARK: debug
            if (current_node == 0) {
                std::cerr << "candidate nodes in last bin contains node id 0" << std::endl;
                exit(1);
            }
            if (!selected.contains(current_node)) {
                citations[current_citation_index] = current_node;
                selected.insert(current_node);
                current_citation_index ++;
            }
        }
    }
    return actual_num_cited;
}

int ABM::MakeCitations(Graph* graph, const std::unordered_map<int, int>& continuous_node_mapping, int current_year, const std::vector<int>& candidate_nodes, int* citations, double* pa_arr, double* fit_arr, double pa_weight, double fit_weight, int current_graph_size, int num_citations) {
    if (num_citations <= 0) {
        return 0;
    }
    if (candidate_nodes.size() <= 0) {
        return 0;
    }

    int actual_num_cited = num_citations;
    if (candidate_nodes.size() < (size_t)num_citations) {
        actual_num_cited = candidate_nodes.size();
    }
    std::vector<std::pair<double, int>> element_index_vec;
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    /* std::uniform_real_distribution<double> wrs_uniform_distribution{std::numeric_limits<double>::min(), 1}; */
    /*
    std::uniform_real_distribution<double> wrs_uniform_distribution{0, 1};
    double pa_sum = 0.0;
    double rec_sum = 0.0;
    double fit_sum = 0.0;
    for(size_t i = 0; i < candidate_nodes.size(); i ++) {
        int continuous_node_id = continuous_node_mapping.at(candidate_nodes.at(i));
        double current_pa = pa_arr[continuous_node_id];
        double current_rec = recency_arr[continuous_node_id];
        double current_fit = fit_arr[continuous_node_id];
        pa_sum += current_pa;
        rec_sum += current_rec;
        fit_sum += current_fit;
    }
    auto cmp = [](const std::pair<double, int> &left, const std::pair<double, int> &right) {
        return left.first > right.first;
    };
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, decltype(cmp)> min_heap(cmp);
    for (size_t i = 0; i < candidate_nodes.size(); i++) {
        int continuous_node_id = continuous_node_mapping.at(candidate_nodes.at(i));
        double current_pa = pa_arr[continuous_node_id] / pa_sum;
        double current_rec = recency_arr[continuous_node_id] / rec_sum;
        double current_fit = fit_arr[continuous_node_id] / fit_sum;

        double weighted_score = (pa_weight * current_pa) + (rec_weight * current_rec) + (fit_weight * current_fit);

        double base = wrs_uniform_distribution(generator);
        double wrs_score = std::pow(base, 1.0 / weighted_score);
        // double log_base = -wrs_exp_distribution(generator);
        // double wrs_score = log_base / weighted_score;

        if ((int)min_heap.size() < actual_num_cited) {
            min_heap.emplace(wrs_score, candidate_nodes.at(i));
        } else if (wrs_score > min_heap.top().first) {
            min_heap.pop();
            min_heap.emplace(wrs_score, candidate_nodes.at(i));
        }
    }

    for (int i = 0; i < actual_num_cited; i++) {
        citations[i] = min_heap.top().second;
        min_heap.pop();
    }
    */
    // /*
    // [pa score node i, fit score node i]   *   [ pa weight]
    // [fit score row vector]      [ fit weight]
    // node i final score = pa score * pa weight + fit score * fit weight
    Eigen::MatrixXd current_scores(candidate_nodes.size(), 2);
    Eigen::Vector2d current_weights(pa_weight, fit_weight);
    double pa_sum = 0.0;
    double fit_sum = 0.0;
    std::vector<double> raw_pa_arr;
    std::vector<double> raw_fit_arr;
    raw_pa_arr.reserve(candidate_nodes.size());
    raw_fit_arr.reserve(candidate_nodes.size());
    for(size_t i = 0; i < candidate_nodes.size(); i ++) {
        int continuous_node_id = continuous_node_mapping.at(candidate_nodes.at(i));
        double current_pa = pa_arr[continuous_node_id];
        double current_fit = fit_arr[continuous_node_id];
        pa_sum += current_pa;
        fit_sum += current_fit;
        raw_pa_arr.push_back(current_pa);
        raw_fit_arr.push_back(current_fit);
    }
    for(size_t i = 0; i < candidate_nodes.size(); i ++) {
        current_scores(i, 0) = raw_pa_arr[i] / pa_sum;
        current_scores(i, 1) = raw_fit_arr[i] / fit_sum;
    }
    Eigen::MatrixXd current_weighted_scores = current_scores * current_weights;
    auto current_wrs_uniform = [&] () {return wrs_uniform_distribution(generator);};
    Eigen::ArrayXd current_bases = Eigen::ArrayXd::NullaryExpr(candidate_nodes.size(), current_wrs_uniform);
    Eigen::ArrayXd weighted_random_sampling_results = current_bases.pow(1.0 / current_weighted_scores.array());


    for(size_t i = 0; i < candidate_nodes.size(); i ++) {
        if (candidate_nodes.at(i) == 0) {
            std::cerr << "0 is included in the candidate node set" << std::endl;
            exit(1);
        }
        /* if (candidate_nodes.at(i) >= 14694932) { */
        /*     std::cerr << "node id greater than or equal to 14694932 is included in the candidate node set" << std::endl; */
        /*     exit(1); */
        /* } */
        element_index_vec.push_back({weighted_random_sampling_results(i), candidate_nodes.at(i)});
    }
    std::ranges::shuffle(element_index_vec, generator);
    /* std::sort(element_index_vec.begin(), element_index_vec.end(), [](auto& left, auto& right){ */
    std::partial_sort(element_index_vec.begin(), element_index_vec.begin() + actual_num_cited, element_index_vec.end(), [](auto& left, auto& right){
        return left.first > right.first; // read
    });
    for (int i = 0; i < actual_num_cited; i ++) {
        citations[i] = element_index_vec[i].second;
        // MARK: debug
        if (element_index_vec[i].second == 0) {
            std::cerr << "trying to cite node id 0" << std::endl;
            exit(1);
        }
        /* if (candidate_nodes.at(i) >= 14694932) { */
        /*     std::cerr << "node id greater than or equal to 14694932 is included in the candidate node set" << std::endl; */
        /*     exit(1); */
        /* } */
    }
    //*/

    return actual_num_cited;
}

std::vector<int> ABM::GetGeneratorNodes(Graph* graph, const std::unordered_map<int, int>& reverse_continuous_node_mapping) {
    std::vector<int> generator_nodes;
    std::uniform_int_distribution<int> generator_uniform_distribution{0, (int)(graph->GetNodeSet().size() - 1)};
    int num_generator_nodes = 1;
    pcg_extras::seed_seq_from<std::random_device> rand_dev;
    pcg32 generator(rand_dev);
    for(int i = 0; i < num_generator_nodes; i ++) {
        int continuous_generator_node = generator_uniform_distribution(generator);
        int generator_node = reverse_continuous_node_mapping.at(continuous_generator_node);
        generator_nodes.push_back(generator_node);
    }
    return generator_nodes;
}

std::unordered_map<int, int> ABM::BinOutdegrees(const std::unordered_map<int, std::vector<int>>& binned_neighborhood, int total_outdegree, std::unordered_map<int, double> binned_recency_probabilities) {
    // initially let's say total_outdegree (requested) = 100
    // [0.5, 0.2, 0.2, 0.1] recency probability that's the same for all agents
    // split 100 into proportions
    // [50, 20, 20, 10]
    // first check sum(E) == 100
    // [(6,7,9), (10,9,2), (19,1), (5 * )] actual neighborhood with node ids
    // [50 - 3, 20 - 3, 20 - 3, 10 - 1] -> "unfulfilled quota"
    // [47, 17, 17, 9] -> "unfulfilled quota"
    // look left and right to see if any bins are available
    // end up citing [3, 3, 2, x] things from each bin
    std::unordered_map<int, int> target_outdegree_per_bin_map;
    int remaining_outdegree = total_outdegree;
    for(int bin_index = 0; bin_index < this->num_bins; bin_index ++) {
        if (remaining_outdegree == 0) {
            break;
        }
        double bin_probability = binned_recency_probabilities[bin_index];
        int current_bin_outdegree = std::round(total_outdegree * bin_probability);
        current_bin_outdegree = std::min(current_bin_outdegree, remaining_outdegree);
        target_outdegree_per_bin_map[bin_index] = current_bin_outdegree;
        remaining_outdegree -= current_bin_outdegree;
    }
    // MARK: debug
    /* for(int bin_index = 0; bin_index < remaining_outdegree; bin_index ++) { */
    /*     target_outdegree_per_bin_map[bin_index] += 1; */
    /* } */
    int current_remaining_bin_index = 0;
    while(remaining_outdegree > 0) {
        target_outdegree_per_bin_map[current_remaining_bin_index] += 1;
        current_remaining_bin_index ++;
        remaining_outdegree --;
        if(current_remaining_bin_index >= this->num_bins) {
            current_remaining_bin_index = 0;
        }
    }
    int outdegree_per_bin_sum = 0;
    for(const auto& current_bin : target_outdegree_per_bin_map) {
        int bin_index = current_bin.first;
        outdegree_per_bin_sum += target_outdegree_per_bin_map[bin_index];
    }
    if (outdegree_per_bin_sum != total_outdegree) {
        std::cerr << "first check: outdegree per bin sum doesn't match total outdegree" << std::endl;
        std::cerr << "total outdgree is " << std::to_string(total_outdegree) << std::endl;
        for(int bin_index = 0; bin_index < this->num_bins; bin_index ++) {
            std::cerr << "bin: " << std::to_string(bin_index) << " has " << std::to_string(target_outdegree_per_bin_map[bin_index]) << std::endl;
        }
        exit(-1);
    }
    for(int bin_index = 0; bin_index < this->num_bins; bin_index ++) {
        int current_uncited_num_nodes = target_outdegree_per_bin_map[bin_index] - binned_neighborhood.at(bin_index).size();
        if (current_uncited_num_nodes > 0) {
            for(int sweep_index = bin_index - 1; sweep_index >= 0 && current_uncited_num_nodes > 0; sweep_index --) {
                if (target_outdegree_per_bin_map[sweep_index] < (int)binned_neighborhood.at(sweep_index).size()) {
                    int current_citable = std::min(current_uncited_num_nodes, (int)binned_neighborhood.at(sweep_index).size() - target_outdegree_per_bin_map[sweep_index]);
                    current_uncited_num_nodes -= current_citable;
                    target_outdegree_per_bin_map[sweep_index] += current_citable;
                    target_outdegree_per_bin_map[bin_index] -= current_citable;
                }
            }
            for(int sweep_index = bin_index + 1; sweep_index < this->num_bins && current_uncited_num_nodes > 0; sweep_index ++) {
                if (target_outdegree_per_bin_map[sweep_index] < (int)binned_neighborhood.at(sweep_index).size()) {
                    int current_citable = std::min(current_uncited_num_nodes, (int)binned_neighborhood.at(sweep_index).size() - target_outdegree_per_bin_map[sweep_index]);
                    current_uncited_num_nodes -= current_citable;
                    target_outdegree_per_bin_map[sweep_index] += current_citable;
                    target_outdegree_per_bin_map[bin_index] -= current_citable;
                }
            }
        }
    }
    /*
    std::unordered_map<int, int> uncited_nodes;
    for(int bin_index = 0; bin_index < this->num_bins; bin_index ++) {
        uncited_nodes[bin_index] = target_outdegree_per_bin_map[bin_index] - binned_neighborhood.at(bin_index).size();
    }
    for(int bin_index = 0; bin_index < this->num_bins; bin_index ++) {
        int current_uncited_num_nodes = uncited_nodes[bin_index];
        if (current_uncited_num_nodes > 0) {
            for(int sweep_index = bin_index - 1; sweep_index >= 0 && uncited_nodes[bin_index] > 0; sweep_index --) {
                if (target_outdegree_per_bin_map[sweep_index] < (int)binned_neighborhood.at(sweep_index).size()) {
                    int current_citable = std::min(uncited_nodes[bin_index], (int)binned_neighborhood.at(sweep_index).size() - target_outdegree_per_bin_map[sweep_index]);
                    uncited_nodes[bin_index] -= current_citable;
                    target_outdegree_per_bin_map[sweep_index] += current_citable;
                    target_outdegree_per_bin_map[bin_index] -= current_citable;
                }
            }
            for(int sweep_index = bin_index + 1; sweep_index < this->num_bins && uncited_nodes[bin_index] > 0; sweep_index ++) {
                if (target_outdegree_per_bin_map[sweep_index] < (int)binned_neighborhood.at(sweep_index).size()) {
                    int current_citable = std::min(uncited_nodes[bin_index], (int)binned_neighborhood.at(sweep_index).size() - target_outdegree_per_bin_map[sweep_index]);
                    uncited_nodes[bin_index] -= current_citable;
                    target_outdegree_per_bin_map[sweep_index] += current_citable;
                    target_outdegree_per_bin_map[bin_index] -= current_citable;
                }
            }
        }
    }
    */
    outdegree_per_bin_sum = 0;
    int neighborhood_per_bin_sum = 0;
    for(const auto& current_bin : target_outdegree_per_bin_map) {
        int bin_index = current_bin.first;
        outdegree_per_bin_sum += target_outdegree_per_bin_map[bin_index];
        neighborhood_per_bin_sum += binned_neighborhood.at(bin_index).size();
    }
    if (outdegree_per_bin_sum <= neighborhood_per_bin_sum && outdegree_per_bin_sum != total_outdegree) {
        std::cerr << "last check: outdegree per bin sum doesn't match total outdegree" << std::endl;
        std::cerr << "total outdgree is " << std::to_string(total_outdegree) << std::endl;
        std::cerr << "total neighborhood size is " << std::to_string(neighborhood_per_bin_sum) << std::endl;
        for(int bin_index = 0; bin_index < this->num_bins; bin_index ++) {
            std::cerr << "bin: " << std::to_string(bin_index) << " has " << std::to_string(target_outdegree_per_bin_map[bin_index]) << std::endl;
        }
        exit(-1);
    }
    return target_outdegree_per_bin_map;
}

int ABM::GetBinIndex(int year_diff) {
    /* 5,10,25 with explicit but ignored 1*/
    /* [1, 5) */
    /* [5, 10) */
    /* [10, 25) */
    /* [25, inf) */
    for (int i = 0; i < this->num_bins - 1; i ++) {
        if (this->bin_boundaries.at(i) <= year_diff && year_diff < this->bin_boundaries.at(i + 1)) {
            return i;
        }
    }
    return this->bin_boundaries.size();
}

int ABM::GetBinIndex(Graph* graph, int current_node, int current_year) {
    int current_diff = current_year - graph->GetIntAttribute("year", current_node);
    return this->GetBinIndex(current_diff);
}

std::unordered_map<int, std::vector<int>> ABM::BinNeighborhood(Graph* graph, int current_year, std::vector<int> n_hop_list) {
    std::unordered_map<int, std::vector<int>> binned_neighborhood;
    for(int i = 0; i < this->num_bins; i ++) {
        binned_neighborhood[i] = {};
    }
    for(size_t i = 0; i < n_hop_list.size(); i ++) {
        int current_node = n_hop_list[i];
        int current_bin = GetBinIndex(graph, current_node, current_year);
        binned_neighborhood[current_bin].push_back(current_node);
    }
    // MARK: debug
    for(int i = 0; i < this->num_bins; i ++) {
        for(size_t j = 0; j < binned_neighborhood[i].size(); j ++) {
            if (binned_neighborhood[i][j] == 0) {
                std::cerr << "bin index " + std::to_string(i) + " has element 0 at index " + std::to_string(j) << std::endl;
                std::cerr << "bin index " + std::to_string(i) + " has size " + std::to_string(binned_neighborhood[i].size()) << std::endl;
                exit(1);
            }
            /* if (binned_neighborhood[i][j] >= 14694932) { */
            /*     std::cerr << "node id greater than or equal to 14694932 is included in a binned neighborhood" << std::endl; */
            /*     std::cerr << "bin index " + std::to_string(i) + " has size " + std::to_string(binned_neighborhood[i].size()) << std::endl; */
            /*     exit(1); */
            /* } */
        }
    }
    return binned_neighborhood;
}
std::unordered_map<int, std::vector<int>> ABM::GetNeighborhoodMap(Graph* graph, int current_year, const std::vector<int>& generator_nodes, int num_hops) {
    if (this->use_alpha) {
        // create distance 1 and distance 2 neighborhoods
        return this->GetOneAndTwoDistanceNeighborhoods(graph, current_year, generator_nodes, num_hops);
    } else {
        return this->GetNHopNeighborhood(graph, current_year, generator_nodes, num_hops);
    }
}

std::unordered_map<int, std::vector<int>> ABM::GetOneAndTwoDistanceNeighborhoods(Graph* graph, int current_year, const std::vector<int>& generator_nodes, int num_hops) {
    std::unordered_map<int, std::vector<int>> n_hop_map;
    if (this->neighborhood_sample == -1) {
        std::set<int> visited;
        for(size_t i = 0; i < generator_nodes.size(); i ++) {
            int generator_node = generator_nodes.at(i);
            std::queue<std::pair<int, int>> to_visit;
            to_visit.push({generator_node, 0});
            visited.insert(generator_node);
            while(!to_visit.empty()) {
                std::pair<int, int> current_pair = to_visit.front();
                to_visit.pop();
                int current_node = current_pair.first;
                int current_distance = current_pair.second;
                if (current_distance > 0) {
                    n_hop_map[current_distance].push_back(current_node);
                }
                if (current_distance < num_hops) {
                    if(graph->GetOutDegree(current_node) > 0) {
                        for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(current_node)) {
                            if(!visited.contains(outgoing_neighbor)) {
                                visited.insert(outgoing_neighbor);
                                to_visit.push({outgoing_neighbor, current_distance + 1});
                            }
                        }
                    }
                    if(graph->GetInDegree(current_node) > 0) {
                        for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(current_node)) {
                            if(!visited.contains(incoming_neighbor)) {
                                visited.insert(incoming_neighbor);
                                to_visit.push({incoming_neighbor, current_distance + 1});
                            }
                        }
                    }
                }
            }
        }
    } else {
        // NOTE: only supports randomly sampling from up to 2-hop
        n_hop_map[1].reserve(this->neighborhood_sample);
        n_hop_map[2].reserve(this->neighborhood_sample);
        size_t max_neighborhood_size = this->neighborhood_sample;
        std::set<int> visited;
        pcg_extras::seed_seq_from<std::random_device> rand_dev;
        pcg32 generator(rand_dev);
        for(size_t i = 0; i < generator_nodes.size(); i ++) {
            // get the 1-hop first
            int generator_node = generator_nodes.at(i);
            visited.insert(generator_node);
            std::vector<int> current_one_hop_neighborhood;
            if(graph->GetOutDegree(generator_node) > 0) {
                for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(generator_node)) {
                    if (!visited.contains(outgoing_neighbor)) {
                        current_one_hop_neighborhood.push_back(outgoing_neighbor);
                        visited.insert(outgoing_neighbor);
                    }
                }
            }
            if(graph->GetInDegree(generator_node) > 0) {
                for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(generator_node)) {
                    if (!visited.contains(incoming_neighbor)) {
                        current_one_hop_neighborhood.push_back(incoming_neighbor);
                        visited.insert(incoming_neighbor);
                    }
                }
            }
            // pick random nodes to get to 2-hop
            if (current_one_hop_neighborhood.size() > max_neighborhood_size) {
                std::vector<int> sampled_one_hop_neighborhood;
                std::sample(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), std::back_inserter(sampled_one_hop_neighborhood), max_neighborhood_size, generator);
                current_one_hop_neighborhood = sampled_one_hop_neighborhood;
            }
            // until here should be fast
            std::copy(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), std::back_inserter(n_hop_map[1]));
            std::shuffle(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), generator);
            for(size_t j = 0; j < current_one_hop_neighborhood.size(); j ++) {
                int current_two_hop_size = 0;
                if (graph->GetOutDegree(current_one_hop_neighborhood[j]) > 0) {
                    current_two_hop_size += graph->GetForwardAdjMap().at(current_one_hop_neighborhood[j]).size();
                }
                if (graph->GetInDegree(current_one_hop_neighborhood[j]) > 0) {
                    current_two_hop_size += graph->GetBackwardAdjMap().at(current_one_hop_neighborhood[j]).size();
                }
                // worst case for one node we might only grab the outgoing edges
                if (n_hop_map[2].size() + current_two_hop_size < max_neighborhood_size) {
                    if (graph->GetOutDegree(current_one_hop_neighborhood[j]) > 0) {
                        for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(current_one_hop_neighborhood[j])) {
                            if (!visited.contains(outgoing_neighbor)) {
                                visited.insert(outgoing_neighbor);
                                n_hop_map[2].push_back(outgoing_neighbor);
                            }
                        }
                    }
                    if (graph->GetInDegree(current_one_hop_neighborhood[j]) > 0) {
                        for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(current_one_hop_neighborhood[j])) {
                            if (!visited.contains(incoming_neighbor)) {
                                visited.insert(incoming_neighbor);
                                n_hop_map[2].push_back(incoming_neighbor);
                            }
                        }
                    }
                } else {
                    std::vector<int> to_be_sampled_neighborhood;
                    if (graph->GetOutDegree(current_one_hop_neighborhood[j]) > 0) {
                        for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(current_one_hop_neighborhood[j])) {
                            if (!visited.contains(outgoing_neighbor)) {
                                to_be_sampled_neighborhood.push_back(outgoing_neighbor);
                            }
                        }
                    }
                    if (graph->GetInDegree(current_one_hop_neighborhood[j]) > 0) {
                        for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(current_one_hop_neighborhood[j])) {
                            if (!visited.contains(incoming_neighbor)) {
                                to_be_sampled_neighborhood.push_back(incoming_neighbor);
                            }
                        }
                    }
                    std::sample(to_be_sampled_neighborhood.begin(), to_be_sampled_neighborhood.end(), std::back_inserter(n_hop_map[2]), max_neighborhood_size - n_hop_map[2].size(), generator);
                    if (n_hop_map[2].size() == max_neighborhood_size) {
                        return n_hop_map;
                    }
                }
            }
        }
    }
    /* for(int k = 1; k < 3; k ++) { */
    /*     for(int l = 1; l < n_hop_map[k].size(); l ++) { */
            /* if (n_hop_map[k][l] > 14694932) { */
            /*     std::cerr << std::to_string(k) + " distance neighborhood contained a node greater than 14694932 at the end" << std::endl; */
            /*     exit(1); */
            /* } */
        /* } */
    /* } */
    return n_hop_map;
}

std::unordered_map<int, std::vector<int>> ABM::GetNHopNeighborhood(Graph* graph, int current_year, const std::vector<int>& generator_nodes, int num_hops) {
    std::unordered_map<int, std::vector<int>> n_hop_map;
    std::vector<int> n_hop_neighborhood;
    if (this->neighborhood_sample == -1) {
        std::set<int> visited;
        for(size_t i = 0; i < generator_nodes.size(); i ++) {
            int generator_node = generator_nodes.at(i);
            std::queue<std::pair<int, int>> to_visit;
            to_visit.push({generator_node, 0});
            visited.insert(generator_node);
            while(!to_visit.empty()) {
                std::pair<int, int> current_pair = to_visit.front();
                to_visit.pop();
                int current_node = current_pair.first;
                int current_distance = current_pair.second;
                if (current_distance > 0) {
                    n_hop_neighborhood.push_back(current_node);
                }
                if (current_distance < num_hops) {
                    if(graph->GetOutDegree(current_node) > 0) {
                        for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(current_node)) {
                            if(!visited.contains(outgoing_neighbor)) {
                                visited.insert(outgoing_neighbor);
                                to_visit.push({outgoing_neighbor, current_distance + 1});
                            }
                        }
                    }
                    if(graph->GetInDegree(current_node) > 0) {
                        for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(current_node)) {
                            if(!visited.contains(incoming_neighbor)) {
                                visited.insert(incoming_neighbor);
                                to_visit.push({incoming_neighbor, current_distance + 1});
                            }
                        }
                    }
                }
            }
        }
    } else {
        // NOTE: only supports randomly sampling from up to 2-hop
        n_hop_neighborhood.reserve(this->neighborhood_sample);
        size_t max_neighborhood_size = this->neighborhood_sample;
        std::set<int> visited;
        pcg_extras::seed_seq_from<std::random_device> rand_dev;
        pcg32 generator(rand_dev);
        for(size_t i = 0; i < generator_nodes.size(); i ++) {
            // get the 1-hop first
            int generator_node = generator_nodes.at(i);
            visited.insert(generator_node);
            std::vector<int> current_one_hop_neighborhood;
            if(graph->GetOutDegree(generator_node) > 0) {
                for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(generator_node)) {
                    if (!visited.contains(outgoing_neighbor)) {
                        current_one_hop_neighborhood.push_back(outgoing_neighbor);
                        visited.insert(outgoing_neighbor);
                    }
                }
            }
            if(graph->GetInDegree(generator_node) > 0) {
                for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(generator_node)) {
                    if (!visited.contains(incoming_neighbor)) {
                        current_one_hop_neighborhood.push_back(incoming_neighbor);
                        visited.insert(incoming_neighbor);
                    }
                }
            }
            // pick random nodes to get to 2-hop
            if (current_one_hop_neighborhood.size() > max_neighborhood_size) {
                std::vector<int> sampled_one_hop_neighborhood;
                std::sample(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), std::back_inserter(sampled_one_hop_neighborhood), max_neighborhood_size, generator);
                current_one_hop_neighborhood = sampled_one_hop_neighborhood;
            }
            // until here should be fast
            std::copy(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), std::back_inserter(n_hop_neighborhood));
            std::shuffle(current_one_hop_neighborhood.begin(), current_one_hop_neighborhood.end(), generator);
            for(size_t j = 0; j < current_one_hop_neighborhood.size(); j ++) {
                int current_two_hop_size = 0;
                if (graph->GetOutDegree(current_one_hop_neighborhood[j]) > 0) {
                    current_two_hop_size += graph->GetForwardAdjMap().at(current_one_hop_neighborhood[j]).size();
                }
                if (graph->GetInDegree(current_one_hop_neighborhood[j]) > 0) {
                    current_two_hop_size += graph->GetBackwardAdjMap().at(current_one_hop_neighborhood[j]).size();
                }
                // worst case for one node we might only grab the outgoing edges
                if (n_hop_neighborhood.size() + current_two_hop_size < max_neighborhood_size) {
                    if (graph->GetOutDegree(current_one_hop_neighborhood[j]) > 0) {
                        for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(current_one_hop_neighborhood[j])) {
                            if (!visited.contains(outgoing_neighbor)) {
                                visited.insert(outgoing_neighbor);
                                n_hop_neighborhood.push_back(outgoing_neighbor);
                            }
                        }
                    }
                    if (graph->GetInDegree(current_one_hop_neighborhood[j]) > 0) {
                        for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(current_one_hop_neighborhood[j])) {
                            if (!visited.contains(incoming_neighbor)) {
                                visited.insert(incoming_neighbor);
                                n_hop_neighborhood.push_back(incoming_neighbor);
                            }
                        }
                    }
                } else {
                    std::vector<int> to_be_sampled_neighborhood;
                    if (graph->GetOutDegree(current_one_hop_neighborhood[j]) > 0) {
                        for(auto const& outgoing_neighbor : graph->GetForwardAdjMap().at(current_one_hop_neighborhood[j])) {
                            if (!visited.contains(outgoing_neighbor)) {
                                to_be_sampled_neighborhood.push_back(outgoing_neighbor);
                            }
                        }
                    }
                    if (graph->GetInDegree(current_one_hop_neighborhood[j]) > 0) {
                        for(auto const& incoming_neighbor : graph->GetBackwardAdjMap().at(current_one_hop_neighborhood[j])) {
                            if (!visited.contains(incoming_neighbor)) {
                                to_be_sampled_neighborhood.push_back(incoming_neighbor);
                            }
                        }
                    }
                    std::sample(to_be_sampled_neighborhood.begin(), to_be_sampled_neighborhood.end(), std::back_inserter(n_hop_neighborhood), max_neighborhood_size - n_hop_neighborhood.size(), generator);
                    if (n_hop_neighborhood.size() == max_neighborhood_size) {
                        n_hop_map[1] = n_hop_neighborhood;
                        return n_hop_map;
                    }
                }
            }
        }
    }
    n_hop_map[1] = n_hop_neighborhood;
    return n_hop_map;
}

void ABM::LogTime(int current_year, std::string label, int time_elapsed) {
    this->timing_file_handle << std::to_string(current_year) << "," << label << "," << std::to_string(time_elapsed) << "\n";
    std::flush(this->timing_file_handle);
}

void ABM::LogTime(int current_year, std::string label) {
    std::chrono::time_point<std::chrono::steady_clock> current_time = std::chrono::steady_clock::now();
    auto milliseconds_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - this->prev_time);
    this->timing_map[current_year].push_back({label, milliseconds_elapsed.count()});
    this->prev_time = current_time;
    this->timing_file_handle << std::to_string(current_year) << "," << label << "," << std::to_string(milliseconds_elapsed.count()) << "\n";
    std::flush(this->timing_file_handle);
}


std::chrono::time_point<std::chrono::steady_clock> ABM::LocalLogTime(std::vector<std::pair<std::string, int>>& local_parallel_stage_time_vec, std::chrono::time_point<std::chrono::steady_clock> local_prev_time, std::string label) {
    std::chrono::time_point<std::chrono::steady_clock> local_current_time = std::chrono::steady_clock::now();
    auto milliseconds_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(local_current_time - local_prev_time);
    local_parallel_stage_time_vec.push_back({label, milliseconds_elapsed.count()});
    return local_current_time;
}

void ABM::WriteTimingFile(int start_year, int end_year) {
    for(int i = start_year; i < end_year; i ++) {
        for(size_t j = 0; j < this->timing_map[i].size(); j ++) {
            this->timing_file_handle << std::to_string(i) << "," << (this->timing_map[i][j]).first << "," << std::to_string((this->timing_map[i][j]).second) << "\n";
        }
    }
}

std::unordered_map<int, int> ABM::GetNumCitationsPerNeighborhood(double alpha, int total_num_citations_neighborhood, const std::unordered_map<int, std::vector<int>>& n_hop_map) {
    std::unordered_map<int, int> num_citations_per_neighborhood;
    if (this->use_alpha) {
        num_citations_per_neighborhood[1] = std::min((int)(total_num_citations_neighborhood * alpha), (int)n_hop_map.at(1).size());
        num_citations_per_neighborhood[2] = std::min(total_num_citations_neighborhood - num_citations_per_neighborhood[1], (int)n_hop_map.at(2).size());
    } else {
        num_citations_per_neighborhood[1] = std::min(total_num_citations_neighborhood, (int)n_hop_map.at(1).size());
    }
    return num_citations_per_neighborhood;
}


void ABM::InitializeBinBoundaries() {
    std::string bin_boundary_string = this->recency_bins;
    std::stringstream ss(bin_boundary_string);
    std::string current_value;
    int element_no = 0;
    while(std::getline(ss, current_value, ',')) {
        /* if (element_no > 0) { */
        int current_bin_value = std::stoi(current_value);
        if (current_bin_value != -1) {
            this->bin_boundaries.push_back(current_bin_value);
        }
        /* } */
        element_no ++;
    }
    this->num_bins = element_no;
}

bool ABM::ValidateBinBoundaries() {
    this->WriteToLogFile(std::to_string(this->bin_boundaries.size()) + " bins have been created", Log::info);
    if (this->bin_boundaries.size() == 0) {
        this->WriteToLogFile("At least one bin is required", Log::error);
        return false;
    }
    if (this->bin_boundaries.at(0) != 1) {
        this->WriteToLogFile("The first bin must start with year 1", Log::error);
        return false;
    }
    std::string recency_bin_string;
    for(size_t i = 0; i < this->bin_boundaries.size() - 1; i ++) {
        recency_bin_string += ("[" + std::to_string(bin_boundaries.at(i)) + "," + std::to_string(this->bin_boundaries.at(i + 1)) + "), ");
    }
    recency_bin_string += ("[" + std::to_string(this->bin_boundaries.at(this->bin_boundaries.size() - 1)) + ",infinity)");
    this->WriteToLogFile("Here are the bins: " + recency_bin_string, Log::info);
    return true;
}

bool ABM::ValidateArguments() {
    if (!this->ValidateArgument("Environment", "edgelist", this->edgelist, "NOTFOUND")) {
        return false;
    }
    if (!this->ValidateArgument("Environment", "nodelist", this->nodelist, "NOTFOUND")) {
        return false;
    }
    if (!this->ValidateArgument("Environment", "growth_rate", this->growth_rate, -42)) {
        return false;
    }
    if (!this->ValidateArgument("Environment", "num_cycles", this->num_cycles, -42)) {
        return false;
    }
    if (!this->ValidateArgument("Environment", "out_degree_bag", this->out_degree_bag, "NOTFOUND")) {
        return false;
    }
    if (!this->ValidateArgument("Environment", "recency_table", this->recency_table, "NOTFOUND")) {
        return false;
    }
    if (this->planted_nodes == "NOTFOUND") {
        this->WriteToLogFile("No agents will be planted", Log::info);
    } else {
        this->WriteToLogFile("planted_nodes: " + this->planted_nodes, Log::info);
    }
    if (!this->ValidateArgument("Agent", "fully_random_citations", this->fully_random_citations, -42)) {
        return false;
    }
    if (this->preferential_weight == -42) {
        this->WriteToLogFile("Required parameter 'preferential_weight' was not found in the 'Agent' section", Log::error);
        return false;
    } else if (this->preferential_weight == -1) {
        this->WriteToLogFile("preferential_weight: randomized", Log::info);
    } else {
        this->WriteToLogFile("preferential_weight: " + std::to_string(this->preferential_weight), Log::info);
    }
    if (this->fitness_weight == -42) {
        this->WriteToLogFile("Required parameter 'fitness_weight' was not found in the 'Agent' section", Log::error);
        return false;
    } else if (this->fitness_weight == -1) {
        this->WriteToLogFile("fitness_weight: randomized", Log::info);
    } else {
        this->WriteToLogFile("fitness_weight: " + std::to_string(this->fitness_weight), Log::info);
    }
    if (!this->ValidateArgument("Agent", "fitness_value_min", this->fitness_value_min, -42)) {
        return false;
    }
    if (!this->ValidateArgument("Agent", "fitness_value_max", this->fitness_value_max, -42)) {
        return false;
    }
    if (!this->ValidateArgument("Agent", "same_year_citations", this->same_year_citations, -42)) {
        return false;
    }
    if (!this->ValidateArgument("Agent", "neighborhood_sample", this->neighborhood_sample, -42)) {
        return false;
    }
    if (this->use_alpha) {
        if (this->alpha == -42) {
            this->WriteToLogFile("Required parameter 'alpha' was not found in the 'Agent' section while 'use_alpha' was true", Log::error);
            return false;
        } else if (this->alpha == -1) {
            this->WriteToLogFile("alpha: randomized", Log::info);
        } else {
            this->WriteToLogFile("alpha: " + std::to_string(this->alpha), Log::info);
        }
    } else {
        if (this->alpha == -42) {
            this->WriteToLogFile("Alpha ignored. Agents will not split the neighborhood into 1 and 2 distance neighborhoods", Log::info);
        } else {
            this->WriteToLogFile("'use_alpha' is false but a value for 'alpha' was provided", Log::error);
            return false;
        }
    }
    if (this->start_from_checkpoint) {
        this->WriteToLogFile("Starting from a checkpoint. Make sure the edgelist and nodelist provided are the result of a previous simulation", Log::info);
    } else {
        this->WriteToLogFile("Not using checkpointing. Starting simulation from the first year.", Log::info);
    }
    if (!this->ValidateArgument("General", "output_file", this->output_file, "NOTFOUND")) {
        return false;
    }
    if (!this->ValidateArgument("General", "recency_bins", this->recency_bins, "NOTFOUND")) {
        return false;
    }
    if (!this->ValidateArgument("General", "auxiliary_information_file", this->auxiliary_information_file, "NOTFOUND")) {
        return false;
    }
    if (!this->ValidateArgument("General", "log_file", this->log_file, "NOTFOUND")) {
        return false;
    }
    if (!this->ValidateArgument("General", "num_processors", this->num_processors, -42)) {
        return false;
    }
    if (!this->ValidateArgument("General", "log_level", this->log_level, -42)) {
        return false;
    }
    return true;
}


int ABM::main() {
    /* std::cerr << "running with asserts" << std::endl; */
    if (!this->ValidateBinBoundaries()) {
        return 1;
    }
    return 1;
    Graph* graph = new Graph(this->edgelist, this->nodelist, this->start_from_checkpoint);
    this->WriteToLogFile("loaded graph", Log::info);
    /* if (!this->start_from_checkpoint) { */
    /*     this->InitializeSeedFitness(graph); */
    /* } */
    /* this->InitializeFitness(graph); */

    /* node ids to continous integer from 0 */
    std::unordered_map<int, int> continuous_node_mapping = this->BuildContinuousNodeMapping(graph);

    /* continous integer from 0 to node ids*/
    std::unordered_map<int, int> reverse_continuous_node_mapping = this->ReverseMapping(continuous_node_mapping);

    int start_year = this->GetMaxYear(graph) + 1;
    int next_node_id = this->GetMaxNode(graph) + 1;
    int initial_next_node_id = next_node_id;

    /* get input to score arrays based on continuous_node_mapping */
    int initial_graph_size = graph->GetNodeSet().size();
    int final_graph_size = this->GetFinalGraphSize(graph);
    this->WriteToLogFile("final graph size is " + std::to_string(final_graph_size), Log::info);
    int* in_degree_arr = new int[final_graph_size];
    int* fitness_arr = new int[final_graph_size];
    double* pa_arr = new double[final_graph_size];
    double* fit_arr = new double[final_graph_size];
    double* random_weight_arr = new double[final_graph_size];
    double* current_score_arr = new double[final_graph_size];

    // the first new agent node has index 0 but is actually index initial_graph_size in the continuous mapping
    double* pa_weight_arr = new double[final_graph_size - initial_graph_size];
    double* fit_weight_arr = new double[final_graph_size - initial_graph_size];
    int* out_degree_arr = new int[final_graph_size - initial_graph_size];
    double* alpha_arr = new double[final_graph_size - initial_graph_size];

    this->PopulateWeightArrs(pa_weight_arr, fit_weight_arr, final_graph_size - initial_graph_size);
    this->PopulateAlphaArr(alpha_arr, final_graph_size - initial_graph_size);
    this->PopulateOutDegreeArr(out_degree_arr, final_graph_size - initial_graph_size);


    std::vector<int> new_nodes_vec;
    std::vector<std::pair<int, int>> new_edges_vec;
    std::set<int> same_year_source_nodes;
    std::unordered_map<int, double> tanh_cached_results;
    for(int i = 0; i < 1000; i ++) {
        tanh_cached_results[i] = this->peak_constant * std::tanh((pow(i, 3)/this->delay_constant)*(1/this->peak_constant));
        /* cached_results[i] = (i * this->gamma) + 1; */
    }
    std::unordered_map<int, double> exp_cached_results;
    for(int i = 0; i < 1000; i ++) {
        exp_cached_results[i] = std::max(pow(i, this->gamma), 1.0) + 1;
        /* cached_results[i] = (i * this->gamma) + 1; */
    }
    Eigen::setNbThreads(1);
    for (int current_year = start_year; current_year < start_year + this->num_cycles; current_year ++) {
        this->prev_time = std::chrono::steady_clock::now();
        int current_graph_size = graph->GetNodeSet().size();
        this->WriteToLogFile("current year is: " + std::to_string(current_year) + " and the graph is " + std::to_string(current_graph_size) + " nodes large", Log::info);
        this->FillInDegreeArr(graph, continuous_node_mapping, in_degree_arr);
        this->LogTime(current_year, "Fill in-degree array");
        this->FillFitnessArr(graph, continuous_node_mapping, current_year, fitness_arr);
        this->LogTime(current_year, "Fill fitness array");
        std::unordered_map<int, double> binned_recency_probabilities = GetBinnedRecencyProbabilities(graph, current_year); // MARK: probbaly don't need to re-generate
        this->CalculateExpScores(exp_cached_results, in_degree_arr, pa_arr, current_graph_size);
        this->LogTime(current_year, "Process in-degree array");
        /* this->CalculateTanhScores(tanh_cached_results, fitness_arr, fit_arr, current_graph_size); */
        this->CalculateExpScores(exp_cached_results, fitness_arr, fit_arr, current_graph_size);
        this->LogTime(current_year, "Process fitness array");

        /* initialize new nodes */
        int num_new_nodes = std::ceil(current_graph_size * this->growth_rate);
        this->WriteToLogFile("making " + std::to_string(num_new_nodes) + " nodes this year", Log::info);
        for(int i = 0; i < num_new_nodes; i ++) {
            continuous_node_mapping[next_node_id] = current_graph_size + i;
            reverse_continuous_node_mapping[current_graph_size + i] = next_node_id;
            new_nodes_vec.push_back(next_node_id);
            graph->SetIntAttribute("year", next_node_id, current_year);
            graph->SetStringAttribute("type", next_node_id, "agent");
            next_node_id ++;
        }
        this->LogTime(current_year, "Create new node ids");
        this->FillSameYearSourceNodes(same_year_source_nodes, new_nodes_vec.size());
        this->LogTime(current_year, "Pick same year nodes");

        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            int new_node = new_nodes_vec[i];
            std::vector<int> generator_nodes = this->GetGeneratorNodes(graph, reverse_continuous_node_mapping);
            this->UpdateGraphAttributesGeneratorNodes(graph, new_node, generator_nodes);
        }
        this->LogTime(current_year, "Pick generator nodes");
        std::vector<int> sampled_neighborhood_sizes_map(new_nodes_vec.size());
        std::vector<int> fully_random_citations_map(new_nodes_vec.size());
        std::vector<std::vector<int>> sampled_binned_neighborhood_sizes_map(new_nodes_vec.size());

        std::vector<std::pair<std::string, int>> parallel_stage_time_vec;
        /* #pragma omp parallel for reduction(custom_merge_vec_int: new_edges_vec) reduction(merge_str_int_pair_vecs: parallel_stage_time_vec) */
        /* #pragma omp parallel for reduction(merge_int_pair_vecs: new_edges_vec) reduction(merge_str_int_pair_vecs: parallel_stage_time_vec) */
        #pragma omp parallel for reduction(custom_merge_vec_int: new_edges_vec) reduction(merge_str_int_pair_vecs: parallel_stage_time_vec)
        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            /* std::cerr << "starting new node" << std::endl; */
            std::chrono::time_point<std::chrono::steady_clock> local_prev_time = std::chrono::steady_clock::now();
            std::vector<std::pair<int, int>> local_new_edges_vec;
            std::vector<std::pair<std::string, int>> local_parallel_stage_time_vec;


            int citations[250]; // out-degree assumed to be max 249 MARK: macro 250 or better would be parse outdegree bag and set to max outdegree
            int new_node = new_nodes_vec[i];
            // continuous_node_mapping = node id -> 0..n but guaranteed 0 .. initial graph size are seed nodes
            // initial graphsize .. n are agent nodes
            int weight_arr_index = continuous_node_mapping[new_node] - initial_graph_size;
            double pa_weight = pa_weight_arr[weight_arr_index];
            double fit_weight = fit_weight_arr[weight_arr_index];
            double alpha = alpha_arr[weight_arr_index];
            std::vector<int> generator_nodes = this->GetGraphAttributesGeneratorNodes(graph, new_node);
            int num_hops = 2;
            // if use alpha then map has keys 1 and 2
            // if use alpha false then map has only key 1
            std::unordered_map<int, std::vector<int>> n_hop_map = this->GetNeighborhoodMap(graph, current_year, generator_nodes, num_hops);
            // out-degree assigned is D
            // G=1 for generator node
            // S=0 or 1 sometimes if i'm a same year source
            // 5% or 0.05 is default proprotion of random so R = (D * 0.05)
            // citing from neighborhood based on scores is N = D - G - S - R
            // if use alpha true then there's 2 neighborhoods so N * alpha for 1-hop N * (1- alpha) for 2-hop
            // if use alpha false then there's 1 neighborhood so cite N things from there

            int num_generator_node_citation = generator_nodes.size(); // should be 1 for now
            int same_year_citation = same_year_source_nodes.count(i); // could be 0 or 1
            int num_fully_random_cited_reserved = std::floor(this->fully_random_citations * out_degree_arr[weight_arr_index]); // e.g., 5% of out-degree. some small number
            int total_num_citations_neighborhood = out_degree_arr[weight_arr_index] - num_generator_node_citation - same_year_citation - num_fully_random_cited_reserved;
            std::unordered_map<int, int> num_citations_per_neighborhood = this->GetNumCitationsPerNeighborhood(alpha, total_num_citations_neighborhood, n_hop_map);

            int num_actually_cited = 0;
            if (same_year_citation) {
                num_actually_cited += this->MakeSameYearCitations(same_year_source_nodes, new_nodes_vec.size(), reverse_continuous_node_mapping, citations, current_graph_size);
            }
            /* std::cerr << "num actually cited at " + std::to_string(num_actually_cited) + " after same year citations" << std::endl; */
            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make same year citations");

            for(size_t current_neighborhood_index = 1; current_neighborhood_index < n_hop_map.size() + 1; current_neighborhood_index ++) { // 2 iter if use alpha true
                std::unordered_map<int, std::vector<int>> binned_neighborhood = this->BinNeighborhood(graph, current_year, n_hop_map.at(current_neighborhood_index));
                local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "bin neighborhood");

                std::unordered_map<int, int> outdegree_per_bin_map = this->BinOutdegrees(binned_neighborhood, num_citations_per_neighborhood.at(current_neighborhood_index), binned_recency_probabilities);
                for(int bin_index = 0; bin_index < this->num_bins - 1; bin_index ++) { // if there's only 1 bin then this is always false
                    num_actually_cited += this->MakeCitations(graph, continuous_node_mapping, current_year, binned_neighborhood[bin_index], citations + num_actually_cited, pa_arr, fit_arr, pa_weight, fit_weight, current_graph_size, outdegree_per_bin_map[bin_index]);
                    /* std::cerr << "num actually cited at " + std::to_string(num_actually_cited) + " after citing from bin index" + std::to_string(bin_index) << std::endl; */
                }
                num_actually_cited += this->MakeUniformRandomCitations(graph, continuous_node_mapping, current_year, binned_neighborhood[this->num_bins - 1], citations + num_actually_cited, current_graph_size, outdegree_per_bin_map[this->num_bins - 1]);
                /* std::cerr << "num actually cited at " + std::to_string(num_actually_cited) + " after citing from bin index" + std::to_string(this->num_bins - 1) << std::endl; */
            }

            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make neighborhood citations");

            int num_fully_random_cited = out_degree_arr[weight_arr_index] - num_actually_cited - generator_nodes.size(); // finalized later
            fully_random_citations_map[i] = num_fully_random_cited;
            num_actually_cited += this->MakeUniformRandomCitationsFromGraph(graph, reverse_continuous_node_mapping, generator_nodes, citations, num_actually_cited, num_fully_random_cited);
            /* std::cerr << "num actually cited at " + std::to_string(num_actually_cited) + " after citing from graph" << std::endl; */
            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "make random citations");

            for(size_t j = 0; j < generator_nodes.size(); j ++) {
                local_new_edges_vec.push_back({new_node, generator_nodes[j]});
            }
            for(int j = 0; j < num_actually_cited; j ++) {
                local_new_edges_vec.push_back({new_node, citations[j]});
                /* if ((j > 0 && citations[j] >= initial_next_node_id) || citations[j] == 0) { */
                if (citations[j] == 0) {
                    std::cerr << "citing node id " + std::to_string(citations[j]) << std::endl;
                    std::cerr << "initial next node id was " + std::to_string(initial_next_node_id) << std::endl;
                    std::cerr << "it was at index " + std::to_string(j) << std::endl;
                    std::cerr << "there are " + std::to_string(this->num_bins) + " bins total" << std::endl;
                    exit(1);
                }
            }
            new_edges_vec.insert(new_edges_vec.end(), local_new_edges_vec.begin(), local_new_edges_vec.end());
            local_prev_time = this->LocalLogTime(local_parallel_stage_time_vec, local_prev_time, "record edges");
            parallel_stage_time_vec.insert(parallel_stage_time_vec.end(), local_parallel_stage_time_vec.begin(), local_parallel_stage_time_vec.end());
            /* std::cerr << "finished node" << std::endl; */
        } // end of omp parallel loop
        std::map<std::string, int> per_stage_time_map;
        for(size_t i = 0; i < parallel_stage_time_vec.size(); i ++) {
            per_stage_time_map[parallel_stage_time_vec[i].first] +=  parallel_stage_time_vec[i].second;
        }
        this->LogTime(current_year, "find neighborhood", per_stage_time_map["retrieve neighborhood"]);
        this->LogTime(current_year, "bin neighborhood", per_stage_time_map["bin neighborhood"]);
        this->LogTime(current_year, "make same year citations", per_stage_time_map["make same year citations"]);
        this->LogTime(current_year, "make neighborhood citations", per_stage_time_map["make neighborhood citations"]);
        this->LogTime(current_year, "make random citations", per_stage_time_map["make random citations"]);
        this->LogTime(current_year, "record edges", per_stage_time_map["record edges"]);
        for(size_t i = 0; i < new_edges_vec.size(); i ++) {
            int new_node = new_edges_vec[i].first;
            int destination_id = new_edges_vec[i].second;
            graph->AddEdge({new_node, destination_id});
        }
        this->LogTime(current_year, "Add edges to graph");

        for(size_t i = 0; i < new_nodes_vec.size(); i ++) {
            int new_node = new_nodes_vec[i];
            graph->SetIntAttribute("sampled_neighborhood_size", new_node, sampled_neighborhood_sizes_map[i]);
            graph->SetIntAttribute("fully_random_citations", new_node, fully_random_citations_map[i]);
        }

        this->LogTime(current_year, "Update graph attributes (neighborhood sizes)");
        this->AssignPeakFitnessValues(graph, new_nodes_vec);
        this->AssignFitnessLagDuration(graph, new_nodes_vec);
        this->AssignFitnessPeakDuration(graph, new_nodes_vec);
        this->PlantNodes(graph, new_nodes_vec, current_year - start_year + 1);
        this->LogTime(current_year, "Assign fitness values to new nodes");
        new_nodes_vec.clear();
        new_edges_vec.clear();
        same_year_source_nodes.clear();
    }

    this->WriteToLogFile("finished sim", Log::info);
    graph->WriteGraph(this->output_file);

    this->UpdateGraphAttributesWeights(graph, initial_next_node_id, pa_weight_arr, fit_weight_arr, final_graph_size - initial_graph_size);
    this->UpdateGraphAttributesOutDegrees(graph, initial_next_node_id, out_degree_arr, final_graph_size - initial_graph_size);

    for(auto const& node_id : graph->GetNodeSet()) {
        graph->SetIntAttribute("in_degree", node_id, graph->GetInDegree(node_id));
        graph->SetIntAttribute("out_degree", node_id, graph->GetOutDegree(node_id));
    }

    graph->WriteAttributes(this->auxiliary_information_file);
    this->WriteToLogFile("wrote graph", Log::info);
    delete[] in_degree_arr;
    delete[] fitness_arr;
    delete[] pa_arr;
    delete[] fit_arr;
    delete[] pa_weight_arr;
    delete[] fit_weight_arr;
    delete[] out_degree_arr;
    delete[] random_weight_arr;
    delete[] current_score_arr;
    delete[] alpha_arr;

    delete graph;
    return 0;
}
