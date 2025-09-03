#include "graph.h"


Graph::Graph(std::string edgelist, std::string nodelist, bool start_from_checkpoint): edgelist(edgelist), nodelist(nodelist), start_from_checkpoint(start_from_checkpoint) {
    this->ParseEdgelist();
    this->ParseNodelist();
}

void Graph::ParseEdgelist() {
    char delimiter = Graph::GetDelimiter(this->edgelist);
    std::ifstream input_edgelist(this->edgelist);
    std::string line;
    int line_no = 0;
    while(std::getline(input_edgelist, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        if(line_no != 0) {
            int integer_citing = std::stoi(current_line[0]);
            int integer_cited = std::stoi(current_line[1]);
            this->AddEdge({integer_citing, integer_cited});
        }
        line_no ++;
    }
}

void Graph::ParseNodelist() {
    char delimiter = Graph::GetDelimiter(this->nodelist);
    std::unordered_map<std::string, int> header_to_index_map = Graph::GetHeaderToIndexMap(delimiter, this->nodelist);
    std::ifstream input_nodelist(this->nodelist);
    std::string line;
    int line_no = 0;
    while(std::getline(input_nodelist, line)) {
        std::stringstream ss(line);
        std::string current_value;
        std::vector<std::string> current_line;
        while(std::getline(ss, current_value, delimiter)) {
            current_line.push_back(current_value);
        }
        if (line_no != 0) {
            int integer_node = std::stoi(current_line[header_to_index_map["node_id"]]);
            int integer_year = std::stoi(current_line[header_to_index_map["year"]]);
            this->SetIntAttribute("year", integer_node, integer_year);
            if (this->start_from_checkpoint) {
                std::string type_string = current_line[header_to_index_map["type"]];
                this->SetStringAttribute("type", integer_node, type_string);
                double pa_weight = std::stod(current_line[header_to_index_map["pa_weight"]]);
                this->SetDoubleAttribute("preferential_attachment_weight", integer_node, pa_weight);
                double fit_weight = std::stod(current_line[header_to_index_map["fit_weight"]]);
                this->SetDoubleAttribute("fitness_weight", integer_node, fit_weight);
                int fit_lag_duration = std::stoi(current_line[header_to_index_map["fit_lag_duration"]]);
                this->SetIntAttribute("fitness_lag_duration", integer_node, fit_lag_duration);
                int fit_peak_value = std::stoi(current_line[header_to_index_map["fit_peak_value"]]);
                this->SetIntAttribute("fitness_peak_value", integer_node, fit_peak_value);
                int fit_peak_duration = std::stoi(current_line[header_to_index_map["fit_peak_duration"]]);
                this->SetIntAttribute("fitness_peak_duration", integer_node, fit_peak_duration);
                int assigned_out_degree = std::stoi(current_line[header_to_index_map["assigned_out_degree"]]);
                this->SetIntAttribute("assigned_out_degree", integer_node, assigned_out_degree);
                int planted_nodes_line_number = std::stoi(current_line[header_to_index_map["planted_nodes_line_number"]]);
                this->SetIntAttribute("planted_nodes_line_number", integer_node, planted_nodes_line_number);
                int sampled_neighborhood_size = std::stoi(current_line[header_to_index_map["sampled_neighborhood_size"]]);
                this->SetIntAttribute("sampled_neighborhood_size", integer_node, sampled_neighborhood_size);
                std::string generator_node_string = current_line[header_to_index_map["generator_node_string"]];
                this->SetStringAttribute("generator_node_string", integer_node, generator_node_string);
                int fully_random_citations = std::stoi(current_line[header_to_index_map["fully_random_citations"]]);
                this->SetIntAttribute("fully_random_citations", integer_node, fully_random_citations);
            } else {
                int fitness_lag_uniform = 0; // MARK: hard coded to be static fitness
                int fitness_peak_uniform = 1000; // MARK: hard coded to be static fitness
                int fitness_power = 1;
                this->SetStringAttribute("type", integer_node, "seed");
                this->SetIntAttribute("fitness_lag_duration", integer_node, fitness_lag_uniform);
                this->SetIntAttribute("fitness_peak_duration", integer_node, fitness_peak_uniform);
                this->SetIntAttribute("fitness_peak_value", integer_node, fitness_power);
            }
        }
        line_no ++;
    }
}

void Graph::SetIntAttribute(std::string attribute_key, int node, int attribute_value) {
    this->int_attribute_map[attribute_key][node] = attribute_value;
}

int Graph::GetIntAttribute(std::string attribute_key, int node) const {
    return this->int_attribute_map.at(attribute_key).at(node);
}

void Graph::SetStringAttribute(std::string attribute_key, int node, std::string attribute_value) {
    this->string_attribute_map[attribute_key][node] = attribute_value;
}

std::string Graph::GetStringAttribute(std::string attribute_key, int node) const {
    return this->string_attribute_map.at(attribute_key).at(node);
}

void Graph::SetDoubleAttribute(std::string attribute_key, int node, double attribute_value) {
    this->double_attribute_map[attribute_key][node] = attribute_value;
}

double Graph::GetDoubleAttribute(std::string attribute_key, int node) const {
    return this->double_attribute_map.at(attribute_key).at(node);
}

bool Graph::HasIntAttribute(std::string attribute_key, int node) const {
    return this->int_attribute_map.contains(attribute_key) && this->int_attribute_map.at(attribute_key).contains(node);
}

void Graph::AddEdge(std::pair<int, int> edge) {
    this->forward_adj_map[edge.first].push_back(edge.second);
    this->backward_adj_map[edge.second].push_back(edge.first);
    this->AddNode(edge.first);
    this->AddNode(edge.second);
}

int Graph::GetInDegree(int node) const {
    if (this->backward_adj_map.contains(node)) {
        return this->backward_adj_map.at(node).size();
    }
    return 0;
}

int Graph::GetOutDegree(int node) const {
    if (this->forward_adj_map.contains(node)) {
        return this->forward_adj_map.at(node).size();
    }
    return 0;
}


void Graph::AddNode(int u) {
    this->node_set.insert(u);
}

const std::set<int>& Graph::GetNodeSet() const {
    return this->node_set;
}
const std::unordered_map<int, std::vector<int>>& Graph::GetForwardAdjMap() const {
    return this->forward_adj_map;
}

const std::unordered_map<int, std::vector<int>>& Graph::GetBackwardAdjMap() const {
    return this->backward_adj_map;
}

void Graph::PrintGraph() const {
    for(auto const& [u,u_neighbors] : this->GetForwardAdjMap()) {
        for(const int& v : u_neighbors) {
            /* if (u < v) { */
            std::cout << u << "-" << v << std::endl;
            /* } */
        }
    }
}

void Graph::WriteGraph(std::string output_file) const {
    std::ofstream output_filehandle(output_file);
    output_filehandle << "source,target" << std::endl;
    for(auto const& [u,u_neighbors] : this->GetForwardAdjMap()) {
        for(const int& v : u_neighbors) {
            /* if (u < v) { */
            output_filehandle << u << "," << v << std::endl;
            /* } */
        }
    }
    output_filehandle.close();
}

void Graph::WriteAttributes(std::string auxiliary_information_file) const {
    std::ofstream auxiliary_information_filehandle(auxiliary_information_file);
    auxiliary_information_filehandle << "node_id,type,year,pa_weight,fit_weight,fit_lag_duration,fit_peak_value,fit_peak_duration,in_degree,out_degree,assigned_out_degree,planted_nodes_line_number,generator_node_string,sampled_neighborhood_size,fully_random_citations\n";
    for(const auto& node_id : this->GetNodeSet()) {
        std::string node_type = this->GetStringAttribute("type", node_id);
        int year = this->GetIntAttribute("year", node_id);
        double pa_weight = -1;
        double fit_weight = -1;
        int fit_lag_duration = this->GetIntAttribute("fitness_lag_duration", node_id);
        int fit_peak_value = this->GetIntAttribute("fitness_peak_value", node_id);
        int fit_peak_duration = this->GetIntAttribute("fitness_peak_duration", node_id);
        int out_degree = this->GetIntAttribute("out_degree", node_id);
        int assigned_out_degree  = -1;
        int in_degree = this->GetIntAttribute("in_degree", node_id);
        int planted_nodes_line_number = -1;
        std::string generator_node_string  = "no_generators";
        int neighborhood_size = -1;
        int fully_random_citations = -1;
        if(node_type == "agent") {
            pa_weight = this->GetDoubleAttribute("preferential_attachment_weight", node_id);
            fit_weight = this->GetDoubleAttribute("fitness_weight", node_id);
            assigned_out_degree = this->GetIntAttribute("assigned_out_degree", node_id);
            generator_node_string = this->GetStringAttribute("generator_node_string", node_id);
            if(this->HasIntAttribute("planted_nodes_line_number", node_id)) {
                planted_nodes_line_number = this->GetIntAttribute("planted_nodes_line_number", node_id);
            }
            neighborhood_size = this->GetIntAttribute("sampled_neighborhood_size", node_id);
            fully_random_citations = this->GetIntAttribute("fully_random_citations", node_id);
        }
        auxiliary_information_filehandle << node_id << "," << node_type << "," << year << "," << pa_weight << "," << fit_weight << "," << fit_lag_duration << "," << fit_peak_value << "," << fit_peak_duration << "," << in_degree << "," << out_degree << "," << assigned_out_degree << "," << planted_nodes_line_number << "," << generator_node_string << "," << neighborhood_size << "," << fully_random_citations << "\n";
    }
    auxiliary_information_filehandle.close();
}
