#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <algorithm>
#include <numeric>
#include <unordered_set>

using namespace std;

struct ProblemSpace {
    double L; // Length of the cube
    std::vector<double> puzzle_piece_positions; // Positions of the puzzle pieces
    std::vector<double> cubesat_positions; // Positions of the CubeSats
};


typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_weight_t, double>> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;

double distance(const vector<double>& p1, const vector<double>& p2) {
    if (p1.size() != p2.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    double sum = 0;
    for (int i = 0; i < p1.size(); ++i) {
        sum += pow(p1[i] - p2[i], 2);
    }
    return sqrt(sum);
}




// Define the objective function to be optimized
double objective_function(const vector<double>& x) {
    return -(x[0]*x[0] + x[1]*x[1]);
}

// Define the Cross-Entropy Method algorithm
std::vector<double> cross_entropy_method(int n_iterations, int n_samples, double elite_percent, int n_parameters, double mean_init, double std_init, double learning_rate) {
    // Initialize the mean and standard deviation of the sampling distribution
    std::vector<double> mean(n_parameters, mean_init);
    std::vector<double> std(n_parameters, std_init);

    // Define a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0, 1);

    // Iterate over the specified number of iterations
    for (int i = 0; i < n_iterations; i++) {
        // Generate samples from the sampling distribution
        std::vector<std::vector<double>> samples(n_samples, std::vector<double>(n_parameters));
        for (int j = 0; j < n_samples; j++) {
            for (int k = 0; k < n_parameters; k++) {
                samples[j][k] = mean[k] + std[k] * dist(gen);
            }
        }

        // Evaluate the objective function for each sample
        std::vector<double> values(n_samples);
        for (int j = 0; j < n_samples; j++) {
            values[j] = objective_function(samples[j]);
        }

        // Determine the elite samples
        int n_elite = n_samples * elite_percent;
        std::vector<double> elite_values(n_elite);
        std::vector<std::vector<double>> elite_samples(n_elite, std::vector<double>(n_parameters));
        std::partial_sort_copy(values.begin(), values.end(), elite_values.begin(), elite_values.end(), std::greater<double>());
        for (int j = 0; j < n_elite; j++) {
            auto it = std::find(values.begin(), values.end(), elite_values[j]);
            int index = std::distance(values.begin(), it);
            elite_samples[j] = samples[index];
        }

        // Update the mean and standard deviation of the sampling distribution
        std::vector<double> new_mean(n_parameters);
        for (int j = 0; j < n_parameters; j++) {
            new_mean[j] = std::accumulate(elite_samples.begin(), elite_samples.end(), 0.0, [j](double acc, const std::vector<double>& x){ return acc + x[j]; }) / n_elite;
        }
        std::vector<double> new_std(n_parameters);
        for (int j = 0; j < n_parameters; j++) {
            new_std[j] = std::sqrt(std::accumulate(elite_samples.begin(), elite_samples.end(), 0.0, [j, &new_mean, &elite_samples](double acc, const std::vector<double>& x){
        return acc + std::pow(x[j] - new_mean[j], 2.0);
        }) / (n_elite - 1));
            // Update the mean and standard deviation using the learning rate
    for (int j = 0; j < n_parameters; j++) {
        mean[j] = (1 - learning_rate) * mean[j] + learning_rate * new_mean[j];
        std[j] = (1 - learning_rate) * std[j] + learning_rate * new_std[j];
    }
}

return mean;
    }

vector<Vertex> prm(const ProblemSpace& problem_space, double distance_threshold, const vector<vector<double>>& puzzle_piece_positions) {
    Graph graph;
    vector<Vertex> puzzle_piece_vertices;
    vector<Vertex> cubesat_vertices;

    // Add puzzle piece vertices to the graph
    for (auto pos : puzzle_piece_positions) {
        Vertex v = boost::add_vertex(graph);
        graph[v].pos = pos;
        puzzle_piece_vertices.push_back(v);
    }

    // Add CubeSat vertices to the graph
    for (auto pos : problem_space.cubesat_positions) {
        Vertex v = boost::add_vertex(graph);
        graph[v].pos = pos;
        cubesat_vertices.push_back(v);
    }

    // Connect the vertices in the graph based on distance threshold
    for (int i = 0; i < puzzle_piece_vertices.size(); i++) {
        for (int j = i + 1; j < puzzle_piece_vertices.size(); j++) {
            double d = distance(graph[puzzle_piece_vertices[i]].pos, graph[puzzle_piece_vertices[j]].pos);
            if (d <= distance_threshold) {
                Edge e;
                bool success;
                tie(e, success) = boost::add_edge(puzzle_piece_vertices[i], puzzle_piece_vertices[j], graph);
                if (success) {
                    graph[e].weight = d;
                }
            }
        }
    }

    for (int i = 0; i < cubesat_vertices.size(); i++) {
        for (int j = 0; j < puzzle_piece_vertices.size(); j++) {
            double d = distance(graph[cubesat_vertices[i]].pos, graph[puzzle_piece_vertices[j]].pos);
            if (d <= distance_threshold) {
                Edge e;
                bool success;
                tie(e, success) = boost::add_edge(cubesat_vertices[i], puzzle_piece_vertices[j], graph);
                if (success) {
                    graph[e].weight = d;
                }
            }
        }
    }

    return puzzle_piece_vertices;
}


    // PRM algorithm steps
vector<Vertex> candidate_nodes;
for (int i = 0; i < 1000; ++i) { // Sample 1000 candidate nodes
    Vertex v = boost::add_vertex(graph);
    boost::put(boost::vertex_index, graph, v, puzzle_piece_vertices.size() + cubesat_vertices.size() + candidate_nodes.size());
    candidate_nodes.push_back(v);
    // Connect the candidate node to nearby vertices in the graph
    for (auto pp_v : puzzle_piece_vertices) {
        double dist = distance(problem_space.puzzle_piece_positions[boost::get(boost::vertex_index, graph, pp_v)], problem_space.puzzle_piece_positions[boost::get(boost::vertex_index, graph, v)]);
        if (dist <= distance_threshold) {
            Edge e;
            bool success;
            tie(e, success) = boost::add_edge(v, pp_v, graph);
            boost::put(boost::edge_weight, e, dist);
        }
    }
    for (auto cs_v : cubesat_vertices) {
        double dist = distance(problem_space.cubesat_positions[boost::get(boost::vertex_index, graph, cs_v)], problem_space.puzzle_piece_positions[boost::get(boost::vertex_index, graph, v)]);
        if (dist <= distance_threshold) {
            Edge e;
            bool success;
            tie(e, success) = boost::add_edge(v, cs_v, graph);
            boost::put(boost::edge_weight, e, dist);
        }
    }
}



    // Find a path from the start node to the goal node using IDA*, A*, or Dijkstra
    Vertex start = puzzle_piece_vertices[0];
    Vertex goal = cubesat_vertices[0];
    vector<Vertex> path;
    bool success = false;

    // Try IDA*
    double ida_star_search(Graph& graph, Vertex current, Vertex goal, vector<Vertex>& path, double cost, double bound) {
    double f = cost + distance(current, goal);
    if (f > bound) {
        return f;
    }
    if (current == goal) {
        path.push_back(current);
        return 0;
    }
    double min_cost = numeric_limits<double>::max();
    for (auto neighbor : boost::adjacent_vertices(current, graph)) {
        vector<Vertex>::iterator it = find(path.begin(), path.end(), neighbor);
        if (it == path.end()) {
            path.push_back(neighbor);
            double neighbor_cost = ida_star_search(graph, neighbor, goal, path, cost + boost::get(boost::edge_weight, graph, boost::edge(current, neighbor, graph).first), bound);
            if (neighbor_cost == 0) {
                return 0;
            }
            if (neighbor_cost < min_cost) {
                min_cost = neighbor_cost;
            }
            path.pop_back();
        }
    }
    return min_cost;
}

    bool ida_star_search(Graph& graph, Vertex start, Vertex goal, vector<Vertex>& path) {
    double bound = distance(start, goal);
    while (true) {
        double cost = 0;
        vector<Vertex> temp_path;
        temp_path.push_back(start);
        double result = ida_star_search(graph, start, goal, temp_path, cost, bound);
        if (result == 0) {
            path = temp_path;
            return true;
        }
        if (result == numeric_limits<double>::max()) {
            return false;
        }
        bound = result;
    }
}
    success = ida_star_search(graph, start, goal, path);
    if (!success) {
        // Try A*
    bool a_star_search(Graph& graph, Vertex start, Vertex goal, vector<Vertex>& path) {
        vector<Vertex> predecessors(boost::num_vertices(graph));
        vector<double> distances(boost::num_vertices(graph));
        unordered_map<Vertex, bool> visited;
        boost::heap::fibonacci_heap<pair<double, Vertex>, boost::heap::compare<greater<pair<double, Vertex>>>> frontier;
        distances[start] = 0;
        frontier.push(make_pair(distance(start, goal), start));
        while (!frontier.empty()) {
            Vertex current = frontier.top().second;
            frontier.pop();
            if (current == goal) {
                path.push_back(current);
                for (Vertex v = goal; v != start; v = predecessors[v]) {
                    path.push_back(v);
                }
                path.push_back(start);
                reverse(path.begin(), path.end());
                return true;
            }
            visited[current] = true;
            for (auto neighbor : boost::adjacent_vertices(current, graph)) {
                double new_distance = distances[current] + boost::get(boost::edge_weight, graph, boost::edge(current, neighbor, graph).first);
                if (!visited[neighbor] || new_distance < distances[neighbor]) {
                    distances[neighbor] = new_distance;
                    predecessors[neighbor] = current;
                    frontier.push(make_pair(new_distance + distance(neighbor, goal), neighbor));
                }
            }
        }
        return false;
    }


        success = a_star_search(graph, start, goal, path);
        if (!success) {
            // Try Dijkstra
            bool dijkstra_search(Graph& graph, Vertex start, Vertex goal, vector<Vertex>& path) {
    vector<Vertex> predecessors(boost::num_vertices(graph));
    vector<double> distances(boost::num_vertices(graph), std::numeric_limits<double>::max());
    unordered_map<Vertex, bool> visited;
    boost::heap::fibonacci_heap<pair<double, Vertex>, boost::heap::compare<greater<pair<double, Vertex>>>> frontier;
    distances[start] = 0;
    frontier.push(make_pair(0, start));

    while (!frontier.empty()) {
        Vertex current = frontier.top().second;
        frontier.pop();

        if (current == goal) {
            path.push_back(current);
            for (Vertex v = goal; v != start; v = predecessors[v]) {
                path.push_back(v);
            }
            path.push_back(start);
            reverse(path.begin(), path.end());
            return true;
        }

        if (visited[current]) continue;
        visited[current] = true;

        for (auto neighbor : boost::adjacent_vertices(current, graph)) {
            double new_distance = distances[current] + boost::get(boost::edge_weight, graph, boost::edge(current, neighbor, graph).first);
            if (!visited[neighbor] && new_distance < distances[neighbor]) {
                distances[neighbor] = new_distance;
                predecessors[neighbor] = current;
                frontier.push(make_pair(new_distance, neighbor));
            }
        }
    }
    return false;
}

            success = dijkstra_search(graph, start, goal, path);
        }
    }

    return path;


