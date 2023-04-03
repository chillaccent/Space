#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <vector>
#include <cmath>

using namespace std;

struct ProblemSpace {
    double L; // Length of the cube
    vector<vector<double>> puzzle_piece_positions; // Positions of the puzzle pieces
    vector<vector<double>> cubesat_positions; // Positions of the CubeSats
};

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_weight_t, double>> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;

double distance(const vector<double>& p1, const vector<double>& p2) {
    double sum = 0;
    for (int i = 0; i < p1.size(); ++i) {
        sum += pow(p1[i] - p2[i], 2);
    }
    return sqrt(sum);
}
vector<Vertex> prm(const ProblemSpace& problem_space, double distance_threshold) {
    Graph graph;
    vector<Vertex> puzzle_piece_vertices;
    vector<Vertex> cubesat_vertices;

    // Add puzzle piece vertices to the graph
    for (auto pos : problem_space.puzzle_piece_positions) {
        Vertex v = boost::add_vertex(graph);
        boost::put(boost::vertex_index, graph, v, puzzle_piece_vertices.size());
        puzzle_piece_vertices.push_back(v);
    }

    // Add CubeSat vertices to the graph
    for (auto pos : problem_space.cubesat_positions) {
        Vertex v = boost::add_vertex(graph);
        boost::put(boost::vertex_index, graph, v, puzzle_piece_vertices.size() + cubesat_vertices.size());
        cubesat_vertices.push_back(v);
    }
    // Add edges between puzzle piece and CubeSat vertices that are within a certain distance threshold
    for (int i = 0; i < puzzle_piece_vertices.size(); ++i) {
        for (int j = 0; j < cubesat_vertices.size(); ++j) {
            double dist = distance(problem_space.puzzle_piece_positions[i], problem_space.cubesat_positions[j]);
            if (dist <= distance_threshold) {
                Edge e;
                bool success;
                tie(e, success) = boost::add_edge(puzzle_piece_vertices[i], cubesat_vertices[j], graph);
                boost::put(boost::edge_weight, graph, e, dist);
            }
        }
    }
    // Add edges between puzzle piece vertices that are within a certain distance threshold
    for (int i = 0; i < puzzle_piece_vertices.size(); ++i) {
        for (int j = i + 1; j < puzzle_piece_vertices.size(); ++j) {
            double dist = distance(problem_space.puzzle_piece_positions[i], problem_space.puzzle_piece_positions[j]);
            if (dist <= distance_threshold) {
                Edge e;
                bool success;
                tie(e, success) = boost::add_edge(puzzle_piece_vertices[i], puzzle_piece_vertices[j], graph);
                boost::put(boost::edge_weight, graph, e, dist);
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
    vector<bool> visited(boost::num_vertices(graph), false);
    priority_queue<pair<double, Vertex>, vector<pair<double, Vertex>>, greater<pair<double, Vertex>>> frontier;
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
    vector<double> distances(boost::num_vertices(graph));
    vector<bool> visited(boost::num_vertices(graph), false);
    priority_queue<pair<double, Vertex>, vector<pair<double, Vertex>>, greater<pair<double, Vertex>>> frontier;
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
        visited[current] = true;
        for (auto neighbor : boost::adjacent_vertices(current, graph)) {
            double new_distance = distances[current] + boost::get(boost::edge_weight, graph, boost::edge(current, neighbor, graph).first);
            if (!visited[neighbor] || new_distance < distances[neighbor]) {
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


