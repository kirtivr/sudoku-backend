/*
 * C++ Program to solve a sudoku and analyze:
 *   1. Whether the solution is unique.
 *   2. How difficult is the sudoku.
 * 
 * Before changes to make priority queue work correctly, time taken was 673 ms.
 * After changes to make priority queue work correctly, time taken was 450 ms.
 * After precomputing some stuff, time taken was 180ms. This was a big win, and quite unexpected to be honest, since the expected
 * asymptotic long pole was elsewhere.
 * Removing -fsanitize=address from the compilation option, time taken was 65ms. The sanitizer has a 3X overhead. Whoosh.
 */

#include <stdio.h>
#include <math.h>
#include <dirent.h>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <optional>
#include <iostream>
#include <fstream>
#include <thread>
#include <filesystem>
#include <queue>
#include <algorithm>
#include <random>
#include <iterator>
#include <chrono>
#include <sstream>
namespace fs = std::filesystem;

typedef struct CommandLineArgs {
    // Folder with all the sudokus which have to be analyzed.
    std::string input_folder;
    std::string output_folder;
    // number of randomized iterations per sudoku.
    uint32_t num_iterations;
    // Extent to which sudoku solving attempts should be randomized.
    uint32_t randomness;
} CommandLineArgs;

std::optional<CommandLineArgs> parse_command_line_args(const std::vector<std::string>& args) {
    CommandLineArgs parsed_args;
    std::string input_folder = args[0];
    if (!std::filesystem::exists(input_folder)) {
        std::cout << "Input folder path " << input_folder << "is invalid\n";
        return std::nullopt;
    }
    std::string output_folder = args[1];
    if (!std::filesystem::exists(output_folder)) {
        std::cout << "Output folder path " << output_folder << "is invalid\n";
        return std::nullopt;
    }
    uint32_t iterations_per_sudoku = 1;
    uint32_t randomness = 0;
    if (args.size() > 2) {
        try {
            int input = std::stoi(args[2]);
            if (input < 0) {
                std::cout << "Iterations per sudoku " << args[2] << " must be positive.\n";
                return std::nullopt;                
            }
            iterations_per_sudoku = (uint32_t)input;
        } catch (const std::exception &ex) {
            std::cout << "Iterations per sudoku " << args[2] << " must be numeric.\n";
            return std::nullopt;
        }
    }
    if (args.size() > 3) {
        try {
            int input = std::stoi(args[3]);
            if (input < 0) {
                std::cout << "Iterations per sudoku " << args[3] << " must be positive.\n";
                return std::nullopt;                
            }
            randomness = (uint32_t)input;
        } catch (const std::exception &ex) {
            std::cout << "Randomness " << args[3] << " must be numeric.\n";
            return std::nullopt;
        }
    }
    return CommandLineArgs {
        .input_folder = input_folder,
        .output_folder = output_folder,
        .num_iterations = iterations_per_sudoku,
        .randomness = randomness
    };
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::set<T>& set)
{
    if (set.empty())
        return out << "[]";
    out << "[ " << *set.begin();
    std::for_each(std::next(set.begin()), set.end(), [&out](const T& element)
    {
        out << ", " << element;
    });
    return out << " ]";
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec)
{
    if (vec.empty())
        return out << "[]";
    out << "[ " << *vec.begin();
    std::for_each(std::next(vec.begin()), vec.end(), [&out](const T& element)
    {
        out << ", " << element;
    });
    return out << " ]";
}

uint32_t generate_random_number_upto(uint32_t high) {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, high); // define the range
    int generated = distr(gen);
    return generated;
}

static std::vector<std::string> splitstring(std::string given, std::string token)
{
    std::size_t pos_start = 0, pos_end, delim_len = token.length();
    std::string chunk;
    std::vector<std::string> out;

    while ((pos_end = given.find(token, pos_start)) != std::string::npos)
    {
        chunk = given.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        out.push_back(chunk);
    }

    out.push_back(given.substr(pos_start));
    return out;
}

std::vector<fs::path> ReadDirectory(std::string directory_path)
{
    DIR *ptr = opendir(directory_path.c_str());
    if (ptr == NULL)
    {
        printf("Could not open given directory: %s", directory_path.c_str());
        return {};
    }

    struct dirent *entry;
    std::vector<fs::path> directory_contents;
    fs::path dir (directory_path);
    while ((entry = readdir(ptr)) != NULL)
    {
        if (strncmp(entry->d_name, ".", 256) == 0 || strncmp(entry->d_name, "..", 256) == 0)
        {
            continue;
        }
        fs::path file(entry->d_name);
        fs::path full_path = dir / file;
        directory_contents.push_back(full_path);
    }

    closedir(ptr);
    return directory_contents;
}

typedef struct CellCandidate {
    uint32_t row; // 0 - 3/8.
    uint32_t col; // 0 - 3/8.
    uint32_t house; // 0 - 3/8.
    std::set<uint32_t> candidates;

    bool operator==(const CellCandidate& rhs) const {
        return candidates.size() == rhs.candidates.size();
    }

    bool operator>(const CellCandidate& rhs) const {
        return candidates.size() > rhs.candidates.size();
    }

    bool operator<(const CellCandidate& rhs) const {
        return candidates.size() < rhs.candidates.size();
    }
} CellCandidate;

typedef struct CellCandidatePtr {
    CellCandidate* ptr;

    bool operator==(const CellCandidatePtr& rhs) const {
        return ptr->candidates.size() == rhs.ptr->candidates.size();
    }

    // We will add this data structure to a max-heap, so reverse the < and > meaning to get a min heap.
    bool operator>(const CellCandidatePtr& rhs) const {
        return ptr->candidates.size() < rhs.ptr->candidates.size();
    }

    bool operator<(const CellCandidatePtr& rhs) const {
        return ptr->candidates.size() > rhs.ptr->candidates.size();
    }

    CellCandidate* operator ->() const {
        return ptr;
    }
} CellCandidatePtr;

class StepDifficulty
{
    public:
        StepDifficulty(uint32_t sudoku_id): sudoku_id(sudoku_id) {}

        void RecordStep(const std::vector<CellCandidate*>& top_candidates, const CellCandidate* selected) {
            number_of_steps += 1;
            //printf("\tDebug: note a fixed cost per step %u has been added for the picked cell and sampling median\n", cost_of_doing_a_step);
            //printf("\tDebug: picked cell candidate set size is %lu\n", selected->candidates.size());
            picked_cell_candidates_sum += (selected->candidates.size() + cost_of_doing_a_step);
            for (const auto& x : top_candidates) {
                store_counts.push_back(x->candidates.size());
                //printf("\tDebug: sampling cell candidate set size is %lu\n", x->candidates.size());
            }
            top_cells_num_candidates_median_sum += (median(store_counts) + cost_of_doing_a_step);
            store_counts.clear();
        }

        void Clear() {
            number_of_steps = 0;
            picked_cell_candidates_sum = 0;
            top_cells_num_candidates_median_sum = 0;
        }

        uint32_t get_number_of_steps() const {
            return number_of_steps;
        }

        double get_average_picked_candidate() const {
            return (double)picked_cell_candidates_sum / number_of_steps;
        }

        double get_average_median_candidates() const {
            return (double)top_cells_num_candidates_median_sum / number_of_steps;
        }

    private:
        int median(std::vector<uint32_t> &v) {
            std::size_t n = v.size() / 2;
            std::nth_element(v.begin(), v.begin()+n, v.end());
            return v[n];
        }

        uint32_t sudoku_id;
        uint32_t number_of_steps = 0;
        uint32_t picked_cell_candidates_sum = 0;
        uint32_t top_cells_num_candidates_median_sum = 0;
        uint32_t cost_of_doing_a_step = 3;
        std::vector<uint32_t> store_counts;
};

class SudokuSolution
{
    public:
        SudokuSolution(uint32_t sudoku_id, const std::vector<std::vector<uint32_t>>& givens):
            sudoku_id(sudoku_id), num_successful_outcomes(0), num_failed_outcomes(0), analyzed_difficulty(StepDifficulty(sudoku_id)) {
            uint32_t num_cells = givens.size() * givens[0].size();
            uint32_t num_givens = 0;
            for (std::size_t j = 0; j < givens.size(); j++) {
                for (std::size_t k = 0; k < givens[0].size(); k++) {
                    if (givens[j][k] != 0) {
                        num_givens++;
                    }
                }
            }
            givens_to_empty_ratio = (double)num_givens/(num_cells - num_givens);
            copy_bitmap(givens, get_assignments());
        }

        uint32_t get_sudoku_id() {
            return sudoku_id;
        }

        StepDifficulty& get_difficulty() {
            return analyzed_difficulty;
        }

        void update_assignment(uint32_t row, uint32_t col, uint32_t value) {
            assignments[row][col] = value;
        }

        std::vector<std::vector<uint32_t>>& get_assignments() {
            return assignments;
        }

        std::vector<std::vector<std::vector<uint32_t>>>& get_solved_bitmaps() {
            return solved_bitmaps;
        }

        std::vector<std::string> solved_bitmaps_as_strings() {
            std::vector<std::string> bitmap_str;
            for (std::size_t i = 0; i < solved_bitmaps.size(); i++) {
                std::stringstream out;
                auto& solved_grid = solved_bitmaps[i];
                for (std::size_t j = 0; j < solved_grid.size(); j++) {
                    for (std::size_t k = 0; k < solved_grid[0].size(); k++) {
                        out << solved_grid[j][k];
                    }
                }
                bitmap_str.push_back(out.str());
            }
            return bitmap_str;
        }

        void RecordStep(const std::vector<CellCandidate*>& top_candidates, const CellCandidate* selected) {
            analyzed_difficulty.RecordStep(top_candidates, selected);
            // printf("recorded step num_steps = %u picked_candidate_size_sum = %f top_n_median_sum = %f\n", analyzed_difficulty.get_number_of_steps(), analyzed_difficulty.get_average_picked_candidate(), analyzed_difficulty.get_average_median_candidates());
        }

        void update_outcome (bool success) {
            if (success) {
                num_successful_outcomes += 1;
                solved_bitmaps.push_back({});
                auto& to_update = solved_bitmaps.back();
                copy_bitmap(get_assignments(), to_update);
            } else {
                num_failed_outcomes += 1;
            }
            // printf("outcome updated success = %u failed = %u\n", num_successful_outcomes, num_failed_outcomes);
        }

        uint32_t get_num_successful_outcomes() const {
            return num_successful_outcomes;
        }

        uint32_t get_num_failed_outcomes() const {
            return num_failed_outcomes;
        }

        double get_givens_to_empty_ratio() const {
            return givens_to_empty_ratio;
        }

    private:
        void copy_bitmap(const std::vector<std::vector<uint32_t>>& givens, std::vector<std::vector<uint32_t>>& to_update) {
            for (uint32_t i = 0; i < givens.size(); i++) {
                to_update.push_back({});
                for (uint32_t j = 0; j < givens[0].size(); j++) {
                    to_update[i].push_back(givens[i][j]);
                }
            }
        }

        uint32_t sudoku_id;
        double givens_to_empty_ratio;
        std::vector<std::vector<uint32_t>> assignments;
        std::vector<std::vector<std::vector<uint32_t>>> solved_bitmaps;
        uint32_t num_successful_outcomes;
        uint32_t num_failed_outcomes;
        StepDifficulty analyzed_difficulty;
};

typedef struct Point {
  uint32_t row;
  uint32_t col;

  bool operator==(const Point& rhs) const {
    return this->row == rhs.row && this->col == rhs.col;
  }

  bool operator>(const Point& rhs) const {
    return this->row > rhs.row ||
    (this->row == rhs.row && this->col > rhs.col);
  }

  bool operator<(const Point& rhs) const {
    return rhs > *this;
  }
 
  friend std::ostream &operator<<(std::ostream &os, const Point &p);
} Point;

std::ostream &operator<<(std::ostream &os, const Point &p)
{
    os << p.row << " " << p.col;
    return os;
}

void print_sudoku_assignments(const std::vector<std::vector<uint32_t>>& assignments) {
    for (uint32_t row = 0; row < assignments.size(); row ++) {
        for (uint32_t col = 0; col < assignments[0].size(); col ++) {
            printf("%d\t", assignments[row][col]);
        }
        printf("\n");
    }
    printf("\n");
}

std::set<uint32_t> givens_in_row(uint32_t row, std::vector<std::vector<uint32_t>> givens) {
    std::set<uint32_t> givens_in_row {};
    for (uint32_t col = 0; col < givens[0].size(); col ++) {
        if (givens[row][col]) {
            givens_in_row.insert(givens[row][col]);
        }
    }
    return givens_in_row;
}

std::set<uint32_t> givens_in_col(uint32_t col, std::vector<std::vector<uint32_t>> givens) {
    std::set<uint32_t> givens_in_col {};
    for (uint32_t row = 0; row < givens.size(); row ++) {
        if (givens[row][col]) {
            givens_in_col.insert(givens[row][col]);
        }
    }
    return givens_in_col;
}

bool is_in_bounding_house(uint32_t row, uint32_t col, Point tl, Point br) {
    if (row >= tl.row && col >= tl.col) {
    if (row <= br.row && col <= br.col) {
        return true;
    }
    }
    return false;
}

uint32_t get_house_for_cell(uint32_t row, uint32_t col, std::map<uint32_t, std::pair<Point, Point>>& houses) {
    for (auto it = houses.begin(); it != houses.end(); it++) {
    if (is_in_bounding_house(row, col, it->second.first, it->second.second)) {
        return it->first;
    }
    }
    return -1;
}

std::set<Point> points_in_house(uint32_t row, uint32_t col, std::map<uint32_t, std::pair<Point, Point>>& houses) {
    /*for (auto& x : houses) {
        std::cout << "house " << x.first << " co-ordinates, top left: " << x.second.first << " bottom right: " << x.second.second << std::endl;
    }*/
    uint32_t house = get_house_for_cell(row, col, houses);
    //std::cout << "house for cell " << row << " " << col << " is " << house << std::endl;
    auto& [top_left, bottom_right] = houses[house];
    //std::cout << "top left " << top_left << " br " << bottom_right << std::endl;

    std::set<Point> points;
    for (uint32_t i = top_left.row; i <= bottom_right.row; i++) {
    for (uint32_t j = top_left.col; j <= bottom_right.col; j++) {
        points.insert(Point{.row = i, .col = j});
    }}

    return points;
}

std::set<uint32_t> givens_in_house(uint32_t row, uint32_t col, std::vector<std::vector<uint32_t>> givens, std::map<uint32_t, std::pair<Point, Point>> houses) {
    auto points_house = points_in_house(row, col, houses);

    std::set<uint32_t> givens_in_house;
    for (uint32_t i = 0; i < points_house.size(); i++) {
        for (auto& p : points_house) {
            if (givens[p.row][p.col]) {
                givens_in_house.insert(givens[p.row][p.col]);
            }
        }
    }
    return givens_in_house;
}


std::set<Point> get_adjacent_points(const Point& coord, std::map<uint32_t, std::pair<Point, Point>> houses) {
    std::set<Point> adj_cols, adj_rows, adj_houses;
    uint32_t N =  houses.size();

    for (uint32_t col = 0; col < N; col ++) {
        adj_cols.insert(Point{.row = coord.row, .col = col});
    }
    for (uint32_t row = 0; row < N; row ++) {
        adj_rows.insert(Point{.row = row, .col = coord.col});
    }
    adj_houses = points_in_house(coord.row, coord.col, houses);

    std::set<Point> adj_points;
    adj_points.insert(adj_rows.begin(), adj_rows.end());
    adj_points.insert(adj_cols.begin(), adj_cols.end());
    adj_points.insert(adj_houses.begin(), adj_houses.end());
    return adj_points;
}

template <typename T>
void print_set(std::set<T> to_print) {
    for (auto& x : to_print) {
        std::cout << x << ", ";
    }
    std::cout << std::endl;
}

CellCandidate CreateCellCandidate(uint32_t row, uint32_t col, std::set<uint32_t> possibles, std::map<uint32_t, std::pair<Point, Point>> houses) {
    uint32_t house = get_house_for_cell(row, col, houses);
    return {
        .row = row,
        .col = col,
        .house = house,
        .candidates = std::move(possibles)
    };
}

std::map<Point, CellCandidate>
ComputeInitialCandidates(uint32_t rank,
                         std::vector<std::vector<uint32_t>> givens,
                         std::map<uint32_t, std::pair<Point, Point>> houses) {
    uint32_t rows = rank * rank;
    uint32_t cols = rank * rank;
    std::map<Point, CellCandidate> candidates;
    std::set<uint32_t> possible_assignments = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t j = 0; j < cols; j++) {
            //std::cout << "Computing candidates for [" << i << ", " << j << "]" << std::endl;
            if (givens[i][j] != 0) {
                // If a value is given for a cell, we should not have to consider the candidates for that cell.
                // candidates.insert({Point {.row = i, .col = j}, CreateCellCandidate(i, j, std::set<uint32_t>{givens[i][j]}, houses)});
                continue;
            }
            std::set<uint32_t> all_givens = givens_in_row(i, givens);
            //std::cout << "Givens in row are ";
            //print_set(all_givens);
            std::set<uint32_t> givens_col = givens_in_col(j, givens);
            //std::cout << "Givens in column are ";
            //print_set(givens_col);
            all_givens.insert(givens_col.begin(), givens_col.end());
            std::set<uint32_t> givens_house = givens_in_house(i, j, givens, houses);
            //std::cout << "Givens in house are ";
            //print_set(givens_house);
            all_givens.insert(givens_house.begin(), givens_house.end());

            std::set<uint32_t> possibles;
            std::set_difference(possible_assignments.begin(), possible_assignments.end(),
                                all_givens.begin(), all_givens.end(), std::inserter(possibles, possibles.end()));

            candidates.insert({Point {.row = i, .col = j}, CreateCellCandidate(i, j, std::move(possibles), houses)});
        }
    }

    return candidates;
}

// This class is thread-compatible.
// It is only ever used by one thread.
class ParsedSudoku
{
public:
    uint32_t get_sudoku_id()
    {
        return sudoku_id;
    }
    std::string get_problem_file_path()
    {
        return problem_file_path;
    }
    std::string get_compressed_repr()
    {
        return compressed_repr;
    }
    uint32_t get_rank()
    {
        return rank;
    }
    std::vector<std::vector<uint32_t>> get_givens()
    {
        return givens;
    }
    std::map<uint32_t, std::pair<Point, Point>> get_houses()
    {
        return houses;
    }

    // Delete default ctor.
    ParsedSudoku() = delete;

    static std::optional<std::unique_ptr<ParsedSudoku>> ParsedSudokuFactory(fs::path problem_file_path);

private:
    ParsedSudoku(
        uint32_t sudoku_id, std::string problem_file_path,
        std::string compressed_repr, uint32_t rank,
        std::vector<std::vector<uint32_t>> givens,
        std::map<uint32_t, std::pair<Point, Point>> houses)
    {
        this->sudoku_id = sudoku_id;
        this->problem_file_path = problem_file_path;
        this->compressed_repr = compressed_repr;
        this->rank = rank;
        this->givens = std::move(givens);
        this->houses = std::move(houses);
    }

    uint32_t sudoku_id;
    std::string problem_file_path;
    std::string compressed_repr;
    uint32_t rank;
    std::map<uint32_t, std::pair<Point, Point>> houses;
    std::vector<std::vector<uint32_t>> givens;
};

// Factory to generate a fully initialized ParsedSudoku.
std::optional<std::unique_ptr<ParsedSudoku>> ParsedSudoku::ParsedSudokuFactory(fs::path problem_file_path)
{
    std::string file_contents;
    std::ifstream sudoku_fs(problem_file_path);
    uint32_t sudoku_id;
    if (sudoku_fs.is_open())
    {
        std::getline(sudoku_fs, file_contents);
        sudoku_id = std::stoi(problem_file_path.stem());
    }
    else
    {
        std::cout << "Failed to open sudoku problem file " << problem_file_path << std::endl;
        return std::nullopt;
    }

    // Givens is one std::string of comma separated:
    // valuerowcolumn std::strings.
    // First find the Rank of the sudoku by counting the
    // comma separated values.
    // Then populate the candidate std::vector.
    std::vector<std::string> valuerowcolumn = splitstring(file_contents, ",");
    // 17 is the minimum number of givens for a 3X3 Sudoku with a unique solution.
    uint32_t rank = valuerowcolumn.size() >= 17 ? 3 : 2;

    std::vector<std::vector<uint32_t>> givens = {};
    for (uint32_t i = 0; i < rank * rank; i++)
    {
        givens.push_back({});
        for (uint32_t j = 0; j < rank * rank; j++)
        {
            givens[i].push_back(0);
        }
    }

    for (std::string vrc : valuerowcolumn)
    {
        uint32_t val = vrc[0] - '0';
        uint32_t row = vrc[1] - '0';
        uint32_t col = vrc[2] - '0';
        if (row >= givens.size() || col >= givens[0].size()) {
            std::cout << "Invalid row and column csv " << row << col << std::endl;
            continue;
        }
        givens[row][col] = val;
    }

    //print_sudoku_assignments(givens);

    std::map<uint32_t, std::pair<Point, Point>> houses;
    uint32_t house_nr = 0;
    while (house_nr < rank * rank) {
        uint32_t row = house_nr;
        uint32_t col = 0;
        while (col < rank * rank) {
            Point top_left = Point {.row = row, .col = col};
            Point bottom_right = Point {.row = row + rank - 1, .col = col + rank - 1};
            houses.insert({house_nr, std::make_pair<Point, Point>(std::move(top_left), std::move(bottom_right)) });
            house_nr ++;
            col += rank;
        }
    }

    return std::make_unique<ParsedSudoku>(
        ParsedSudoku(sudoku_id, problem_file_path, file_contents, rank, std::move(givens),
        std::move(houses)));
}

class SudokuAnalysis {
    public:
        SudokuAnalysis(std::vector<SudokuSolution>&& solutions) : N(solutions.size()), solutions(std::move(solutions)) {
            SudokuAnalysis::Compute();
        }

        std::string serialize_to_json_string() {
            std::stringstream o;
            o << "{";
            json_add_field("id", sudoku_id, o);
            json_add_field("givens_to_empty_ratio", givens_to_empty_ratio, o);
            json_add_field("unique_solutions", all_unique_solutions, o);
            json_add_field("avg_number_of_steps", avg_number_of_steps, o);
            json_add_field("avg_failure_to_success_ratio", avg_failure_to_success_ratio, o);
            json_add_field("avg_candidates_for_picked_cell", avg_candidates_for_picked_cell, o);
            json_add_field("avg_median_of_candidates_for_cell", avg_median_of_candidates_for_cell, o);
            o << "}";
            return o.str();
        }

        void DebugAnalysisData() const {
            std::cout << "Sudoku ID: " << sudoku_id << std::endl;
            std::cout << "All unique solutions : " << all_unique_solutions << std::endl;
            std::cout << "Number of steps : " << number_of_steps << std::endl;
            std::cout << "Failure/Success ratio : " << failure_to_success_ratios << std::endl;
            std::cout << "Candidates for picked cell : " << candidates_for_picked_cell << std::endl;
            std::cout << "Median of candidates for cell : " << median_of_candidates_for_cell << std::endl;
        }
    protected:
        std::set<std::string> unique_solutions_in_iteration(SudokuSolution& s) {
            std::vector<std::string> bitmap_str = s.solved_bitmaps_as_strings();
            std::set<std::string> unique_solutions;
            for (std::string& el : bitmap_str) {
                unique_solutions.insert(el);
            }
            return unique_solutions;
        }
        void Compute() {
            for (std::size_t i = 0; i < N; i++) {
                auto& s = solutions[i];
                sudoku_id = s.get_sudoku_id();
                givens_to_empty_ratio = s.get_givens_to_empty_ratio();
                number_of_steps.push_back(s.get_difficulty().get_number_of_steps());
                //double f2s_ratio = s.get_num_failed_outcomes()/s.get_num_successful_outcomes();
                //failure_to_success_ratios.push_back(f2s_ratio);
                candidates_for_picked_cell.push_back(s.get_difficulty().get_average_picked_candidate());
                median_of_candidates_for_cell.push_back(s.get_difficulty().get_average_median_candidates());
                auto it_unique = unique_solutions_in_iteration(s);
                all_unique_solutions.insert(it_unique.begin(), it_unique.end());
            }

            double total_number_of_steps = std::accumulate(number_of_steps.begin(), number_of_steps.end(), 0.0);
            avg_number_of_steps = total_number_of_steps/N;
            double total_f2s_ratio = std::accumulate(failure_to_success_ratios.begin(), failure_to_success_ratios.end(), 0.0);
            avg_failure_to_success_ratio = total_f2s_ratio/N;
            double total_picked_cell = std::accumulate(candidates_for_picked_cell.begin(), candidates_for_picked_cell.end(), 0.0);
            avg_candidates_for_picked_cell = total_picked_cell/N;
            double total_median_candidates = std::accumulate(median_of_candidates_for_cell.begin(), median_of_candidates_for_cell.end(), 0.0);
            avg_median_of_candidates_for_cell = total_median_candidates/N;
        }

        void json_add_field(std::string_view key, auto& value, std::stringstream& s) {
            s << key << " : " << value << ",\n";
        }
        std::set<std::string>& get_all_unique_solutions() {
            return all_unique_solutions;
        }
        double get_givens_to_empty_ratio()  {
            return givens_to_empty_ratio;
        }
        std::vector<double> get_number_of_steps()  {
            return number_of_steps;
        }
        double get_avg_number_of_steps() {
            return avg_number_of_steps;
        }
        std::vector<double> get_success_failure_ratios()  {
            return failure_to_success_ratios;
        }
        double get_avg_success_failure_ratio() {
            return avg_failure_to_success_ratio;
        }
        std::vector<double> get_candidates_for_picked_cell() {
            return candidates_for_picked_cell;
        }
        double get_avg_candidates_for_picked_cell() {
            return avg_candidates_for_picked_cell;
        }
        std::vector<double> get_median_of_candidates_for_cell() {
            return median_of_candidates_for_cell;
        }
        double get_avg_median_of_candidates_for_cell() {
            return avg_median_of_candidates_for_cell;
        }
        std::vector<SudokuSolution>& get_all_solutions() {
            return solutions;
        }

    private:
        std::size_t N;
        uint32_t sudoku_id;
        std::set<std::string> all_unique_solutions;
        std::vector<double> number_of_steps;
        double avg_number_of_steps;
        double givens_to_empty_ratio;
        std::vector<double> failure_to_success_ratios;
        double avg_failure_to_success_ratio;
        std::vector<double> candidates_for_picked_cell;
        double avg_candidates_for_picked_cell;
        std::vector<double> median_of_candidates_for_cell;
        double avg_median_of_candidates_for_cell;
        std::vector<SudokuSolution> solutions;
};

class FancyDecorator : public SudokuAnalysis {
    public:
        FancyDecorator(std::vector<SudokuSolution>&& solutions) : SudokuAnalysis(std::move(solutions)) {}
        
        template <typename ... Args>
        auto call(std::string method_name, Args&& ... args) {
            /*
            // This did not work as expected, as auto needs one resolved datatype.
            if (!computed) {            
                SudokuAnalysis::Compute();
                computed = true;
            }

            switch (method_name) {
                case "get_has_unique_solution":
                    return SudokuAnalysis::get_has_unique_solution(std::forward<Args>(args)...);
                    break;
                case "get_all_unique_solutions":
                    return SudokuAnalysis::get_all_unique_solutions(std::forward<Args>(args)...);
                    break;
                case "get_success_failure_ratios":
                    return SudokuAnalysis::get_success_failure_ratios(std::forward<Args>(args)...);
                    break;
                case "get_avg_success_failure_ratio":
                    return SudokuAnalysis::get_avg_success_failure_ratio(std::forward<Args>(args)...);
                    break;
                case "get_candidates_for_picked_cell":
                    return SudokuAnalysis::get_candidates_for_picked_cell(std::forward<Args>(args)...);
                    break;
                case "get_avg_candidates_for_picked_cell":
                    return SudokuAnalysis::get_avg_candidates_for_picked_cell(std::forward<Args>(args)...);
                    break;
                case "get_median_of_candidates_for_cell":
                    return SudokuAnalysis::get_median_of_candidates_for_cell(std::forward<Args>(args)...);
                    break;
                case "get_avg_median_of_candidates_for_cell":
                    return SudokuAnalysis::get_avg_median_of_candidates_for_cell(std::forward<Args>(args)...);
                    break;
                case "get_givens_to_empty_ratio":
                    return SudokuAnalysis::get_givens_to_empty_ratio(std::forward<Args>(args)...);
                    break;
                case "serialize_to_json_string":
                    return SudokuAnalysis::serialize_to_json_string(std::forward<Args>(args)...);
                    break;
                default:
                    printf("Invalid method name : %s\n", method_name.c_str());
                    break;
            }*/
        }

    private:
        bool computed = false;
};

class Solver
{
    // Strategy
    // 1. Given populated candidates[][], sort the bags by lowest number of candidates.
    //    For the bags, one by one try each equal choice.
    //    If the attempt eventually succeeds, add to the steps_to_success std::vector.
    //    Else add to steps_to_failure std::vector.

    //  2. Update, re-analyze and re-sort candidates[][].
public:
    // Runs in a single thread.
    void SolveSudoku() {
        //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        //std::cout<< "Thread ID is " << std::this_thread::get_id() << std::endl;
        for (uint32_t i = 0; i < num_iterations; i++) {
            std::cout << "Iteration number " << i << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
            SolveSudokuAndRecurse(solver_state.get());
            //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            //std::cout << "Time taken = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
            found_solutions.push_back(solver_state->get_solution());
            solver_state->Reset();
        }

        AnalyzeSolutions();
    }

    Solver(std::unique_ptr<ParsedSudoku> sudoku, fs::path output_folder,
           uint32_t num_iterations, uint32_t randomness)
    {
        this->sudoku = std::move(sudoku);
        this->solution_file_path = output_folder;
        this->num_iterations = num_iterations;
        this->randomness = randomness;
        printf("Initialized solver with o/p %s, iter %u, randomness %u\n", this->solution_file_path.c_str(), this->num_iterations, this->randomness);
        BuildInitialSudokuState();
    }

    Solver() = delete;

private:
    class SudokuStepState {
        public:
            std::vector<std::vector<uint32_t>>& get_assignments() {
                return solution.get_assignments();
            }

            SudokuSolution& get_solution() {
                return solution;
            }

            StepDifficulty& get_difficulty_metrics() {
                return solution.get_difficulty();
            }

            void Add(CellCandidate* cell, uint32_t value) {
               solution.update_assignment(cell->row, cell->col, value);
               UpdateAdjacentHouses(cell, value);
            }

            void Subtract(CellCandidate* cell, uint32_t value) {
                solution.update_assignment(cell->row, cell->col, 0);
                UndoUpdateToAdjacentHouses(cell, value);
            }

            CellCandidate* PopAndGetOptimalNextCell() {
                if (pq.empty()) {
                    //printf("pq is empty\n");
                    //print_sudoku_assignments(get_assignments());
                    return nullptr;
                }
                uint32_t shuffle_max = randomness > pq.size() ? pq.size() : randomness;
                return ShuffleCandidatesAndRecordStep(shuffle_max);
            }

            void PushToPQ(CellCandidate* cell) {
                //printf("pushed [%u, %u]\n", cell->row, cell->col);
                pq.push(CellCandidatePtr{.ptr = cell});
            }

            void Reset() {
                if (!revert_updates.empty()) {
                    printf("Something is wrong, revert_updates should be empty before Reset()\n");
                }
                revert_updates.clear();
                BuildPriorityQueue();
                solution = SudokuSolution(SudokuSolution(sudoku_ptr->get_sudoku_id(), sudoku_ptr->get_givens()));
            }

            SudokuStepState(std::map<Point, CellCandidate> candidates,
                            ParsedSudoku* sudoku, uint32_t randomness) : solution(SudokuSolution(sudoku->get_sudoku_id(), sudoku->get_givens())) {
                this->randomness = randomness;
                pq_map = std::move(candidates);
                BuildPriorityQueue();
                sudoku_ptr = sudoku;
                for (auto it = pq_map.begin(); it != pq_map.end(); it++) {
                    adjacent_points.insert({it->first, get_adjacent_points(it->first, sudoku_ptr->get_houses())});
                }
            }

            SudokuStepState() = delete;

            void print_priority_queue(int N) {
                std::vector<CellCandidatePtr> popped;
                int count = 0;
                if (N == -1) {
                    N = pq.size();
                }
                while (count < N) {
                    auto& top = pq.top();
                    printf("[%u, %u] : ptr: %p size = %lu\n", top.ptr->row, top.ptr->col, (void*)top.ptr, top.ptr->candidates.size());
                    popped.push_back(top);
                    pq.pop();
                    count++;
                }
                for (uint32_t i = 0; i < popped.size(); i++) {
                    for (uint32_t j = i + 1; j < popped.size(); j++) {
                        if (popped[i].ptr == popped[j].ptr) {
                            printf("Heap corruption, element [%u, %u] is on the heap twice.\n", popped[i].ptr->row, popped[i].ptr->col);
                        }
                    }
                }
                for (auto& elem : popped) {
                    pq.push(CellCandidatePtr{.ptr = elem.ptr});
                }
                // std::cout << std::endl;
            }

            // Do not ask why I have to write this. Hint: gdb.
            void __attribute__ ((noinline)) print_map() {
                for (auto it = pq_map.begin(); it != pq_map.end(); it++) {
                    const Point& p = it->first;
                    CellCandidate* cell = &it->second;
                    std::cout << "Key : " << p << " Cell Ptr = " << cell << " Candidates = ";
                    print_set<uint32_t>(cell->candidates);
                }
            }
        private:
            void BuildPriorityQueue() {
                pq = std::priority_queue<CellCandidatePtr>();
                for (auto it = pq_map.begin(); it != pq_map.end(); it++) {
                    CellCandidate* cell = &it->second;
                    //std::cout<< "In BuildPriorityQueue for point " << it->first << " Cell pointer is " << cell << std::endl;
                    pq.push(CellCandidatePtr{.ptr = cell});
                }
                //print_priority_queue(pq.size());
            }
            CellCandidate* ShuffleCandidatesAndRecordStep(uint32_t shuffle_max);
            void UndoUpdateToAdjacentHouses(CellCandidate* cell, uint32_t value_to_be_readded);
            void UpdateAdjacentHouses(CellCandidate* cell, uint32_t value);

            typedef struct Update {
                uint32_t value;
                Point updated;
                bool operator==(const Update& rhs) const {
                    return this->value == rhs.value && this->updated == rhs.updated;
                }
                bool operator>(const Update& rhs) const {
                    return this->value > rhs.value || (this->value == rhs.value && this->updated > rhs.updated);
                }
                bool operator<(const Update& rhs) const {
                    return rhs > *this;
                }
            } Update;
            // For any step, keeps track of the cells updated by a new assignment.
            // Also keep track of who updated which neighbor, and what that update was,
            // so we do not have unexpected updates.
            std::map<Point, std::set<Update>> revert_updates;
            ParsedSudoku* sudoku_ptr;
            // Authoritative owner of all CellCandidate(s).
            std::map<Point, CellCandidate> pq_map;
            std::priority_queue<CellCandidatePtr> pq;
            // Current assignments to the coordinates of the sudoku.
            std::vector<std::vector<uint32_t>> assignments;
            // Adjacent points of a given point.
            std::map<Point, std::set<Point>> adjacent_points;
            // Associated solution object.
            SudokuSolution solution;
            uint32_t randomness;
    };

    void SolveSudokuAndRecurse(SudokuStepState* state);
    void AnalyzeSolutions();
    bool WriteToOutputFile(std::string json_stats);
    void BuildInitialSudokuState() {
        auto cell_candidates = ComputeInitialCandidates(sudoku->get_rank(), sudoku->get_givens(), sudoku->get_houses());
        solver_state = std::make_unique<SudokuStepState>(std::move(cell_candidates), sudoku.get(), randomness);
    }

    std::unique_ptr<ParsedSudoku> sudoku;
    fs::path solution_file_path;
    uint32_t num_iterations;
    uint32_t randomness;
    std::unique_ptr<SudokuStepState> solver_state;
    std::vector<SudokuSolution> found_solutions;
};

void Solver::AnalyzeSolutions() {
    std::unique_ptr<SudokuAnalysis> analyzer = std::make_unique<SudokuAnalysis>(std::move(found_solutions));
    analyzer->DebugAnalysisData();
    auto json = analyzer->serialize_to_json_string();
    std::cout << "Solver statistics as JSON: " << json << std::endl;
    WriteToOutputFile(json);
}

bool Solver::WriteToOutputFile(std::string json_stats) {
    std::ofstream output_file;
    fs::path out_file (std::to_string(sudoku->get_sudoku_id()) + ".json");
    fs::path output_path = solution_file_path / out_file;
    output_file.open(output_path);
    output_file << json_stats;
    output_file.close();
    return true;
}

CellCandidate* Solver::SudokuStepState::ShuffleCandidatesAndRecordStep(uint32_t shuffle_max) {
    std::vector<CellCandidate*> top_candidates;
    top_candidates.reserve(shuffle_max);
    int cell_with_no_candidates = -1;

    for(uint32_t i = 0; i < shuffle_max; i++) {
        CellCandidatePtr el = pq.top();
        if (el->candidates.size() == 0) {
            cell_with_no_candidates = i;
        }
        top_candidates.push_back(el.ptr);
        pq.pop();
    }

    CellCandidate* top = nullptr;
    if (cell_with_no_candidates == -1) {
        uint32_t shuffle_to_idx = shuffle_max > 0 ? shuffle_max - 1 : shuffle_max;
        uint32_t shuffle_pick = generate_random_number_upto(shuffle_to_idx);
        top = top_candidates[shuffle_pick];
    } else {
        top = top_candidates[cell_with_no_candidates];
    }

    //std::cout << "Top cell has ptr " << top << std::endl;
    solution.RecordStep(top_candidates, top);

    for(uint32_t i = 0; i < shuffle_max; i++) {
        if (top_candidates[i] != top) {
            pq.push(CellCandidatePtr(top_candidates[i]));
        }
    }

    return top;
}

void Solver::SudokuStepState::UndoUpdateToAdjacentHouses(CellCandidate* cell, uint32_t value_to_be_readded) {
    Point cell_coord = Point({.row = cell->row, .col = cell->col});
    //std::cout << "Reverting everything for " << cell_coord << std::endl;
    auto it = revert_updates.find(cell_coord);
    if (it != revert_updates.end()) {
        auto& [_, re_add_set] = *it;
        // What were the updates made by cell at cell_coord, and to which points.
        for (auto& [x_value, x_coord] : re_add_set) {
            //std::cout<< "CellCandidate " << cell_coord << " was reverted " << value_to_be_readded << " and adjacent cell " << x_coord << " was updated." << std::endl;
            if (!pq_map.contains(x_coord)) {
               continue;
            }
            CellCandidate& x_cell = pq_map[x_coord];
            if (x_value != value_to_be_readded) {
                std::cout << "ERROR: Unexpected, changed value " << x_value << " should match " << value_to_be_readded << std::endl;
            }
            x_cell.candidates.insert(value_to_be_readded);
        }
    }
    revert_updates.erase(cell_coord);
}

void Solver::SudokuStepState::UpdateAdjacentHouses(CellCandidate* cell, uint32_t value) {
    Point cell_coord = Point{.row = cell->row, .col = cell->col};
    if (revert_updates.contains(cell_coord)) {
        std::cout << "ERROR: Cell at " << cell_coord << " has already been traversed in DFS" << std::endl;
    }

    // This can probably be pre-computed.
    std::set<Point>& adj_points = adjacent_points[cell_coord];
    //auto it = revert_updates.end();

    for (auto& update_coord : adj_points) {
        if (!pq_map.contains(update_coord)) {
            continue;
        }
        auto& update_cell = pq_map[update_coord];
        // Assigned value was a candidate for this cell.
        if (update_cell.candidates.contains(value)) {
            update_cell.candidates.erase(value);
            Update to_insert = Update {.value = value, .updated = update_coord};
            auto it = revert_updates.find(cell_coord);
            //std::cout<< "CellCandidate " << cell_coord << " was assigned " << value << " and adjacent cell " << update_coord << " needs to be updated." << std::endl;
            // Add print statements here.
            if (it == revert_updates.end()) {
                //std::cout<< "CellCandidate " << cell_coord << " was assigned " << value << " and adjacent cell " << update_coord << " was updated." << std::endl;
               it = revert_updates.insert(revert_updates.end(), {cell_coord, std::set<Update>({to_insert})});
            } else {
                //std::cout<< "CellCandidate " << cell_coord << " was assigned " << value << " and adjacent cell " << update_coord << " was updated. it->key = " << it->first << std::endl;
                auto& [_, update_set] = *it;
                update_set.insert(to_insert);
            }
        }
    }

    /*auto it = revert_updates.find(cell_coord);
    if (it != revert_updates.end()) {
        std::cout << "CellCandidate update_set contains: " << std::endl;
        auto& [_, update_set] = *it;
        for (auto& [val, u] : update_set) {
            std::cout << val << ", " << u << std::endl;
        }
        std::cout<< std::endl;
    }*/
}

void Solver::SolveSudokuAndRecurse(SudokuStepState* state) {
    CellCandidate* optimal_next_cell = state->PopAndGetOptimalNextCell();
    //state->print_map();
    if (optimal_next_cell == nullptr) {
        //std::cout << "no optimal next cell\n";
        state->get_solution().update_outcome(true);
        return;
    }
    //printf("Next cell from heap is [%u, %u]\n", optimal_next_cell->row, optimal_next_cell->col);
    if (optimal_next_cell->candidates.size() == 0) {
        // We made an incorrect assignment.
        //std::cout<< "No options for cell " << optimal_next_cell->row << " " << optimal_next_cell->col << std::endl;
        state->PushToPQ(optimal_next_cell);
        state->get_solution().update_outcome(false);
        return;
    }

    std::set<uint32_t> candidates_to_pick = std::set<uint32_t>(optimal_next_cell->candidates);
    /*std::cout<< "\nCandidates for [" << optimal_next_cell->row << ", " << optimal_next_cell->col << "]" << " are ";
    for (uint32_t c : candidates_to_pick) {
        std::cout << c << ",";
    }
    std::cout << std::endl;*/

    if (state->get_assignments()[optimal_next_cell->row][optimal_next_cell->col] != 0) {
        std::cout<< "ERROR: this cell is assigned, already" << std::endl;
        state->PushToPQ(optimal_next_cell);
        state->get_solution().update_outcome(false);
        return;
    }

    // Pick a candidate coloring for the cell.
    for (uint32_t picked_candidate : candidates_to_pick) {
        //std::cout<< "Assigning " << picked_candidate << " to " << "[" << optimal_next_cell->row << ", " << optimal_next_cell->col << "]" << std::endl;
        //printf("[%u, %u] Before adding\n", copy->row, copy->col);
        //state->print_priority_queue(0);
        state->Add(optimal_next_cell, picked_candidate);
        //printf("[%u, %u]Entering recursion\n", optimal_next_cell->row, optimal_next_cell->col);
        //state->print_priority_queue(5);
        SolveSudokuAndRecurse(state);
        // Undo previous update.
        //printf("[%u, %u]After recursion ended, and we need to start reverting\n", optimal_next_cell->row, optimal_next_cell->col);
        //state->print_priority_queue(5);
        state->Subtract(optimal_next_cell, picked_candidate);
        //printf("[%u, %u]After reverting\n", copy->row, copy->col);
        // state->print_priority_queue(-1);
    }

    // We could not find a correct assignment, which means we want to recurse upwards.
    state->PushToPQ(optimal_next_cell);
}

std::vector<std::unique_ptr<ParsedSudoku>> ParseSudokuProblems(std::vector<fs::path> &sudoku_problem_files)
{
    std::vector<std::unique_ptr<ParsedSudoku>> parsed_sudokus = {};
    for (auto &file : sudoku_problem_files)
    {
        // std::cout << "passing file " << file << std::endl;
        auto maybe_parsed_file = ParsedSudoku::ParsedSudokuFactory(file);
        if (maybe_parsed_file != std::nullopt)
        {
            parsed_sudokus.push_back(std::move(*maybe_parsed_file));
        }
    }
    return parsed_sudokus;
}

void analyzeSudokus(CommandLineArgs options)
{
    std::vector<fs::path> sudoku_pathnames = ReadDirectory(options.input_folder);
    std::vector<std::unique_ptr<ParsedSudoku>> parsed_sudokus = ParseSudokuProblems(sudoku_pathnames);

    std::vector<Solver> solvers;
    auto ps_it = parsed_sudokus.begin();
    uint32_t i = 0;
    while (i < 1/*parsed_sudokus.size()*/) {
        auto next_ps_it = ps_it + 1;
        solvers.push_back(Solver(std::move(*ps_it), fs::path(options.output_folder), options.num_iterations, options.randomness));
        ps_it = next_ps_it;
        i++;
    }

    uint32_t concurrent_threads = std::thread::hardware_concurrency();
    auto sol_it = solvers.begin();
    while (sol_it != solvers.end()) {
        uint32_t running_threads = 0;
        std::vector<std::thread> sudoku_threads;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        // For each Sudoku, we create a thread to solve and publish the result of the sudoku
        // as a CSV file into the output folder.
        while (running_threads < concurrent_threads && sol_it != solvers.end()) {
            sudoku_threads.push_back(std::thread(&Solver::SolveSudoku, &(*sol_it)));
            sol_it++;
            running_threads++;
        }
        // Wait for all threads to finish.
        for (auto& thread : sudoku_threads) {
            while (!thread.joinable()) {}
            thread.join();
        }
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time taken = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms" << std::endl;
        std::cout << "Sudokus remaining: " << (solvers.end() - sol_it) << std::endl;
        // Free up memory we are not going to use.
        sudoku_threads.clear();
        solvers.erase(solvers.begin(), solvers.begin() + running_threads);
        sol_it = solvers.begin();
        running_threads = 0;
    }

    std::cout << "back to the main thread with ID " << std::this_thread::get_id() << " threads joined" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc < 3 || argc > 5) {
        std::cerr << "Please provide necessary command line values.\n\n";
        std::cerr << "usage: ./analyzer <input folder> <output folder> <iterations>(optional) <randomness>(optional) \n";
        std::cerr << "Example: ./analyzer input_sudokus_folder/ sudoku_analyzis_folder/ 1 0\n\n";
        return EXIT_FAILURE;
    }
    const std::vector<std::string> args(argv + 1, argv + argc);
    auto parsed_args = parse_command_line_args(args);
    if (parsed_args.has_value()) {
        analyzeSudokus(*parsed_args);
    }
}
