/*
 * C++ Program to solve a sudoku and analyse:
 *   1. Whether the solution is unique.
 *   2. How difficult is the sudoku.
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
namespace fs = std::filesystem;

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

typedef struct StepDifficulty
{
    int sudoku_id;
    int num_of_easy_cells;
    int candidates_per_easy_cell;
    int average_number_of_candidates_per_cell;
} StepDifficulty;

typedef struct Difficulty
{
    int sudoku_id;
    std::vector<StepDifficulty> steps_to_success;
    std::vector<StepDifficulty> steps_to_failure;
} Difficulty;

typedef struct SudokuSolution
{
    int sudoku_id;
    std::vector<std::vector<int>> assignments;
    Difficulty analysed_difficulty;
} SudokuSolution;

typedef struct CellCandidate {
    int row; // 0 - 3/8.
    int col; // 0 - 3/8.
    int house; // 0 - 3/8.
    std::set<int> candidates;

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

typedef struct Point {
  int row;
  int col;

  bool operator==(const Point& rhs) const {
    return this->row == rhs.row && this->col == rhs.col;
  }

  bool operator>(const Point& rhs) const {
    return this->row > rhs.row ||
    (this->row == rhs.row && this->col >= rhs.col);
  }

  bool operator<(const Point& rhs) const {
    return rhs > *this && *this != rhs;
  }
 
  friend std::ostream &operator<<(std::ostream &os, const Point &p);
} Point;

std::ostream &operator<<(std::ostream &os, const Point &p)
{
    os << p.row << " " << p.col;
    return os;
}

void print_sudoku_assignments(std::vector<std::vector<int>>& assignments) {
    for (uint32_t row = 0; row < assignments.size(); row ++) {
        for (uint32_t col = 0; col < assignments[0].size(); col ++) {
            printf("%d\t", assignments[row][col]);
        }
        printf("\n");
    }
}

std::set<int> givens_in_row(int row, std::vector<std::vector<int>> givens) {
    std::set<int> givens_in_row {};
    for (uint32_t col = 0; col < givens[0].size(); col ++) {
        if (givens[row][col]) {
            givens_in_row.insert(givens[row][col]);
        }
    }
    return givens_in_row;
}

std::set<int> givens_in_col(int col, std::vector<std::vector<int>> givens) {
    std::set<int> givens_in_col {};
    for (uint32_t row = 0; row < givens.size(); row ++) {
        if (givens[row][col]) {
            givens_in_col.insert(givens[row][col]);
        }
    }
    return givens_in_col;
}

bool is_in_bounding_house(int row, int col, Point tl, Point br) {
    if (row >= tl.row && col >= tl.col) {
    if (row <= br.row && col <= br.col) {
        return true;
    }
    }
    return false;
}

int get_house_for_cell(int row, int col, std::map<int, std::pair<Point, Point>> houses) {
    for (auto it = houses.begin(); it != houses.end(); it++) {
    if (is_in_bounding_house(row, col, it->second.first, it->second.second)) {
        return it->first;
    }
    }
    return -1;
}

std::set<Point> points_in_house(int row, int col, std::map<int, std::pair<Point, Point>> houses) {
    /*for (auto& x : houses) {
        std::cout << "house " << x.first << " co-ordinates, top left: " << x.second.first << " bottom right: " << x.second.second << std::endl;
    }*/
    int house = get_house_for_cell(row, col, houses);
    //std::cout << "house for cell " << row << " " << col << " is " << house << std::endl;
    auto& [top_left, bottom_right] = houses[house];
    //std::cout << "top left " << top_left << " br " << bottom_right << std::endl;

    std::set<Point> points;
    for (int i = top_left.row; i <= bottom_right.row; i++) {
    for (int j = top_left.col; j <= bottom_right.col; j++) {
        points.insert(Point{.row = i, .col = j});
    }}

    return points;
}

std::set<int> givens_in_house(int row, int col, std::vector<std::vector<int>> givens, std::map<int, std::pair<Point, Point>> houses) {
    auto points_house = points_in_house(row, col, houses);

    std::set<int> givens_in_house;
    for (uint32_t i = 0; i < points_house.size(); i++) {
        for (auto& p : points_house) {
            if (givens[p.row][p.col]) {
                givens_in_house.insert(givens[p.row][p.col]);
            }
        }
    }
    return givens_in_house;
}

CellCandidate CreateCellCandidate(int row, int col, std::set<int> possibles, std::map<int, std::pair<Point, Point>> houses) {
    int house = get_house_for_cell(row, col, houses);
    return {
        .row = row,
        .col = col,
        .house = house,
        .candidates = std::move(possibles)
    };
}

std::map<Point, CellCandidate>
ComputeInitialCandidates(int rank,
                         std::vector<std::vector<int>> givens,
                         std::map<int, std::pair<Point, Point>> houses) {
    int rows = rank * rank;
    int cols = rank * rank;
    std::map<Point, CellCandidate> candidates;
    std::set<int> possible_assignments = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (givens[i][j] != 0) {
                // If a value is given for a cell, we should not have to consider the candidates for that cell.
                // candidates.insert({Point {.row = i, .col = j}, CreateCellCandidate(i, j, std::set<int>{givens[i][j]}, houses)});
                continue;
            }
            std::set<int> all_givens = givens_in_row(i, givens);
            std::set<int> givens_col = givens_in_col(i, givens);
            all_givens.insert(givens_col.begin(), givens_col.end());
            std::set<int> givens_house = givens_in_house(i, j, givens, houses);
            all_givens.insert(givens_house.begin(), givens_house.end());

            std::set<int> possibles;
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
    int get_sudoku_id()
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
    int get_rank()
    {
        return rank;
    }
    std::vector<std::vector<int>> get_givens()
    {
        return givens;
    }
    std::map<int, std::pair<Point, Point>> get_houses()
    {
        return houses;
    }

    // Delete default ctor.
    ParsedSudoku() = delete;

    static std::optional<std::unique_ptr<ParsedSudoku>> ParsedSudokuFactory(fs::path problem_file_path);

private:
    ParsedSudoku(
        int sudoku_id, std::string problem_file_path,
        std::string compressed_repr, int rank,
        std::vector<std::vector<int>> givens,
        std::map<int, std::pair<Point, Point>> houses)
    {
        this->sudoku_id = sudoku_id;
        this->problem_file_path = problem_file_path;
        this->compressed_repr = compressed_repr;
        this->rank = rank;
        this->givens = std::move(givens);
        this->houses = std::move(houses);
    }

    int sudoku_id;
    std::string problem_file_path;
    std::string compressed_repr;
    int rank;
    std::map<int, std::pair<Point, Point>> houses;
    std::vector<std::vector<int>> givens;
};

// Factory to generate a fully initialized ParsedSudoku.
std::optional<std::unique_ptr<ParsedSudoku>> ParsedSudoku::ParsedSudokuFactory(fs::path problem_file_path)
{
    std::string file_contents;
    std::ifstream sudoku_fs(problem_file_path);
    int sudoku_id;
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
    int rank = valuerowcolumn.size() >= 17 ? 3 : 2;

    std::vector<std::vector<int>> givens = {};
    for (int i = 0; i < rank * rank; i++)
    {
        givens.push_back({});
        for (int j = 0; j < rank * rank; j++)
        {
            givens[i].push_back(0);
        }
    }

    for (std::string vrc : valuerowcolumn)
    {
        int val = vrc[0] - '0';
        int row = vrc[1] - '0';
        int col = vrc[2] - '0';
        if (row >= (int)givens.size() || col >= (int)givens[0].size()) {
            std::cout << "Invalid row and column csv " << row << col << std::endl;
            continue;
        }
        givens[row][col] = val;
    }

    std::map<int, std::pair<Point, Point>> houses;
    int house_nr = 0;
    while (house_nr < rank * rank) {
        int row = house_nr;
        int col = 0;
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
        SudokuSolution solution;
        SolveSudokuAndRecurse(solver_state.get(), solution);
    }

    Solver(std::unique_ptr<ParsedSudoku> sudoku, fs::path output_folder)
    {
        this->sudoku = std::move(sudoku);
        this->solution_file_path = output_folder;
        BuildInitialSudokuState();
    }

    Solver() = delete;

private:
    class SudokuStepState {
        public:
            std::vector<std::vector<int>>& get_assignments() {
                return assignments;
            }

            void Add(uint32_t value) {
               CellCandidate* cell = pq.top();
               std::cout<< "Assigning " << value << " to " << "[" << cell->row << ", " << cell->col << "]" << std::endl;
               assignments[cell->row][cell->col] = value;
               UpdateAdjacentHouses(value);
               pq.pop();
            }

            void Subtract(CellCandidate* cell, uint32_t value) {
                UndoUpdateToAdjacentHouses(value);
                // The priority queue should be updated after the adjacent cells have been updated.
                pq.push(cell);
                assignments[cell->row][cell->col] = 0;
            }

            CellCandidate* OptimalNextCell() {
                if (pq.empty()) {
                    return nullptr;
                }
                return pq.top();
            }

            std::vector<std::vector<int>> get_cell_assignments() {
                return assignments;
            }

            SudokuStepState(std::map<Point, CellCandidate> candidates,
                            ParsedSudoku* sudoku) {
                pq_map = std::move(candidates);
                for (auto it = pq_map.begin(); it != pq_map.end(); it++) {
                    CellCandidate* cell = &it->second;
                    pq.push(cell);
                }
                std::cout << "Printing Givens " << " size is " << sudoku->get_givens().size() << " columns is " <<  sudoku->get_givens()[0].size() << std::endl;
                for (uint32_t i = 0; i < sudoku->get_givens().size(); i++) {
                    this->assignments.push_back({});
                    for (auto j = 0; j < sudoku->get_givens()[0].size(); j++) {
                        this->assignments[i].push_back(sudoku->get_givens()[i][j]);
                    }
                }
                sudoku_ptr = sudoku;
            }

            SudokuStepState() = delete;

        private:
            std::set<Point> get_adjacent_points(CellCandidate& cell);
            void UndoUpdateToAdjacentHouses(int value_to_be_readded);
            void UpdateAdjacentHouses(uint32_t value);

            // Authoritative owner of all CellCandidate(s).
            std::map<Point, CellCandidate> pq_map;
            std::priority_queue<CellCandidate*> pq;
            // Current assignments to the coordinates of the sudoku.
            std::vector<std::vector<int>> assignments;
            // For any step, keeps track of the cells updated by a new cell assignment.
            std::set<CellCandidate*> cells_updated;
            ParsedSudoku* sudoku_ptr;
    };

    // TODO: Add and verify difficulty metrics.
    // This should be a combination of:
    // 1) Number of assignments already done
    // 2) The number of candidates in the optimal next step.
    // In the future we can also consider other cells which are close to the most optimal,
    // let us say, top 10 cells or so and their average difficulty. This would slow down the
    // sudoku program by a factor though.
    void SolveSudokuAndRecurse(SudokuStepState* state, SudokuSolution& solution);

    void BuildInitialSudokuState() {
        auto cell_candidates = ComputeInitialCandidates(sudoku->get_rank(), sudoku->get_givens(), sudoku->get_houses());
        solver_state = std::make_unique<SudokuStepState>(std::move(cell_candidates), sudoku.get());
        printf("built initial sudoku state\n");
    }

    void WriteSudokuSolutionToFile() {}

    std::unique_ptr<ParsedSudoku> sudoku;
    fs::path solution_file_path;
    std::unique_ptr<SudokuStepState> solver_state;
    std::vector<SudokuSolution> found_solutions;
};

std::set<Point> Solver::SudokuStepState::get_adjacent_points(CellCandidate& cell) {
    std::set<Point> adj_cols, adj_rows, adj_houses;
    int N = (int) assignments.size();

    for (int col = 0; col < N; col ++) {
        adj_cols.insert(Point{.row = cell.row, .col = col});
    }
    for (int row = 0; row < N; row ++) {
        adj_rows.insert(Point{.row = row, .col = cell.col});
    }
    adj_houses = points_in_house(cell.row, cell.col, sudoku_ptr->get_houses());

    std::set<Point> adj_points;
    adj_points.insert(adj_rows.begin(), adj_rows.end());
    adj_points.insert(adj_cols.begin(), adj_cols.end());
    adj_points.insert(adj_houses.begin(), adj_houses.end());
    return adj_points;
}

void Solver::SudokuStepState::UndoUpdateToAdjacentHouses(int value_to_be_readded) {
    for (CellCandidate* cell : cells_updated) {
        cell->candidates.insert(value_to_be_readded);
    }
    cells_updated.clear();
}

void Solver::SudokuStepState::UpdateAdjacentHouses(uint32_t value) {
    CellCandidate& cell = *pq.top();
    //std::cout << "getting adjacent points for " << cell.row << " " << cell.col << std::endl;
    std::set<Point> adj_points = get_adjacent_points(cell);
    //std::cout << "got adjacent points";
    for (auto& update_coord : adj_points) {
        auto& update_cell = pq_map[update_coord];
        // Assigned value was a candidate for this cell.
        if (update_cell.candidates.contains(value)) {
            update_cell.candidates.erase(value);
            cells_updated.insert(&update_cell);
        }
    }
}

void Solver::SolveSudokuAndRecurse(SudokuStepState* state, SudokuSolution& solution) {
    CellCandidate* optimal_next_cell = state->OptimalNextCell();
    if (optimal_next_cell == nullptr) {
        solution.assignments = state->get_cell_assignments();
        for (uint32_t i = 0; i < solution.assignments.size(); i++) {
            for (uint32_t j = 0; j < solution.assignments[0].size(); j++) {
                printf("%d\t", solution.assignments[i][j]);
            }
            printf("\n");
        }
        found_solutions.push_back(solution);
        return;
    }

    if (optimal_next_cell->candidates.size() == 0) {
        // We made an incorrect assignment.
        std::cout<< "No options for cell " << optimal_next_cell->row << " " << optimal_next_cell->col << std::endl;
        return;
    }

    print_sudoku_assignments(state->get_assignments());

    // Pick a candidate coloring for the cell.
    std::set<int> candidates_to_pick = std::set<int>(optimal_next_cell->candidates);
    for (int picked_candidate : candidates_to_pick) {
        // std::cout << "Picked candidate " << picked_candidate << std::endl;
        state->Add(picked_candidate);
        CellCandidate* copy = optimal_next_cell;
        SolveSudokuAndRecurse(state, solution);
        // Undo previous update.
        state->Subtract(copy, picked_candidate);
    }
}

class SolveSudokusInFolder
{
public:
    void AnalyseSudokus(std::string input_folder, std::string output_folder)
    {
        std::vector<fs::path> sudoku_pathnames = ReadDirectory(input_folder);
        std::vector<std::unique_ptr<ParsedSudoku>> parsed_sudokus = ParseSudokuProblems(sudoku_pathnames);
        std::cout << "starting threads parsed_sudokus.size() = " << parsed_sudokus.size() << std::endl;

        // For each Sudoku, we create a thread to solve and publish the result of the sudoku
        // as a CSV file into the output folder.
        std::vector<std::thread> sudoku_threads;
        std::vector<Solver> solvers;
//        for (unsigned int i = 0; i < parsed_sudokus.size(); i++)
//        {
        solvers.push_back(Solver(std::move(parsed_sudokus[0]), fs::path(output_folder)));
            //std::thread new_thread = std::thread(&Solver::SolveSudoku, &solvers[i]);
            //new_thread.join();
            // sudoku_threads.push_back(std::move(new_thread));
//        }
        solvers[0].SolveSudoku();
        /*            std::this_thread::sleep_for(std::chrono::seconds(5));
                    for (auto& thread : sudoku_threads) {
                        while (!thread.joinable()) {
                            continue;
                        }
                        thread.join();
                    }*/

        std::cout << "threads joined" << std::endl;
    }

private:
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
};

int main()
{
    SolveSudokusInFolder solver;
    solver.AnalyseSudokus("parsed_sudokus", "sudoku_solutions");
}
