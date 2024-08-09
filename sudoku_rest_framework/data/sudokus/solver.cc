/*
 * C++ Program to solve a sudoku and analyse:
 *   1. Whether the solution is unique.
 *   2. How difficult is the sudoku.
 * 
 * Before changes to make priority queue work correctly, time taken was 673204[µs].
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
    uint32_t sudoku_id;
    uint32_t num_of_easy_cells;
    uint32_t candidates_per_easy_cell;
    uint32_t average_number_of_candidates_per_cell;
} StepDifficulty;

typedef struct Difficulty
{
    uint32_t sudoku_id;
    std::vector<StepDifficulty> steps_to_success;
    std::vector<StepDifficulty> steps_to_failure;
} Difficulty;

typedef struct SudokuSolution
{
    uint32_t sudoku_id;
    std::vector<std::vector<uint32_t>> assignments;
    Difficulty analysed_difficulty;
} SudokuSolution;

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
} CellCandidatePtr;

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

void print_sudoku_assignments(std::vector<std::vector<uint32_t>>& assignments) {
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

uint32_t get_house_for_cell(uint32_t row, uint32_t col, std::map<uint32_t, std::pair<Point, Point>> houses) {
    for (auto it = houses.begin(); it != houses.end(); it++) {
    if (is_in_bounding_house(row, col, it->second.first, it->second.second)) {
        return it->first;
    }
    }
    return -1;
}

std::set<Point> points_in_house(uint32_t row, uint32_t col, std::map<uint32_t, std::pair<Point, Point>> houses) {
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
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        SudokuSolution solution;
        SolveSudokuAndRecurse(solver_state.get(), solution);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time taken = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

        if (found_solutions.size() > 0) {
            printf("Found solution: \n");
            for (auto& solution : found_solutions) {
                for (uint32_t i = 0; i < solution.assignments.size(); i++) {
                    for (uint32_t j = 0; j < solution.assignments[0].size(); j++) {
                        printf("%d\t", solution.assignments[i][j]);
                    }
                    printf("\n");
                }
                printf("\n\n");
            } 
        }

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
            std::vector<std::vector<uint32_t>>& get_assignments() {
                return assignments;
            }

            void Add(CellCandidate* cell, uint32_t value) {
               assignments[cell->row][cell->col] = value;
               UpdateAdjacentHouses(cell, value);
            }

            void Subtract(CellCandidate* cell, uint32_t value) {
                UndoUpdateToAdjacentHouses(cell, value);
                assignments[cell->row][cell->col] = 0;
            }

            CellCandidate* PopAndGetOptimalNextCell() {
                if (pq.empty()) {
                    return nullptr;
                }
                CellCandidate* top = pq.top().ptr;
                pq.pop();
                return top;
            }

            void PushToPQ(CellCandidate* cell) {
                pq.push(CellCandidatePtr{.ptr = cell});
            }

            std::vector<std::vector<uint32_t>> get_cell_assignments() {
                return assignments;
            }

            SudokuStepState(std::map<Point, CellCandidate> candidates,
                            ParsedSudoku* sudoku) {
                pq_map = std::move(candidates);
                for (auto it = pq_map.begin(); it != pq_map.end(); it++) {
                    CellCandidate* cell = &it->second;
                    pq.push(CellCandidatePtr{.ptr = cell});
                }
                // std::cout << "Printing Givens " << " size is " << sudoku->get_givens().size() << " columns is " <<  sudoku->get_givens()[0].size() << std::endl;
                for (uint32_t i = 0; i < sudoku->get_givens().size(); i++) {
                    this->assignments.push_back({});
                    for (uint32_t j = 0; j < sudoku->get_givens()[0].size(); j++) {
                        this->assignments[i].push_back(sudoku->get_givens()[i][j]);
                    }
                }
                sudoku_ptr = sudoku;
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
                    printf("[%u, %u] : size = %lu\n", top.ptr->row, top.ptr->col, top.ptr->candidates.size());
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
                std::cout << std::endl;
            }
        private:
            std::set<Point> get_adjacent_points(CellCandidate& cell);
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
    }

    void WriteSudokuSolutionToFile() {}

    std::unique_ptr<ParsedSudoku> sudoku;
    fs::path solution_file_path;
    std::unique_ptr<SudokuStepState> solver_state;
    std::vector<SudokuSolution> found_solutions;
};

void Solver::SudokuStepState::UndoUpdateToAdjacentHouses(CellCandidate* cell, uint32_t value_to_be_readded) {
    Point cell_coord = Point({.row = cell->row, .col = cell->col});
    //std::cout << "Reverting everything for " << cell_coord << std::endl;
    auto it = revert_updates.find(cell_coord);
    if (it != revert_updates.end()) {
        auto& [_, re_add_set] = *it;
        // What were the updates made by cell at cell_coord, and to which points.
        for (auto& [x_value, x_coord] : re_add_set) {
            //std::cout<< "CellCandidate " << cell_coord << " was reverted " << value_to_be_readded << " and adjacent cell " << x_coord << " was updated." << std::endl;
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
    std::set<Point> adj_points = get_adjacent_points(*cell);
    //auto it = revert_updates.end();

    for (auto& update_coord : adj_points) {
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

void Solver::SolveSudokuAndRecurse(SudokuStepState* state, SudokuSolution& solution) {
    CellCandidate* optimal_next_cell = state->PopAndGetOptimalNextCell();
    //printf("Next cell from heap is [%u, %u]\n", optimal_next_cell->row, optimal_next_cell->col);
    if (optimal_next_cell == nullptr) {
        solution.assignments = state->get_cell_assignments();
        return;
    }

    if (optimal_next_cell->candidates.size() == 0) {
        // We made an incorrect assignment.
        //std::cout<< "No options for cell " << optimal_next_cell->row << " " << optimal_next_cell->col << std::endl;
        return;
    }

    std::set<uint32_t> candidates_to_pick = std::set<uint32_t>(optimal_next_cell->candidates);
    /*std::cout<< "\nCandidates for [" << optimal_next_cell->row << ", " << optimal_next_cell->col << "]" << " are ";
    for (uint32_t c : candidates_to_pick) {
        std::cout << c << ",";
    }
    std::cout << std::endl;*/

    if (state->get_cell_assignments()[optimal_next_cell->row][optimal_next_cell->col] != 0) {
        std::cout<< "EEEEEEEEERRRRRRRRRRRRROOOOOOOOOOORRRRRRRRRRRRRRR this cell is assigned, already" << std::endl;
    }

    //print_sudoku_assignments(state->get_assignments());

    // Pick a candidate coloring for the cell.
    for (uint32_t picked_candidate : candidates_to_pick) {
        //std::cout<< "Assigning " << picked_candidate << " to " << "[" << optimal_next_cell->row << ", " << optimal_next_cell->col << "]" << std::endl;
        //printf("[%u, %u] Before adding\n", copy->row, copy->col);
        //state->print_priority_queue(0);
        state->Add(optimal_next_cell, picked_candidate);
        //printf("[%u, %u]Entering recursion\n", optimal_next_cell->row, optimal_next_cell->col);
        //state->print_priority_queue(5);
        SolveSudokuAndRecurse(state, solution);
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

std::set<Point> Solver::SudokuStepState::get_adjacent_points(CellCandidate& cell) {
    std::set<Point> adj_cols, adj_rows, adj_houses;
    uint32_t N =  assignments.size();

    for (uint32_t col = 0; col < N; col ++) {
        adj_cols.insert(Point{.row = cell.row, .col = col});
    }
    for (uint32_t row = 0; row < N; row ++) {
        adj_rows.insert(Point{.row = row, .col = cell.col});
    }
    adj_houses = points_in_house(cell.row, cell.col, sudoku_ptr->get_houses());

    std::set<Point> adj_points;
    adj_points.insert(adj_rows.begin(), adj_rows.end());
    adj_points.insert(adj_cols.begin(), adj_cols.end());
    adj_points.insert(adj_houses.begin(), adj_houses.end());
    return adj_points;
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
//        for (unsigned uint32_t i = 0; i < parsed_sudokus.size(); i++)
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
