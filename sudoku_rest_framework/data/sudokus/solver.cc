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
namespace fs = std::filesystem;

std::optional<std::vector<fs::path>> ReadDirectory(std::string directory_path)
{
    DIR *ptr = opendir(directory_path.c_str());
    if (ptr == NULL)
    {
        printf("Could not open given directory: %s", directory_path.c_str());
        return std::nullopt;
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
    std::vector<std::vector<int>> solution;
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
} Point;

std::set<int> givens_in_row(int row, std::vector<std::vector<int>> givens) {
    std::set<int> givens_in_row {};
    for (uint32_t col = 0; col < givens[0].size(); col ++) {
        givens_in_row.insert(givens[row][col]);
    }
    return givens_in_row;
}

std::set<int> givens_in_col(int col, std::vector<std::vector<int>> givens) {
    std::set<int> givens_in_col {};
    for (uint32_t row = 0; row < givens.size(); row ++) {
        givens_in_col.insert(givens[row][col]);
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

std::set<int> givens_in_house(int row, int col, std::vector<std::vector<int>> givens, std::map<int, std::pair<Point, Point>> houses) {
    // 4 Quadrants.
    // If rank is odd:
    //  [0, 0] to [rank * rank // 2, rank * rank // 2]
    //  [rank * rank // 2 + 1, 0] to [rank * rank - 1, rank * rank // 2]
    //  [0, rank * rank // 2 + 1] to [rank * rank // 2, rank * rank - 1]
    //  [rank * rank // 2 + 1, [rank * rank // 2 + 1] to [rank * rank - 1, rank * rank - 1]
    int house = get_house_for_cell(row, col, houses);
    auto [top_left, bottom_right] = houses[house];

    std::set<int> givens_in_house {};
    for (int row = top_left.row; row <= bottom_right.row; row++) {
    for (int col = top_left.col; row <= bottom_right.col; col++) {
        givens_in_house.insert(givens[row][col]);
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

std::map<Point, CellCandidate> ComputeInitialCandidates(int rank, std::vector<std::vector<int>> givens, std::map<int, std::pair<Point, Point>> houses) {
    int rows = rank * rank;
    int cols = rank * rank;
    std::map<Point, CellCandidate> candidates;
    std::set<int> possible_assignments = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (givens[i][j] != 0) {
                candidates.insert({Point {.row = i, .col = j}, CreateCellCandidate(i, j, std::set<int>{givens[i][j]}, houses)});
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
    std::map<Point, CellCandidate> get_candidates()
    {
        return candidates;
    }

    // Delete default ctor.
    ParsedSudoku() = delete;

    static std::optional<std::unique_ptr<ParsedSudoku>> ParsedSudokuFactory(fs::path problem_file_path);

private:
    ParsedSudoku(
        int sudoku_id, std::string problem_file_path,
        std::string compressed_repr, int rank,
        std::vector<std::vector<int>> givens,
        std::map<Point, CellCandidate> candidates,
        std::map<int, std::pair<Point, Point>> houses)
    {
        this->sudoku_id = sudoku_id;
        this->problem_file_path = problem_file_path;
        this->compressed_repr = compressed_repr;
        this->givens = std::move(givens);
        this->houses = std::move(houses);
        this->candidates = std::move(candidates);
    }

    int sudoku_id;
    std::string problem_file_path;
    std::string compressed_repr;
    int rank;
    std::map<Point, CellCandidate> candidates;
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
    // rowcolumnvalue std::strings.
    // First find the Rank of the sudoku by counting the
    // comma separated values.
    // Then populate the candidate std::vector.
    std::vector<std::string> rowcolumnvalue = splitstring(file_contents, ",");
    int rank = rowcolumnvalue.size() >= 17 ? 3 : 2;

    std::vector<std::vector<int>> givens = {};
    for (int i = 0; i < rank * rank; i++)
    {
        givens.push_back({});
        for (int j = 0; j < rank * rank; j++)
        {
            givens[i].push_back(0);
        }
    }

    std::cout << "rank = " << rank;

    for (std::string rcv : rowcolumnvalue)
    {
        int row = rcv[0] - '0';
        int col = rcv[1] - '0';
        int val = rcv[2] - '0';
        if (row >= (int)givens.size() || col >= (int)givens[0].size()) {
            std::cout << "Invalid row and column csv " << row << col << std::endl;
            continue;
        }
        givens[row][col] = val;
    }

    int mid = rank % 2 == 0 ? rank * rank / 2 - 1 : rank * rank / 2;
    int end = rank * rank - 1;
    std::map<int, std::pair<Point, Point>> houses = {
        {0, std::make_pair<Point, Point>({.row = 0, .col = 0}, {.row = mid, .col = mid})},
        {1, std::make_pair<Point, Point>({.row = 0, .col = mid + 1}, {.row = mid, .col = end})},
        {2, std::make_pair<Point, Point>({.row = mid + 1, .col = 0}, {.row = end, .col = mid})},
        {3, std::make_pair<Point, Point>({.row = mid + 1, .col = mid + 1}, {.row = end, .col = end})}
    };
    auto cell_candidates = ComputeInitialCandidates(rank, givens, houses);
    return std::make_unique<ParsedSudoku>(
        ParsedSudoku(sudoku_id, problem_file_path, file_contents, rank, std::move(givens),
        std::move(cell_candidates), std::move(houses)));
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
        std::copy(sudoku->givens.begin(), sudoku->givens.end(), std::back_inserter(solution.solution));  
        SolveSudokuAndRecurse(step_zero_pq, solution);
    }

    Solver(std::unique_ptr<ParsedSudoku> sudoku, std::unique_ptr<SudokuSolution> solution, fs::path output_folder)
    {
        this->sudoku = std::move(sudoku);
        this->solution = std::move(solution);
        this->solution_file_path = output_folder;
        BuildInitialPriorityQueue();
    }

    Solver() = delete;

private:
    void UpdateCandidateCells() {
        
    }

    uint32_t generate_random_number_in_range(uint32_t low, uint32_t high, int picked_already[]) {
        std::random_device rd; // obtain a random number from hardware
        std::mt19937 gen(rd()); // seed the generator
        std::uniform_int_distribution<> distr(0, max_index); // define the range
        int generated = distr(gen);
        while (picked_already[generated] == 0) {
            generated = distr(gen);
        }
        return generated;
    }

    void SolveSudokuAndRecurse(std::priority_queue pq, SudokuSolution solution) {
        if (pq.empty()) {
            solutions.push_back(solution);
            return;
        }

        CandidateCell& optimal_next_cell = pq.top();

        // Pick a candidate coloring for the cell.
        int num_candidates = optimal_next_cell.candidates.size();
        int picked_already[num_candidates] = {0};
        int num_picked = 0;

        while (num_picked < num_candidates) {
            int picked_candidate =
                optimal_next_cell[generate_random_number_in_range(0, optimal_next_cell.candidates.size() - 1, picked_already)];
            picked_already[picked_candidate] = 1;
            pq.pop();
            optimal_next_cell.candidates = {picked_candidate};
            solution.solution[optimal_next_cell.row][optimal_next_cell.col] = picked_candidate;

            // Update the candidate cells in the priority queue.
            UpdateCandidateCells();
            SolveSudokuAndRecurse(pq, solution);

            // Undo previous update.
            UpdateCandidateCells();
            solution.solution[optimal_next_cell.row][optimal_next_cell.col] = 0;
        }
    }

    void BuildInitialPriorityQueue() {
        for (auto it = sudoku->get_candidates().begin();
             it != sudoku->get_candidates().end(); it++) {
                // Does this make a copy? .
                step_zero_pq.push(it->second);
             }
    }

    void WriteSudokuSolutionToFile() {}

    std::unique_ptr<ParsedSudoku> sudoku;
    fs::path solution_file_path;
    std::vector<SudokuSolution> found_solutions;
    std::priority_queue<CellCandidate> step_zero_pq;
};

class SolveSudokusInFolder
{
public:
    void AnalyseSudokus(std::string input_folder, std::string output_folder)
    {
        std::vector<fs::path> sudoku_pathnames = GetSudokuPathnames(input_folder);
        std::vector<std::unique_ptr<ParsedSudoku>> parsed_sudokus = ParseSudokuProblems(sudoku_pathnames);
        std::cout << "starting threads parsed_sudokus.size() = " << parsed_sudokus.size() << std::endl;

        // For each Sudoku, we create a thread to solve and publish the result of the sudoku
        // as a CSV file into the output folder.
        std::vector<std::thread> sudoku_threads;
        std::vector<Solver> solvers;
        for (unsigned int i = 0; i < parsed_sudokus.size(); i++)
        {
            solvers.emplace_back(std::move(parsed_sudokus[i]), parsed_sudokus[i]->givens, output_folder);
            std::thread new_thread = std::thread(&Solver::SolveSudoku, &solvers[i]);
            new_thread.join();
            // sudoku_threads.push_back(std::move(new_thread));
        }
        /*            std::this_thread::sleep_for(std::chrono::seconds(5));
                    for (auto& thread : sudoku_threads) {
                        while (!thread.joinable()) {
                            continue;
                        }
                        thread.join();
                    }*/

        std::cout << "threads joined parsed_sudokus.size() = " << parsed_sudokus.size() << std::endl;
    }

private:
    std::vector<fs::path> GetSudokuPathnames(std::string folder_path)
    {
        // std::string docPath = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
        std::optional<std::vector<fs::path>> files = ReadDirectory(folder_path);
        if (files != std::nullopt)
        {
            return *files;
        }
        return {};
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
};

int main()
{
    SolveSudokusInFolder solver;
    solver.AnalyseSudokus("parsed_sudokus", "sudoku_solutions");
}
