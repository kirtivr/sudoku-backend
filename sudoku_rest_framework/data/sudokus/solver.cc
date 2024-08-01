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
#include <optional>
#include <iostream>
#include <fstream>
#include <thread>
#include <filesystem>
#include <queue>
#include <algorithm>
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
    int box; // 0 - 3/8.
    std::set<int> candidates;
} CellCandidate;

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
    std::map<std::pair<int, int>, CellCandidate> get_candidates()
    {
        return candidates;
    }

    std::set<int> givens_in_row(int row, std::vector<std::vector<int>> givens) {
        std::set<int> givens_in_row {};
        for (int col = 0; col < givens[0].size(); col ++) {
            givens_in_row.insert(givens[row][col]);
        }
        return givens_in_row;
    }

    std::set<int> givens_in_col(int col, std::vector<std::vector<int>> givens) {
        std::set<int> givens_in_col {};
        for (int row = 0; row < givens.size(); row ++) {
            givens_in_col.insert(givens[row][col]);
        }
        return givens_in_col;
    }

    std::set<int> givens_in_box(int row, int col, std::vector<std::vector<int>> givens) {
        // 4 Quadrants.
        // If rank is odd:
        //  [0, 0] to [rank * rank // 2, rank * rank // 2]
        //  [rank * rank // 2 + 1, 0] to [rank * rank - 1, rank * rank // 2]
        //  [0, rank * rank // 2 + 1] to [rank * rank // 2, rank * rank - 1]
        //  [rank * rank // 2 + 1, [rank * rank // 2 + 1] to [rank * rank - 1, rank * rank - 1]

        int rank = math.sqrt(givens.size());
    }

    std::map<std::pair<int, int> CellCandidate> ComputeInitialCandidates(int rank, std::vector<std::vector<int>> givens) {
        int rows = rank * rank;
        int cols = rank * rank;
        std::map<std::pair<int, int> candidates;
        std::set<int> possible_assignments = {1, 2, 3, 4, 5, 6, 7, 8, 9};

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (givens[i][j] != 0) {
                    candidates.insert(std::make_pair<int, int>(i, j), std::set<int>{givens[i][j]});
                    continue;
                }
                std::set<int> all_givens = givens_in_row(i, givens);
                std::set<int> givens_col = givens_in_col(i, givens);
                all_givens.insert(givens_col.begin(), givens_col.end());
                std::set<int> givens_box = givens_in_box(i, j, givens);
                all_givens.insert(givens_box.begin(), givens_box.end());

                std::set<int> possibles;
                std::set_difference(possible_assignments.begin(), possible_assignments.end(),
                                    all_givens.begin(), all_givens.end(), std::inserter(possibles, possibles.end()));

                candidates.insert(std::make_pair<int, int>(i, j), possibles);
            }
        }

        return candidates;
    }

    // Factory to generate a fully initialized ParsedSudoku.
    static std::optional<ParsedSudoku> Parse(fs::path problem_file_path)
    {
        std::string givens;
        std::ifstream sudoku_fs(problem_file_path);
        int sudoku_id;
        if (sudoku_fs.is_open())
        {
            std::getline(sudoku_fs, givens);
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
        std::vector<std::string> rowcolumnvalue = splitstring(givens, ",");
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
            if (row >= givens.size() || col >= givens[0].size()) {
                std::cout << "Invalid row and column csv " << row << col << std::endl;
                continue;
            }
            givens[row][col] = val;
        }

        auto cell_candidates = ComputeInitialCandidates(rank, givens);
        return std::optional<ParsedSudoku>(
            ParsedSudoku(sudoku_id, problem_file_path, givens, rank,
            std::move(givens), std::move(cell_candidates)));
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

    // Delete default ctor.
    ParsedSudoku() = delete;

private:
    ParsedSudoku(
        int sudoku_id, std::string problem_file_path,
        std::string compressed_repr, int rank,
        std::vector<std::vector<int>> givens,
        std::map<std::pair<int, int>, CellCandidate> candidates)
    {
        this->sudoku_id = sudoku_id;
        this->problem_file_path = problem_file_path;
        this->compressed_repr = compressed_repr;
        this->givens = std::move(givens);
        this->candidates = std::move(candidates);
    }

    int sudoku_id;
    std::string problem_file_path;
    std::string compressed_repr;
    int rank;
    std::map<std::pair<int, int>, CellCandidate> candidates;
    std::vector<std::vector<int>> givens;
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
    void SolveSudoku() {

    }

    Solver(ParsedSudoku *sudoku, SudokuSolution *solution, fs::path output_folder)
    {
        sudoku = sudoku;
        solution_file_path = output_folder;
        solution = solution;
    }

    Solver() = delete;

private:
    void WriteSudokuSolutionToFile()
    {
    }

    ParsedSudoku *sudoku = nullptr;
    fs::path solution_file_path;
    SudokuSolution *solution;

public:
    int what = 0;
};

class SolveSudokusInFolder
{
public:
    void AnalyseSudokus(std::string input_folder, std::string output_folder)
    {
        std::vector<fs::path> sudoku_pathnames = GetSudokuPathnames(input_folder);
        std::vector<ParsedSudoku> parsed_sudokus = ParseSudokuProblems(sudoku_pathnames);
        std::vector<SudokuSolution> solutions(parsed_sudokus.size());
        std::cout << "starting threads parsed_sudokus.size() = " << parsed_sudokus.size() << std::endl;

        // For each Sudoku, we create a thread to solve and publish the result of the sudoku
        // as a CSV file into the output folder.
        std::vector<std::thread> sudoku_threads;
        std::vector<Solver> solvers;
        for (unsigned int i = 0; i < parsed_sudokus.size(); i++)
        {
            solvers.emplace_back(&parsed_sudokus[i], &solutions[i], output_folder);
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
        for (unsigned int i = 0; i < parsed_sudokus.size(); i++)
        {
            std::cout << solvers[i].what << std::endl;
        }
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

    std::vector<ParsedSudoku> ParseSudokuProblems(std::vector<fs::path> &sudoku_problem_files)
    {
        std::vector<ParsedSudoku> parsed_sudokus = {};
        for (auto &file : sudoku_problem_files)
        {
            // std::cout << "passing file " << file << std::endl;
            std::optional<ParsedSudoku> maybe_parsed_file = ParsedSudoku::Parse(file);
            if (maybe_parsed_file != std::nullopt)
            {
                parsed_sudokus.push_back(*maybe_parsed_file);
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