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

std::optional<std::vector<std::string>> ReadDirectory(std::string directory_path) {
    DIR* ptr = opendir(directory_path.c_str());

    if (ptr == NULL) {
        printf("Could not open given directory: %s", directory_path.c_str());
        return std::nullopt;
    }

    struct dirent *entry;
    std::vector<std::string> directory_contents;
    while ((entry = readdir(ptr)) != NULL)  {
        if (strncmp(entry->d_name, ".", 256) == 0 || strncmp(entry->d_name, "..", 256) == 0) {
            continue;
        }
        directory_contents.push_back(entry->d_name);
    }

    closedir(ptr);
    return directory_contents;
}

// This class is thread-compatible.
// It is only ever used by one thread.
class ParsedSudoku
{
public:
    int get_sudoku_id() {
        return sudoku_id;
    }
    std::string get_problem_file_path() {
        return problem_file_path;
    }
    std::string get_compressed_repr() {
        return compressed_repr;
    }
    int get_rank() {
        return rank;
    }
    std::vector<std::vector<int>> get_candidates() {
        return candidates;
    }
    
    static std::optional<ParsedSudoku> Parse(std::string problem_file_path) {
        std::string givens;
        std::ifstream sudoku_fs(problem_file_path);
        int sudoku_id;
        if (sudoku_fs.is_open()) {
            std::getline(sudoku_fs, givens);
            sudoku_id = std::stoi(problem_file_path.substr(0, givens.length() - 4).c_str());
        } else {
            return std::nullopt;
        }

        // Givens is one std::string of comma separated:
        // rowcolumnvalue std::strings.
        // First find the Rank of the sudoku by counting the
        // comma separated values.
        // Then populate the candidate std::vector.
        std::vector<std::string> rowcolumnvalue = splitstring(givens, ",");
        int rank = rowcolumnvalue.size() >= 17 ? 3 : 2;

        std::vector<std::vector<int>> candidates = {};
        for (int i = 0; i < rank * rank; i++) {
            candidates.push_back({});
            for (int j = 0; j < rank * rank; j++) {
                candidates[i].push_back(0);
            }
        }

        for (std::string rcv :rowcolumnvalue) {
            int row = int(rcv[0]);
            int col = int(rcv[1]);
            int val = int(rcv[2]);
            candidates[row][col] = val;
        }

        return std::optional<ParsedSudoku>(ParsedSudoku(sudoku_id, problem_file_path, givens, rank, candidates));
    }

    static std::vector<std::string> splitstring(std::string given, std::string token) {
        std::size_t pos_start = 0, pos_end, delim_len = token.length();
        std::string chunk;
        std::vector<std::string> out;

        while ((pos_end = given.find(token, pos_start)) != std::string::npos) {
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
        std::vector<std::vector<int>> candidates) {
            this->sudoku_id = sudoku_id;
            this->problem_file_path = problem_file_path;
            this->compressed_repr = compressed_repr;
            this->candidates = std::move(candidates);
    }

    int sudoku_id;
    std::string problem_file_path;
    std::string compressed_repr;
    int rank;
    std::vector<std::vector<int>> candidates;
};

typedef struct StepDifficulty {
    int sudoku_id;
    int num_of_easy_cells;
    int candidates_per_easy_cell;
    int average_number_of_candidates_per_cell;
} StepDifficulty;

typedef struct Difficulty {
    int sudoku_id;
    std::vector<StepDifficulty> steps_to_success;
    std::vector<StepDifficulty> steps_to_failure;
} Difficulty;

typedef struct SudokuSolution {
    int sudoku_id;
    std::vector<std::vector<int>> solution;
    Difficulty analysed_difficulty;
} SudokuSolution; 

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
            std::cout << ":SolveSudoku. Solution File Path is " << std::flush;
            std::cout.flush();
            std::this_thread::sleep_for(std::chrono::seconds(10));
        }

        Solver(ParsedSudoku* sudoku, SudokuSolution* solution, std::string output_folder) {
            sudoku = sudoku;
            solution_file_path = output_folder;
            solution = solution;
        }

        Solver() = delete;

    private:
        void WriteSudokuSolutionToFile() {

        }

        ParsedSudoku* sudoku = nullptr;
        std::string solution_file_path;
        SudokuSolution* solution;
};

class SolveSudokusInFolder {
    public:
        void AnalyseSudokus(std::string input_folder, std::string output_folder) {
            std::vector<std::string> sudoku_pathnames = GetSudokuPathnames(input_folder);
            //for (auto& path : sudoku_pathnames) {
                //std::cout << "sudoku pathnames are " << path << std::endl;
            //}
            std::vector<ParsedSudoku> parsed_sudokus = ParseSudokuProblems(sudoku_pathnames);
            std::vector<SudokuSolution> solutions(parsed_sudokus.size());
            std::cout << "starting threads" << std::endl;
            // For each Sudoku, we create a thread to solve and publish the result of the sudoku
            // as a CSV file into the output folder.
            std::vector<std::thread> sudoku_threads;
            for (int i = 0; i < parsed_sudokus.size(); i++) {
                Solver solver(&parsed_sudokus[i], &solutions[i], output_folder);
                sudoku_threads.push_back(std::thread(&Solver::SolveSudoku, &solver));
            }
            std::this_thread::sleep_for(std::chrono::seconds(10));
            for (auto& thread : sudoku_threads) {
                if(thread.joinable()) {
                    thread.join();
                } else {
                    std::cout << "Thread not joinable" << std::endl << std::flush;
                }
            }

            std::cout << "threads joined" << std::endl;
        }

    private:
        std::vector<std::string> GetSudokuPathnames(std::string folder_path) {
            // std::string docPath = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
            std::optional<std::vector<std::string>> files = ReadDirectory(folder_path);
            if (files != std::nullopt) {
                return *files;
            }
            return std::vector<std::string> {};
        }

        std::vector<ParsedSudoku> ParseSudokuProblems(std::vector<std::string>& sudoku_problem_files) {
            std::vector<ParsedSudoku> parsed_sudokus = {};
            for (auto& file : sudoku_problem_files) {
                //std::cout << "passing file " << file << std::endl;
                std::optional<ParsedSudoku> maybe_parsed_file = ParsedSudoku::Parse(file);
                if (maybe_parsed_file != std::nullopt) {
                    parsed_sudokus.push_back(*maybe_parsed_file);
                }
            }
            return parsed_sudokus;
        }
};

int main() {
    SolveSudokusInFolder solver;
    solver.AnalyseSudokus("parsed_sudokus", "sudoku_solutions");
}