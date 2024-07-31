/*
 * C# Program to solve a sudoku and analyse:
 *   1. Whether the solution is unique.
 *   2. How difficult is the sudoku.
 */
 
using System;
using System.IO;
using System.Threading;

public static class NullableExtensions
{
    public static void Deconstruct<T>(
        this T? nullable,
        out bool hasValue,
        out T value) where T : struct
    {
        hasValue = nullable.HasValue;
        value = nullable.GetValueOrDefault();
    }
}

// This class is thread-compatible.
// It is only ever used by one thread.
class ParsedSudoku
{
    private int sudoku_id;
    private string problem_file_path;
    private string compressed_repr;
    private int rank;
    private List<List<int>> candidates;

    private ParsedSudoku(
        int sudoku_id, string problem_file_path,
        string compressed_repr, int rank,
        List<List<int>> candidates) {
            sudoku_id = sudoku_id;
            problem_file_path = problem_file_path;
            compressed_repr = compressed_repr;
            candidates = candidates.Clone();
    }

    static NullableExtensions<ParsedSudoku> ParseSudoku(string problem_file_path) {
        string givens;
        int sudoku_id;
        try {
            StreamReader sr = new StreamReader(problem_file_path);
            givens = sr.ReadLine();
            sr.Close();

            if (problem_file_path.Substring(problem_file_path.Length - 4) == ".csv") {
                sudoku_id = int(problem_file_path.Substring(0, problem_file_path.Length - 4));                
            } else {
                sudoku_id = int(problem_file_path);
            }
        } catch(Exception e) {
            Console.WriteLine("Exception: " + e.Message);
            return null;
        }

        // Givens is one string of comma separated:
        // rowcolumnvalue strings.
        // First find the Rank of the sudoku by counting the
        // comma separated values.
        // Then populate the candidate list.
        string[] rowcolumnvalue = givens.Split(',');
        int rank = rowcolumnvalue.Length >= 17 ? 3 : 2;

        int index = 0;
        List<List<int>> candidates = new List<List<int>>();
        for (int i = 0; i < rank * rank; i++) {
            candidates.Add(new List<int>());
            for (int j = 0; j < rank * rank; j++) {
                candidates[i].Add(0);
            }
        }

        for (string rcv in rowcolumnvalue) {
            int row = int(rcv[0]);
            int col = int(rcv[1]);
            int val = int(rcv[2]);
            candidates[row][col] = val;
        }

        return NullableExtensions<ParsedSudoku>(
            sudoku_id, problem_file_path, givens, rank, candidates
        );
    }
}

class StepDifficulty {
    public int sudoku_id;
    public int num_of_easy_cells;
    public int candidates_per_easy_cell;
    public int average_number_of_candidates_per_cell;
}

class Difficulty {
    public int sudoku_id;
    public List<StepDifficulty> steps_to_success;
    public List<StepDifficulty> steps_to_failure;
}

// This class is thread-compatible.
// It is only instantiated by one thread.
class SudokuSolution {
    public int sudoku_id;
    public int[,] solution;
    public Difficulty analysed_difficulty;
} 

class Solver
{
    // Strategy
    // 1. Given populated candidates[][], sort the bags by lowest number of candidates.
    //    For the bags, one by one try each equal choice.
    //    If the attempt eventually succeeds, add to the steps_to_success List.
    //    Else add to steps_to_failure list.

    //  2. Update, re-analyze and re-sort candidates[][].
    public static void SolveSudoku() {
        Console.WriteLine("In SolveSudoku. Solution File Path is " + solution_file_path);
    }

    private void WriteSudokuSolutionToFile() {

    }

    public Solver(ParsedSudoku sudoku, string output_folder) {
        this.sudoku = sudoku;
        this.solution_file_path = output_folder;
    }

    private ParsedSudoku sudoku;
    private solution_file_path;
    private SudokuSolution solution;
}

class SolveSudokusInFolder {
    private List<string> GetSudokuPathnames(string folder_path) {
        // string docPath = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
        try {
            var files = from file in Directory.EnumerateFiles(docPath, "*.txt", SearchOption.AllDirectories) select file;
        } catch (Exception e) {
            Console.WriteLine("Error: " + e.Message);
        }
        return new List<string>(files);
    }

    private List<ParsedSudoku> ParseSudokuProblems(List<string> sudoku_problem_files) {
        List<ParsedSudoku> parsed_sudokus = new List<ParsedSudoku>();
        foreach(string file in sudoku_problem_files) {
            NullableExtensions<ParsedSudoku> maybe_parsed_file = ParsedSudoku::ParseSudoku(file);
            if (maybe_parsed_file.hasValue) {
                parsed_sudokus.Add(maybe_parsed_file.GetValueOrDefault());
            }
        }
        return parsed_sudokus;
    }

    public void AnalyseSudokus(string input_folder, string output_folder) {
        List<string> sudoku_pathnames = GetSudokuPathnames(input_folder);
        List<ParsedSudoku> parsed_sudokus = ParseSudokuProblems(sudoku_pathnames);

        // For each Sudoku, we create a thread to solve and publish the result of the sudoku
        // as a CSV file into the output folder.
        List<Thread> sudoku_threads;
        foreach (ref sudoku in parsed_sudokus) {
            Solver solver(sudoku, output_folder);
            sudoku_threads.Add(new Thread(solver::SolveSudoku()));
        }

        foreach (ref thread: sudoku_threads) {
            thread.Start();            
        }
    }
}

class Runner {
    static void Main()
    {
        SolveSudokusInFolder solver;
        solver.AnalyseSudokus("parsed_sudokus", "sudoku_solutions");
    }
}