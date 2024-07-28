/*
 * C# Program to solve a sudoku and analyse:
 *   1. Whether the solution is unique.
 *   2. How difficult is the sudoku.
 */
 
using System;
using System.IO;

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

    static ParsedSudoku ParseSudoku(string problem_file_path) {
        string givens;
        int sudoku_id;
        try {
            StreamReader sr = new StreamReader(problem_file_path);
            givens = sr.ReadLine();
            sr.Close();
            Console.ReadLine();

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
        List<int> rcv_triplets;
        int rank = 2;        

        for (string rcv in rowcolumnvalue) {
            int row = int(rcv[0]);
            int col = int(rcv[1]);
            int val = int(rcv[2]);
            rcv_triplets.Add(row);
            rcv_triplets.Add(col);
            rcv_triplets.Add(val);
            if (row >= 4 || col >= 4) {
                rank = 3;
            }
        }

        List<List<int>> candidates = new List<List<int>>();
        for (int i = 0; i < rank * rank; i++) {
            candidates.Add(new List<int>());
            for (int j = 0; j < rank * rank; j++) {
                candidates[i].Add(0);
            }
        }

        for (int i = 0; i < rcv_triplets.Length; i+=3) {
            int row = rcv_triplets[i];
            int col = rcv_triplets[i + 1];
            int val = rcv_triplets[i + 2];
            candidates[i][j] = val;
        }

        return ParsedSudoku(
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
    public int rank;
    public int[,] solution;
    public Difficulty analysed_difficulty;
} 

class SudokuSolutionFile
{
    public int sudoku_id;
    public string solution_file_path;
    public SudokuSolution solution;
}

class Solver
{
    public SudokuSolution solve_sudoku() {

    }

    static void Main()
    {
        ParsedSudoku parsedSudoku = ParsedSudoku::ParseSudoku(
            "parsed_sudokus/001.csv"
        );

        
    }
}