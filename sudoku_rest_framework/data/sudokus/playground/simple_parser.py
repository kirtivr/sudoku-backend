import os

OUT_FOLDER = "easy_sudokus"

def analyze_input_file(input_file):
    num_sudokus = 0
    lines = input_file.readlines()

    carry = 0
    for line in lines:
        line = line.strip()
        length = len(line)
        sudokus = (length + carry) // 81;
        carry = (length + carry) % 81;
        num_sudokus += sudokus

    return num_sudokus

def chunk_input_lines(input_file):
    lines = input_file.readlines()

    chunks = []
    carry = ''
    for line in lines:
        line = line.strip()
        if len(line) == 0:
            continue
        effective_line = carry + line
        #print(line)
        if len(effective_line) < 81:
            print("Unexpected input length " + str(len(effective_line)))
            break
        chunks.append(effective_line[:81])
        length = len(effective_line)
        right_unused = length - 81
        carry = line[81:]

    return chunks

def generate_sudoku_string(line):
    if len(line) != 81:
        print("Invalid sudoku line")

    sudoku_string = ''
    first_st = True
    for row in range(0, 9):
        for column in range(0, 9):
            c = line[row * 9 + column]
            if c == '.' or c == '0':
                continue
            if not first_st:
                sudoku_string += ','
            else:
                first_st = False
            sudoku_string = sudoku_string + c + str(row) + str(column)

    return sudoku_string

def write_to_file(sudoku, serial):
    outfile = os.path.join(OUT_FOLDER, str(serial) + '.csv')

    with open(outfile, "w+") as out:
        out.write(sudoku)

def generate_sudoku_strings(input_file):
    sudoku_text = chunk_input_lines(input_file)
    num = 0
    for st in sudoku_text:
        compressed = generate_sudoku_string(st)
        num += 1
        write_to_file(compressed, num)

def parseAll(sudoku_file):
    with open(sudoku_file, "r") as input_file:
        #num_sudokus = analyze_input_file(input_file)
        print(generate_sudoku_strings(input_file))

if __name__ == '__main__':
    print(parseAll("easy_sudokus_"))
