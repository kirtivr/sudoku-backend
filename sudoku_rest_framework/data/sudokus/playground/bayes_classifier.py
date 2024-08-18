import sklearn
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
import os
import json
import glob

EASY_TRAINING_FOLDER = "solutions_easy_training/"
HARD_TRAINING_FOLDER = "solutions_hard_training/"
OUT_FOLDER = "labeled_by_classifier/"

def load_json(input_file):
    return json.load(input_file)

def add_label(sudoku_data, key, value):
    sudoku_data[key] = value

def write_to_file(sudoku, serial):
    outfile = os.path.join(OUT_FOLDER, str(serial) + '.csv')

    with open(outfile, "w+") as out:
        out.write(sudoku)

def extract_labels(sudoku_data, label_key):
    if len(sudoku_data) == 0:
        return []
    
    labels = []
    for sudoku_obj in sudoku_data:
        labels.append(sudoku_obj[label_key])

    return list(labels)

def extract_keys(sudoku_data, exclude_keys):
    if len(sudoku_data) == 0:
        return []
    
    keys = set()
    for sudoku_obj in sudoku_data:
        for key in sudoku_obj.keys():
            if key not in exclude_keys:
                keys.add(key)
    return list(keys)

def extract_row_values(sudoku_data, keys):
    if len(sudoku_data) == 0:
        return []
    
    row_values = []
    for sudoku_obj in sudoku_data:
        to_insert = []
        for key in keys:
            to_insert.append(sudoku_obj[key])
        row_values.append(to_insert)

    return row_values

def process(sudoku_data):
    label_names = ["easy", "medium", "hard"]
    labels = extract_labels(sudoku_data, "difficulty")
    print(f"labels = {labels}")
    exclude_key = ["difficulty", "unique_solutions", "original_sudoku", "id"]
    feature_names = extract_keys(sudoku_data, exclude_key)
    print(f"feature_names = {feature_names}")
    feature_rows = extract_row_values(sudoku_data, feature_names)
    print(f"first row = {feature_rows[0]}")
    train, test, train_labels, test_labels = train_test_split(feature_rows,labels,test_size=0.40,random_state=42)
    gaussian_nb = GaussianNB()
    NB_Clf = gaussian_nb.fit(train,train_labels)
    # Making predictions on the test set
    Preds_NBClf = NB_Clf.predict(test)
    # Print the predictions
    print(Preds_NBClf)
    print('Accuracy:',accuracy_score(test_labels, Preds_NBClf))

def parseAll(training_folder, label):
    parsed_training_set = []
    for training_file in glob.glob(training_folder + label + '/*.json'):
        #print("globbed file " + training_file)
        with open(training_file, "r") as input_file:
            sudoku_obj = load_json(input_file)
            add_label(sudoku_obj, "difficulty", label)
            parsed_training_set.append(sudoku_obj)

    process(parsed_training_set)

if __name__ == '__main__':
    parseAll("solutions_", "easy")
    #print(parseAll("solutions_", "medium"))
    #print(parseAll("solutions_", "hard"))
