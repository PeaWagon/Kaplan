
import os

from kaplan.gac import run_kaplan

# directory for this test file
test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'testfiles')

def test_run_kaplan():
    test_ga = os.path.join(test_dir, "example_ga_input_file.txt")
    test_mol = os.path.join(test_dir, "example_mol_input_file.txt")
    run_kaplan(test_ga, test_mol)
    test_ga2 = os.path.join(test_dir, "example2_ga_input_file.txt")
    test_mol2 = os.path.join(test_dir, "example2_mol_input_file.txt")
    run_kaplan(test_ga2, test_mol2)

if __name__ == "__main__":
    test_run_kaplan()
