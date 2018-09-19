
import os, sys

import src


print(os.getcwd())

def test_mutate():
    src.mutations.mutate([1,3,4], 1, 2)
#    assert not problem

if __name__ == "__main__":
    test_mutate()

