
import kaplan

def test_mutate():
    # test empty input_list
    result1 = kaplan.mutations.mutate([], 1, 5)
    assert result1 == []
    # test no mutations
    result2 = kaplan.mutations.mutate([1,2,3], 0, 0)
    assert result2 == [1,2,3]
    # test one mutation
    result3 = kaplan.mutations.mutate([1,2,3], 1, 1)
    assert len(result3) == 3
    assert sum((result3[0] == 1, result3[1] == 2, result3[2] == 3)) >= 2
    # test mutations equal to length of input_list
    result4 = kaplan.mutations.mutate([4,5,6], 3, 3)
    assert len(result4) == 3

if __name__ == "__main__":
    test_mutate()

