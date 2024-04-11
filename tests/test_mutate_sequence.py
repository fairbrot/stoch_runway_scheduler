import random
from stoch_runway_scheduler.populate import heuristic_move, mutate_sequence

def test_heuristic_move():
    k = 5 # length of sequences
    for i in range(3):
        seq = random.sample(range(20), k=k)
        for j in range(5):
            new_seq = heuristic_move(seq)
            assert len(new_seq) == k
            assert set(new_seq) == set(seq)


def test_mutate_sequence():
    k = 8 # length of sequences
    for i in range(3):
        seq = random.sample(range(20), k=k)
        for j in range(5):
            new_seq = mutate_sequence(seq)
            assert len(new_seq) == k
            assert set(new_seq) == set(seq)
            assert sum(i != j for i, j in zip(seq, new_seq)) <= min(4, k)