import random
from stoch_runway_scheduler.populate import mutate_sequence

def test_mutate_sequence():
    k = 5 # length of sequences
    for i in range(3):
        seq = random.sample(range(20), k=k)
        for j in range(5):
            new_seq = mutate_sequence(seq)
            assert len(new_seq) == k
            assert set(new_seq) == set(seq)