from typing import List, Tuple
from copy import copy
import logging
import math
import random
import itertools
import time
import numpy as np
from .utils import FlightInfo, SequenceInfo

def Repopulate_VNS(GA_PopList: List[List[int]], GA_Info: List[SequenceInfo],
                    GA_PopSize: int, Opt_List: List[SequenceInfo], Opt_Size: int,
                    VNS_counter: int, VNS_limit: int, tot_mut: int):
    """
    Increases population of sequences GA_PopList (and associated GA_Info) to be at least GA_Popsize
    by modifying best current sequence.

    Information on all sequences (including old ones) is reset.

    Arguments:
    ---------
    GA_PopList: current list of sequences in population
    GA_Info: information associated with each sequence
    GA_PopSize: required size of population
    Opt_List: current list of infomation on most "optimal" sequences
    Opt_Size: required size of Opt_List
    VNS_counter: counter (m in paper)
    VNS_limit: m_mut in paper
    tot_mut: total number of times mutate_sequence has been run
    """

    # JF Question: this function seems to keep best previous `Opt_Size` sequences from GA_PopList and Opt_List
    #              and then create a completely new GA_PopList with GA_PopSize sequences.abs(
    #              Also, 
    #
    # JF Note: I don't fully understand logic behind VNS_counter and Opt_List

    # VNS_counter counts how many non-improving heuristic moves we've made since the last reset

    # Logging
    stepthrough_logger = logging.getLogger('stepthrough')
    step_summ_logger = logging.getLogger('step_summ')
    step_new_logger = logging.getLogger('step_new')

    stepthrough_logger.info('Repopulating...')
    step_summ_logger.info('Repopulating...')
    step_new_logger.info('Repopulating...')
    stepthrough_logger.info('Here are the sequences and their costs so far:')
    step_summ_logger.info('Here are the sequences and their costs so far:')
    for j, info in enumerate(GA_Info):
        flight_msg = str(j) + ',' + str(info.v) +',' + str(info.sequence) + '\n'
        stepthrough_logger.info(flight_msg)
        step_summ_logger.info(flight_msg)

    #AC_remaining = len(Arr_Pool) + len(Arr_NotReady)
    
    # JF Note: previously used the line below
    #          to get length of sequence, 
    #          but we should just be able to extract this from
    #          the current best sequence.
    #          This does create a new error in genetic however.
    # no_ACs = min(Max_LookAhead, AC_remaining)
    assert len(Opt_List) != 0
    no_ACs = len(Opt_List[0].sequence)

    queue_probs = [0] * no_ACs

    # Get best average cost (Best_in_pop, Best_in_opt)
    # This block is only used to increment VNS_counter
    # which keeps track of non-improving heuristic moves
    # This checks whether or not the best sequence
    # was already in the population last time we entered 4C
    # If it was we increment VNS_counter (m in pape)
    # otherwise we reset it to zero (steps a-b in 4C)
    # JF Question: I don't understand how the below achieves this aim
    if len(Opt_List) > 0 and len(GA_Info) > 0:

        GA_Info.sort(key=lambda x: x.v)
        Best_in_pop = GA_Info[0].v
        Opt_List.sort(key=lambda x: x.v)
        Best_in_opt = Opt_List[0].v

        if Best_in_opt < Best_in_pop:
            VNS_counter += 1
            step_new_logger.info('VNS_counter increased to %d', VNS_counter)
        else:
            VNS_counter = 0

    elif len(Opt_List) > 0 and len(GA_Info) == 0:
        VNS_counter += 1
        step_new_logger.info('VNS_counter increased to %d', VNS_counter)


    # Get best n=Opt_Size sequences from GA_Info and Opt_List
    New_Opt_List=[]
    for info in GA_Info:
        New_Opt_List.append(copy(info))
    for info in Opt_List:
        New_Opt_List.append(copy(info))

    New_Opt_List.sort(key=lambda x: x.v)

    while len(New_Opt_List) > Opt_Size:
        New_Opt_List.pop(len(New_Opt_List)-1)

    Best_Seq = New_Opt_List[0].sequence

    # Reset all sequence information in new list
    Opt_Seqs = []
    for info in New_Opt_List:
        info.n_traj = 0 # HMMM
        info.v = 0
        info.queue_probs = queue_probs
        info.w = 0
        Opt_Seqs.append(info.sequence[:])

    # Step c in 4C - apply mutate to ceate a new base sequence
    if VNS_counter >= VNS_limit:
        tot_mut += 1 # total mutations - only for logging purposes
        step_new_logger.info('Mutation performed!')
        # Perturb the optimal sequence
        Opt_Seq = Best_Seq[:] # Opt_List[0][0]
        Best_Seq = mutate_sequence(Opt_Seq)
        VNS_counter = 0         # Reset counter

    New_PopList = []

    # Apply heuristic move operator to generate enough sequences
    c = 0
    # JF Question: we should have len(Opt_Seqs) == Opt_Size because of code above (unless len(GA_Info) + len(Opt_List) < Opt_Size to begin with)
    while len(New_PopList) < GA_PopSize or len(Opt_Seqs) < Opt_Size:

        if c < 25: # no_ACs >= 6:
            new_seq = heuristic_move(Best_Seq) # Apply a change to the Best_Seq sequence
        elif c < 50: #else:
            new_seq = random.sample(Best_Seq, k=len(Best_Seq)) # random shuffle
        else:
            # In this case, give up finding new sequences
            break

        if new_seq not in New_PopList and new_seq not in Opt_Seqs:
            if len(New_PopList) < GA_PopSize:
                New_PopList.append(new_seq)
                c = 0
            else:
                Opt_Seqs.append(new_seq)
                New_Opt_List.append(SequenceInfo(new_seq, 0, 0, queue_probs, 0))
                c = 0
        else:
            c += 1

    GA_PopList = New_PopList[:]

    GA_Info = []
    for j in range(len(GA_PopList)):
        GA_Info.append(SequenceInfo(GA_PopList[j][:],0,0,queue_probs,0))

    GA_PopList_sorted=GA_PopList[:]
    GA_PopList_sorted.sort()

    Opt_List = New_Opt_List[:]

    # JF Question: I'm not really sure how this would happen
    if len(Opt_List) == 0:
        Opt_List = GA_Info[:]

    # More logging
    msg = 'Here is the new pop list after adding new sequences and sorting in sequence order:'
    stepthrough_logger.info(msg)
    stepthrough_logger.info(msg)
    for j in range(len(GA_PopList)):
        seq_msg = str(j)+','+str(GA_PopList_sorted[j])
        stepthrough_logger.info(seq_msg)
        step_summ_logger.info(seq_msg)

    step_new_logger.info('Here is the new Opt_List: ')
    for seq_info in Opt_List:
        step_new_logger.info(seq_info)
    step_new_logger.info('Here is the new Pop_List:')
    for seq in GA_PopList:
        step_new_logger.info(seq)

    # JF Question: why is Opt_Seq returned but not Opt_Seqs?
    return GA_PopList, GA_Info, Opt_List, VNS_counter, tot_mut

# 
def Populate(Ac_Info: List[FlightInfo], base_seq: List[int], 
            Arr_Pool: List[int], Arr_NotReady: List[int], 
            GA_PopSize: int, Max_SeqLength: int) -> Tuple[List[List[int]], List[SequenceInfo]]:
    """
    Creates a new population of sequences - used at beginning of algorithm and step 4A

    Arguments
    ---------
    Ac_Info: main list flight information/statuses
    base_seq: initial sequence from which others are derived (flight indices)
    Arr_Pool: indices of flights currently in pool
    Arr_NotReady: indices of flights not yet in the pool
    GA_PopSize: maximum number of sequences in population at one time (called S in paper)
    Max_SeqLength: maximum length of generated sequences (called l in paper)

    Returns
    -------
    GA_PopList: List of sequences
    GA_Info: List of information about sequences
    """

    # JF Question: should we add an assert that len(base_seq) <= Max_SeqLength? Yes

    # Number of flights not yet released from pool
    AC_remaining = len(Arr_Pool) + len(Arr_NotReady)
    no_ACs = min(Max_SeqLength, AC_remaining) # Lengths of all sequences that will be generated
    
    eta_list = [info.eta for AC, info in enumerate(Ac_Info)] # ETAs for each aircraft

    # In this case we extend sequence to be the right length
    # We extend by adding flights with closest ETA which are not already
    # in the sequence
    if len(base_seq) < no_ACs:
        base_seq = extend_sequence(base_seq, eta_list, Arr_Pool, Arr_NotReady, no_ACs)

    # Used in SimHeur as criterion for when flight is released from pool
    queue_probs = [0] * AC_remaining # redefined below

    # number of all possible sequences of remaining flights
    # only needed when there are only a small number of possible sequences left
    max_size = math.factorial(AC_remaining)

    # In this case, just add all remaining possible sequences
    if GA_PopSize >= max_size:
        remaining_seq = Arr_Pool + Arr_NotReady
        GA_PopList = [list(seq) for seq in itertools.permutations(remaining_seq)]
        # JF Question: Not clear why these all share same queue_probs list - this is not the case below where a new prob list is created for each sequence
        GA_Info = [SequenceInfo(poplist[:], 0, 0, queue_probs, 0) for poplist in GA_PopList]

    else:
        GA_PopList = []
        GA_Info = []
        no_seqs = 0
        c = 0 # Counter for number of times we duplicate the same sequence - this could be avoided by counting number of mutations possible with mutate_sequence (a quadratic number)

        while no_seqs < GA_PopSize:
            new_seq = heuristic_move(base_seq)

            # JF: could potentially make this more efficient by checking we are not
            # generating the same pos and pos2 in mutate_sequence
            # Would need some refactoring
            if new_seq not in GA_PopList:
                GA_PopList.append(new_seq)
                queue_probs = [0] * AC_remaining
                GA_Info.append(SequenceInfo(new_seq[:], 0, 0, queue_probs, 0))
                no_seqs += 1
            else:
                # If already in population of sequences,
                # Potentially create a whole new sequence by shuffling
                c += 1
                if c >= 100: # This condition means we retry above first until we have failed enough times
                    new_seq = random.sample(base_seq, k=len(base_seq)) # Create random permutation of sequence
                    if new_seq not in GA_PopList:
                        GA_PopList.append(new_seq)
                        queue_probs = [0] * AC_remaining
                        GA_Info.append(SequenceInfo(new_seq[:], 0, 0, queue_probs, 0))
                        no_seqs += 1

    return GA_PopList, GA_Info

def mutate_sequence(base_seq: List[int]) -> List[int]:
        no_ACs = len(base_seq)
        p = min(4, no_ACs) # no. of ACs to shuffle around
        r = no_ACs - p + 1 # no. of possible start positions

        L = r*(r+1)/2 # total weight
        triang_probs = tuple((r - i)/L for i in range(r))

        # Randomly select start of subsequence to shuffle
        # Use a discrete triangular distribution to put a priority on
        # moving flights near front of sequence
        pos = np.random.choice(range(r), p=triang_probs)

        new_seq = base_seq[:]
        seq_slice = random.sample(new_seq[pos:(pos+p)], k = p)
        new_seq[pos:(pos+p)] = seq_slice

        return new_seq

def heuristic_move(base_seq: List) -> List:
    """
    Takes a sequence and creates a new permuted one according algorithm described
    in Appendix A of "A New Simheuristic Approach for Stochastic Runway Scheduling"
    by Shone et al. (2024).

    Algorithm randomly selects an element in sequence and randomly moves it forward
    or backwards (by up to 3 indices), prioritising elements near the front of the sequence.
    
    In cases where either first or last element is selected, the function may return the same sequence.

    Arguments
    ---------
    base_seq: sequence to be permutated

    """
    no_ACs = len(base_seq)
    L = (no_ACs * (no_ACs + 1)) // 2
    triang_probs = tuple((no_ACs - i)/L for i in range(no_ACs))
    
    # Select a plane to randomly move in the sequence
    # Use a discrete triangular distribution to put a priority on
    # moving flights near front of sequence
    pos = np.random.choice(range(no_ACs), p=triang_probs)

    step_size = random.choice((1,2,3)) # No. of places to move (1, 2 or 3)
    z2 = random.random() # Determine whether to move up or down
    if z2 < 0.5: # Move backwards
        pos2 = max(0, pos - step_size)
    else: # Move forwards
        pos2 = min(no_ACs - 1, pos + step_size)

    new_seq = base_seq[:] # copies base sequence
    AC = new_seq.pop(pos) # Remove aircraft at position pos from sequence
    new_seq.insert(pos2, AC) # Reinsert aircraft at pos2
    return new_seq

def extend_sequence(base_seq: List[int], eta_list: List[float],
                    Arr_Pool: List[int], Arr_NotReady: List[int],
                    seq_length: int) -> List[int]:
    """
    Extends sequence to be of required length.

    Arguments
    ---------
    base_seq: sequence to extend
    eta_list: list of ETAs for all flights (ETA[i] is ETA of flight i)
    Arr_Pool: list of flights currently in Pool
    Arr_NotReady: list of flights not landed but not yet in pool
    seq_length: required length of sequence
    """
    assert len(Arr_Pool) + len(Arr_NotReady) >= seq_length

    new_seq = base_seq[:]
    ArrTime_Sorted = Arr_Pool + Arr_NotReady
    ArrTime_Sorted.sort(key = lambda ac: eta_list[ac])

    InScope_ACs = ArrTime_Sorted[:seq_length]
    for AC in new_seq:
        if AC not in InScope_ACs:
            # JF Question: removes flights that have landed? Shouldn't this already have been done? Possibly redundant
            # Perhaps we should have an assertion here as I don't think this should happen
            new_seq.remove(AC) # JF so we remove flights in base sequence?

    for ac in InScope_ACs:
        if ac not in new_seq:
            new_seq.append(ac)
        if len(new_seq) >= seq_length:
            break

    return new_seq