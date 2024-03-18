from typing import List
import math
import random
import itertools
import time

def Repopulate_VNS(GA_PopList,GA_Info,Arr_Pool,Arr_NotReady,GA_PopSize,Opt_Seq,OptCost,Opt_List,Opt_Size, Max_LookAhead: int ,VNS_counter,VNS_limit,tot_mut, stepthrough: int, step_summ: int, step_new: int):

    #print('Repopulating')

    start_time=time.time()

    if stepthrough==1:
        st.write('Repopulating...'+'\n')
        st.write('Here are the sequences and their costs so far:'+'\n')
        for j in range(len(GA_Info)):
            st.write(str(j)+','+str(GA_Info[j][2])+','+str(GA_Info[j][0])+'\n')
        st.write('\n')
    if step_summ==1:
        st2.write('Repopulating...'+'\n')
        st2.write('Here are the sequences and their costs so far:'+'\n')
        for j in range(len(GA_Info)):
            st2.write(str(j)+','+str(GA_Info[j][2])+','+str(GA_Info[j][0])+'\n')
        st2.write('\n')
    if step_new==1:
        st3.write('Repopulating...'+'\n')

    AC_remaining=len(Arr_Pool)+len(Arr_NotReady)
    no_ACs=min(Max_LookAhead,AC_remaining)

    queue_probs=[0]*no_ACs

    # print('len(Opt_List): '+str(len(Opt_List)))
    # print('len(GA_Info): '+str(len(GA_Info)))

    if len(Opt_List)>0 and len(GA_Info)>0:

        GA_Info.sort(key=lambda x: x[2])
        Best_in_pop=GA_Info[0][2]
        Opt_List.sort(key=lambda x: x[2])
        Best_in_opt=Opt_List[0][2]

        if Best_in_opt<Best_in_pop:
            VNS_counter+=1
            if step_new==1:
                st3.write('VNS_counter increased to '+str(VNS_counter)+'\n')
        else:
            VNS_counter=0

    elif len(Opt_List)>0 and len(GA_Info)==0:

        VNS_counter+=1
        if step_new==1:
            st3.write('VNS_counter increased to '+str(VNS_counter)+'\n')

    New_Opt_List=[]
    for j in range(len(GA_Info)):
        New_Opt_List.append(GA_Info[j][:])
    for j in range(len(Opt_List)):
        New_Opt_List.append(Opt_List[j][:])

    New_Opt_List.sort(key=lambda x: x[2])

    while len(New_Opt_List)>Opt_Size:
        New_Opt_List.pop(len(New_Opt_List)-1)

    Best_Seq=New_Opt_List[0][0][:]

    Opt_Seqs=[]
    for i in range(len(New_Opt_List)):
        New_Opt_List[i][1]=0 #HMMM
        New_Opt_List[i][2]=0
        New_Opt_List[i][3]=queue_probs
        New_Opt_List[i][4]=0
        #print('New_Opt_List[i][0]: '+str(New_Opt_List[i][0]))
        Opt_Seqs.append(New_Opt_List[i][0][:])

    # GA_Info.sort(key=lambda x: x[2])
    # Opt_List.sort(key=lambda x: x[2])
    # #Tabu_List.sort(key=lambda x: x[2])

    # AC_remaining=len(Arr_Pool)+len(Arr_NotReady)
    # no_ACs=min(Max_LookAhead,AC_remaining)

    # queue_probs=[0]*no_ACs

    # if len(Opt_List)<Opt_Size:
    # 	j=0
    # 	while len(Opt_List)<Opt_Size and j<=len(GA_Info)-1:
    # 		Opt_List.append([GA_Info[j][0][:],0,0,queue_probs,0])
    # 	j+=1
    # 	VNS_counter=0

    # else:

    # 	max_j=min(len(GA_Info),len(Opt_List))

    # 	j=0
    # 	while j<=max_j:

    # 	Best_Seq=GA_Info[0][0][:]
    # 	BestCost=GA_Info[0][2]
    # 	OptCost=Opt_List[len(Opt_List)-1][2]

    # 	if BestCost<OptCost:
    # 		if len(Opt_List)==Opt_Size:
    # 			Opt_List.remove(Opt_List[Opt_Size-1])
    # 		Opt_List.append([Best_Seq[:],0,0,queue_probs,0])
    # 		VNS_counter=0
    # 	else:
    # 		VNS_counter+=1

    # for i in range(len(Opt_List)):
    # 	Opt_List[i][1]=0 #HMMM

    if VNS_counter>=VNS_limit:

        tot_mut+=1 #total mutations
        if step_new==1:
            st3.write('Mutation performed!'+'\n')
        #print('tot_mut: '+str(tot_mut))

        #Perturb the optimal sequence
        Opt_Seq=Best_Seq[:] #Opt_List[0][0]

        perm_size=min(4,no_ACs) #no. of ACs to shuffle around
        no_start_pos=no_ACs-perm_size+1 #no. of possible start positions

        triangle_dist_size=no_start_pos*(no_start_pos+1)/2

        z=random.random()*triangle_dist_size
        totp=0
        for ii in range(no_start_pos):
            if z<(totp+no_start_pos-ii):
                pos=ii
                break
            totp+=no_start_pos-ii

        remove_perm=[]
        for ii in range(perm_size):
            AC=Opt_Seq[pos]
            Opt_Seq.remove(AC)
            remove_perm.append(AC)

        old_perm=remove_perm[:]
        random.shuffle(remove_perm)

        for ii in range(perm_size):
            AC=remove_perm[ii]
            Opt_Seq.insert(pos,AC)

        Best_Seq=Opt_Seq[:]

        VNS_counter=0

    New_PopList=[]

    #no_ACs=min(Max_LookAhead,AC_remaining)

    c=0
    while len(New_PopList)<GA_PopSize or len(Opt_Seqs)<Opt_Size:

        if c<25: #no_ACs>=6:

            #Apply a change to the Best_Seq sequence

            triangle_dist_size=no_ACs*(no_ACs+1)/2
            z=random.random()*triangle_dist_size
            totp=0
            for ii in range(no_ACs):
                if z<(totp+no_ACs-ii):
                    pos=ii
                    break
                totp+=no_ACs-ii

            new_seq=Best_Seq[:]

            AC=new_seq[pos]
            new_seq.remove(AC)
            z1=int(random.random()*3)+1 #no. of places to move
            z2=random.random() #determine whether to move up or down
            if z2<0.5:
                z1=min(z1,pos)
                pos2=pos-z1
            else:
                z1=min(z1,no_ACs-1-pos)
                pos2=pos+z1

            new_seq.insert(pos2,AC)

        elif c<50: #else:

            new_seq=Best_Seq[:]
            random.shuffle(new_seq)

        else:

            break

        # if AC_remaining==30 and len(New_PopList)==0:
        # 	new_seq=[0,1,3,6,8,7,10,12,2,4,9,11,5,13,14,15,20,16,17,21,22,24,28,23,25,26,27,29,18,19]

        if new_seq not in New_PopList and new_seq not in Opt_Seqs:
            if len(New_PopList)<GA_PopSize:
                New_PopList.append(new_seq)
                c=0
            else:
                Opt_Seqs.append(new_seq)
                New_Opt_List.append([new_seq,0,0,queue_probs,0])
                c=0
        else:
            c+=1

    GA_PopList=New_PopList[:]

    GA_Info=[]
    for j in range(len(GA_PopList)):
        GA_Info.append([GA_PopList[j][:],0,0,queue_probs,0])

    # GA_PopList_sorted=GA_PopList[:]
    # GA_PopList_sorted.sort()

    Opt_List=New_Opt_List[:]

    if len(Opt_List)==0:
        Opt_List=GA_Info[:]

    if stepthrough==1:
        st.write('Here is the new pop list after adding new sequences and sorting in sequence order:'+'\n')
        for j in range(len(GA_PopList)):
            st.write(str(j)+','+str(GA_PopList_sorted[j])+'\n')
        st.write('\n')
    if step_summ==1:
        st2.write('Here is the new pop list after adding new sequences and sorting in sequence order:'+'\n')
        for j in range(len(GA_PopList)):
            st2.write(str(j)+','+str(GA_PopList_sorted[j])+'\n')
        st2.write('\n')
    if step_new==1:
        st3.write('Here is the new Opt_List: '+'\n')
        for j in range(len(Opt_List)):
            st3.write(str(Opt_List[j])+'\n')
        st3.write('\n')
        st3.write('Here is the new Pop_List: '+'\n')
        for j in range(len(GA_PopList)):
            st3.write(str(GA_PopList[j])+'\n')
        st3.write('\n')

    # run_time=(time.time()-start_time)/conv_factor
    # #run_time=0

    #run_time=fixed_repop_elap/conv_factor

    return GA_PopList,GA_Info,Opt_Seq,OptCost,Opt_List,VNS_counter,tot_mut


# 
def Populate(Ac_Info: List, base_seq: List[int], Arr_Pool: List[int] ,Arr_NotReady: List[int], GA_PopSize: int, Max_SeqLength: int, stepthrough: int, step_summ: int, step_new: int):
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
    """

    #print('Populating')
    # For calculating computational time of function
    start_time=time.time()

    if stepthrough==1:
        st.write('\n'+'Populating...'+'\n')
    if step_summ==1:
        st2.write('\n'+'Populating...'+'\n')
    if step_new==1:
        st3.write('\n'+'Populating...'+'\n')

    # Number of flights not yet released from pool
    AC_remaining = len(Arr_Pool)+len(Arr_NotReady)

    no_ACs = min(Max_SeqLength, AC_remaining) # Lengths of all sequences that will be generated
    
    # In this case we extend sequence to be the right length
    # We extend by adding flights with closest ETA which are not already
    # in the sequence
    if len(base_seq) < no_ACs: 
        ArrTime_Sorted=[] # sequence of flights sorted by current ETA
        for AC in Arr_Pool:
            # AC_Info[3] = latest ETA
            ArrTime_Sorted.append([Ac_Info[AC][3], AC])
        for AC in Arr_NotReady:
            ArrTime_Sorted.append([Ac_Info[AC][3], AC])
        ArrTime_Sorted.sort(key=lambda x: x[0])

        InScope_ACs=[]
        i=0
        while len(InScope_ACs)<no_ACs:
            InScope_ACs.append(ArrTime_Sorted[i][1])
            i+=1
        for AC in base_seq:
            if AC not in InScope_ACs:
                base_seq.remove(AC)
        i=0
        while len(base_seq)<no_ACs:
            if InScope_ACs[i] not in base_seq:
                base_seq.append(InScope_ACs[i])
            i+=1

    # Used in SimHeur as criterion for when flight is released from pool
    queue_probs = [0]*AC_remaining # redefined below

    # number of all possible sequences of remaining flights
    # only needed when there are only a small number of possible sequences left
    max_size = math.factorial(AC_remaining)

    # In this case, just add all remaining possible sequences
    if GA_PopSize>=max_size:
        remaining_seq=Arr_Pool[:]
        for AC in Arr_NotReady:
            remaining_seq.append(AC)
        GA_Tuples=list(itertools.permutations(remaining_seq))
        GA_PopList=[]
        for j in range(len(GA_Tuples)):
            GA_PopList.append(list(GA_Tuples[j]))

        GA_Info=[]
        for j in range(len(GA_PopList)):
            GA_Info.append([GA_PopList[j][:],0,0,queue_probs,0])
    else:
        GA_PopList=[]
        GA_Info=[]
        no_seqs=0
        c=0
        chk=0
        # JF: I have corrected the line below - is this right? Max_LookAhead was previously
        # being extracted from a global variable, but it should be a function argument
        # no_ACs=min(Max_LookAhead,AC_remaining)
        no_ACs=min(Max_SeqLength,AC_remaining)

        while no_seqs < GA_PopSize:

            # Select a plane to randomly move in the sequence
            # Use a discrete triangular distribution to put a priority on
            # moving flights near front of sequence
            # See Appendix A in Paper


            # Could possible replace below with numpy.random.choice
            # used to normalise probabilities
            triangle_dist_size = no_ACs*(no_ACs+1)/2 # (1 + 2 + ...+ no_ACs)

            # Sampling a value from triangular dist - pick a point in the sequence
            z=random.random()*triangle_dist_size
            totp=0
            for ii in range(no_ACs):
                if z<(totp+no_ACs-ii):
                    pos=ii
                    break
                totp+=no_ACs-ii

            new_seq=base_seq[:] # copies base sequence

            AC=new_seq[pos] # gets aircraft at randomly selected point
            new_seq.remove(AC)
            z1=int(random.random()*3)+1 # no. of places to move
            z2=random.random() #determine whether to move up or down
            if z2<0.5:
                z1=min(z1,pos)
                pos2=pos-z1
            else:
                z1=min(z1,no_ACs-1-pos)
                pos2=pos+z1
            # Reinserts aircraft in another place
            new_seq.insert(pos2,AC)

            if new_seq not in GA_PopList:
                GA_PopList.append(new_seq)
                queue_probs=[0]*AC_remaining # remove?
                GA_Info.append([new_seq[:],0,0,queue_probs,0])
                no_seqs+=1
            else:
                # If already in population of sequences,
                # Potentially create a whole new sequence by shuffling
                c+=1
                if c>=100: # This condition means we retry above first until we have failed enough times

                    new_seq=base_seq[:]
                    random.shuffle(new_seq)
                    if new_seq not in GA_PopList:
                        GA_PopList.append(new_seq)
                        queue_probs=[0]*AC_remaining # remove?
                        GA_Info.append([new_seq[:],0,0,queue_probs,0])
                        no_seqs+=1
                    #chk=1

    # Only used for logging purposes
    GA_PopList_sorted=GA_PopList[:]
    GA_PopList_sorted.sort()

    if stepthrough==1:
        for j in range(len(GA_PopList)):
            st.write(str(j)+','+str(GA_PopList_sorted[j])+'\n')
        st.write('\n')
    if step_summ==1:
        for j in range(len(GA_PopList)):
            st2.write(str(j)+','+str(GA_PopList_sorted[j])+'\n')
        st2.write('\n')
    if step_new==1:
        # st3.write('Here is the new Opt_List: '+'\n')
        # for j in range(len(Opt_List)):
        # 	st3.write(str(Opt_List[j])+'\n')
        # st3.write('\n')
        st3.write('Here is the new GA_PopList: '+'\n')
        for j in range(len(GA_PopList)):
            st3.write(str(GA_PopList[j])+'\n')
        st3.write('\n')

    # Tabu_List=[base_seq]

    # run_time=(time.time()-start_time)/conv_factor

    #run_time=fixed_pop_elap/conv_factor

    return GA_PopList,GA_Info