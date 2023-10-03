# JF: there are two versions of this function - this one has more arguments than the second and does not seem to be used so commenting out
# def Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ArrTime_Sorted,ServTime,Trav_Time,output,tm,stored_prev_class,queue_complete,AC_remaining,no_scenarios,use_determ):

# 	totcost=0
# 	ee=0

# 	# f2.write('Time '+','+str(tm)+'\n')
# 	# f2.write('Arr Pool is '+','+str(Arr_Pool)+'\n')
# 	# f2.write('Arr NotReady is '+','+str(Arr_NotReady)+'\n')
# 	# f2.write('ArrTime_Sorted is '+','+str(ArrTime_Sorted)+'\n')
# 	# f2.write('Stored prev class is '+str(stored_prev_class)+'\n')
# 	# f2.write('Testing the sequence '+','+str(Anneal_Seq)+'\n'+'\n')

# 	for sc in range(no_scenarios):

# 		ArrTime_sc=ArrTime[sc]
# 		ArrTime_Sorted_sc=ArrTime_Sorted[sc]
# 		ServTime_sc=ServTime[sc]
# 		Trav_Time_sc=Trav_Time[sc]

# 		anneal_prev_class=stored_prev_class
# 		anneal_queue_complete=queue_complete[sc]
# 		anneal_tm=tm
# 		anneal_weather_state=weather_state

# 		for i in range(AC_remaining):
# 			AC=Anneal_Seq[i]
# 			release_time=max(ArrTime_sc[AC],anneal_tm)
# 			cur_class=Ac_Info[AC][1]

# 			if no_scenarios==1 and use_determ==1:
# 				t1=release_time+tau
# 				t2=anneal_queue_complete+Time_Sep[anneal_prev_class][cur_class]/60
# 				AC_FinishTime=max(t1,t2)
# 			else:
# 				if NormalApprox==0:
# 					anneal_weather_state=weather(release_time,wlb,wub)
# 					AC_FinishTime,straight_into_service=GetServTime(Trav_Time_sc[AC],ServTime_sc[AC],release_time,anneal_prev_class,cur_class,anneal_queue_complete,k,ee,anneal_weather_state)
# 				else:
# 					AC_FinishTime,straight_into_service=Normal_GetServ(release_time,anneal_prev_class,anneal_class,anneal_queue_complete,weather_state, Time_Sep, normcdf, w_rho, k)

# 			#f2.write('AC '+str(AC)+','+'anneal_tm '+str(anneal_tm)+','+'release_time'+str(release_time)+','+'anneal_queue_complete'+str(anneal_queue_complete)+','+'anneal_prev_class'+str(anneal_prev_class)+','+'cur_class'+str(cur_class)+','+'Time Sep '+str(Time_Sep[anneal_prev_class][cur_class]/60)+','+'AC_FinishTime '+str(AC_FinishTime)+'\n')

# 			ps_time=Ac_Info[AC][18]
# 			totcost+=getcost(ps_time,ArrTime_sc[AC],Trav_Time_sc[AC],AC_FinishTime,Ac_Info[AC][10],thres1,thres2, lam1, lam2)
# 			anneal_tm=release_time
# 			anneal_prev_class=cur_class
# 			anneal_queue_complete=AC_FinishTime

# 		#f2.write('\n')

# 	return totcost/no_scenarios

# JF: This isn't used anywhere at the moment
# def GA_Heur(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,GA_PopSize,wlb_tm,wub_tm):

#     start_time=time.time()

#     tm=0

#     Arr_Pool=[]
#     Arr_NotReady=[]
#     for AC in range(NoA):
#         Arr_NotReady.append(AC)

#     totserv=0
#     totcost=0
#     prev_class=4
#     queue_complete=0
#     weather_state=0

#     base_seq=[0]*NoA
#     for i in range(NoA):
#         base_seq[i]=ArrTime_Sorted[i][1]

#     # JF: to check with Rob - NoA is called Max_SeqLength in function definition
#     GA_PopList,GA_Info=Populate(Ac_Info, base_seq,Arr_Pool,Arr_NotReady,GA_PopSize,NoA,stepthrough, step_summ, step_new)

#     #GA_Heur_iters=100
    
#     iter_no=0
#     runtime_chk=0
#     while runtime_chk==0:

#         for seq in range(len(GA_PopList)):

#             perm=GA_PopList[seq]
#             perm_cost=0
#             latest_tm=0
#             perm_prev_class=4
#             perm_queue_complete=0
#             perm_weather_state=0
#             j=0

#             while j<NoA:

#                 AC=perm[j]
#                 Ac_Infoi=Ac_Info[AC]
#                 release_time=max(latest_tm,ArrTime[AC][0])
#                 trav_time=Ac_Infoi[6]
#                 perm_class=Ac_Infoi[1]
#                 begin_serv=max(release_time,perm_queue_complete)
#                 perm_weather_state=weather(release_time,wlb_tm,wub_tm) #weather(begin_serv,wlb_tm,wub_tm)

#                 if NormalApprox==0:
#                     if perm_weather_state==1:
#                         ws=1/w_rho
#                     else:
#                         ws=1
#                     rate=ws*k/(Time_Sep[perm_prev_class][perm_class]/60)
#                     serv=ServTime[AC]/rate
#                     # serv=0
#                     # for m in range(k):
#                     # 	serv+=(-1/rate)*math.log(ServTime[AC][m])
#                 else:
#                     Mn=Time_Sep[perm_prev_class][perm_class]/60
#                     if perm_weather_state==1:
#                         Mn*=w_rho
#                     SD=math.sqrt(Mn**2/k)
#                     z=Ac_Info[AC][7]*SD+Mn
#                     # u=int(z*10000)
#                     # serv=normcdf[u]*SD+Mn

#                 t1=release_time+trav_time
#                 t2=perm_queue_complete+serv
#                 perm_queue_complete=max(t1,t2)

#                 perm_cost+=getcost(Ac_Infoi[18],ArrTime[AC][0],trav_time,perm_queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2)

#                 latest_tm=release_time
#                 perm_prev_class=perm_class

#                 j+=1

#             GA_Info[seq][2]=perm_cost

#         GA_Info.sort(key=lambda x: x[2])
#         totcost=GA_Info[0][2]

#         # JF - the Repopulate function doesn't exist
#         GA_PopList,GA_Info,run_time,Tabu_List,Opt_Seq,OptCost,Opt_List=Repopulate(GA_PopList,GA_Info,Arr_Pool,Arr_NotReady,GA_PopSize,Tabu_List,Tabu_Size,Opt_Seq,OptCost,Opt_List,Opt_Size,stepthrough, step_summ, step_new)

#         iter_no+=1
#         curr_time=time.time()
#         if (curr_time-start_time)/conv_factor>=S*t:
#             runtime_chk=1

#     return totcost,iter_no

# # JF: This isn't used anywhere right now
# def SA_Heur(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm):

#     start_time=time.time()

#     tm=0

#     totserv=0
#     totcost=0
#     prev_class=4
#     queue_complete=0
#     weather_state=0

#     Anneal_Seq=[0]*NoA
#     for i in range(NoA):
#         Anneal_Seq[i]=ArrTime_Sorted[i][1]

#     NewCost=Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)
#     OldCost=NewCost
#     OptCost=NewCost
#     Opt_Seq=Anneal_Seq[:]

#     T=1000
#     iter_no=0
#     runtime_chk=0
#     AC_remaining=NoA
#     c=0

#     while runtime_chk==0:

#         #Move one plane in the sequence
#         no_ACs=AC_remaining
#         triangle_dist_size=no_ACs*(no_ACs+1)/2

#         z=random.random()*triangle_dist_size
#         totp=0
#         for ii in range(no_ACs):
#             if z<(totp+no_ACs-ii):
#                 pos=ii
#                 break
#             totp+=no_ACs-ii

#         AC=Anneal_Seq[pos]
#         Anneal_Seq.remove(AC)
#         z1=int(random.random()*3)+1 #no. of places to move
#         z2=random.random() #determine whether to move up or down
#         if z2<0.5:
#             z1=min(z1,pos)
#             pos2=pos-z1
#         else:
#             z1=min(z1,AC_remaining-1-pos)
#             pos2=pos+z1

#         Anneal_Seq.insert(pos2,AC)

#         #End of mutation part

#         NewCost=Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)

#         if NewCost<OldCost:
#             if iter_no==0 or NewCost<OptCost:
#                 OptCost=NewCost
#                 Opt_Seq=Anneal_Seq[:]
#             OldCost=NewCost
#             c=0
#         else:
#             temp=math.exp(-t/T)
#             #temp=0
#             z=random.random()
#             if z<temp:
#                 OldCost=NewCost #accepted move anyway
#                 c=0
#             else:
#                 AC=Anneal_Seq[pos2]
#                 Anneal_Seq.remove(AC)
#                 Anneal_Seq.insert(pos,AC)
#                 c+=1

#         iter_no+=1
#         curr_time=time.time()
#         if (curr_time-start_time)/conv_factor>=S*t:
#             runtime_chk=1

#     return OptCost,iter_no


# def VNS(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm): #variable neighbourhood search

#     start_time=time.time()

#     tm=0
#     Ns=10 #neighbourhood size

#     totserv=0
#     totcost=0
#     prev_class=4
#     queue_complete=0
#     weather_state=0

#     Anneal_Seq=[0]*NoA
#     for i in range(NoA):
#         Anneal_Seq[i]=ArrTime_Sorted[i][1]

#     NewCost=Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)
#     OldCost=NewCost
#     OptCost=NewCost
#     Opt_Seq=Anneal_Seq[:]

#     iter_no=0
#     runtime_chk=0
#     AC_remaining=NoA

#     c=0

#     while runtime_chk==0:

#         if c>=100:

#             #Perturb the optimal sequence
#             Anneal_Seq=Opt_Seq[:]

#             perm_size=min(4,AC_remaining) #no. of ACs to shuffle around
#             no_start_pos=AC_remaining-perm_size+1 #no. of possible start positions

#             triangle_dist_size=no_start_pos*(no_start_pos+1)/2

#             z=random.random()*triangle_dist_size
#             totp=0
#             for ii in range(no_start_pos):
#                 if z<(totp+no_start_pos-ii):
#                     pos=ii
#                     break
#                 totp+=no_start_pos-ii

#             remove_perm=[]
#             for ii in range(perm_size):
#                 AC=Anneal_Seq[pos]
#                 Anneal_Seq.remove(AC)
#                 remove_perm.append(AC)

#             old_perm=remove_perm[:]
#             random.shuffle(remove_perm)

#             for ii in range(perm_size):
#                 AC=remove_perm[ii]
#                 Anneal_Seq.insert(pos,AC)

#             # if len(Anneal_Seq)<NoA:
#             # 	print('** len(Anneal_Seq): '+str(len(Anneal_Seq)))

#             NewCost=Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)
#             if iter_no==0 or NewCost<OptCost:
#                 OptCost=NewCost
#                 Opt_Seq=Anneal_Seq[:]
#             OldCost=NewCost

#             c=0

#         #Move one plane in the sequence

#         N_OptCost=0
#         N_OptSeq=Anneal_Seq[:]

#         for i in range(Ns):

#             no_ACs=AC_remaining
#             triangle_dist_size=no_ACs*(no_ACs+1)/2

#             z=random.random()*triangle_dist_size
#             totp=0
#             for ii in range(no_ACs):
#                 if z<(totp+no_ACs-ii):
#                     pos=ii
#                     break
#                 totp+=no_ACs-ii

#             AC=Anneal_Seq[pos]
#             Anneal_Seq.remove(AC)
#             z1=int(random.random()*3)+1 #no. of places to move
#             z2=random.random() #determine whether to move up or down
#             if z2<0.5:
#                 z1=min(z1,pos)
#                 pos2=pos-z1
#             else:
#                 z1=min(z1,AC_remaining-1-pos)
#                 pos2=pos+z1

#             Anneal_Seq.insert(pos2,AC)

#             #End of mutation part

#             NewCost=Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)

#             if i==0 or NewCost<N_OptCost:
#                 N_OptCost=NewCost
#                 N_OptSeq=Anneal_Seq[:]

#             AC=Anneal_Seq[pos2]
#             Anneal_Seq.remove(AC)
#             Anneal_Seq.insert(pos,AC)

#         if N_OptCost<OldCost:
#             if iter_no==0 or N_OptCost<OptCost:
#                 OptCost=N_OptCost
#                 Opt_Seq=N_OptSeq[:]
#             OldCost=N_OptCost
#             Anneal_Seq=N_OptSeq[:]
#         else:
#             # AC=Anneal_Seq[pos2]
#             # Anneal_Seq.remove(AC)
#             # Anneal_Seq.insert(pos,AC)
#             c+=1

#         iter_no+=1
#         curr_time=time.time()
#         if (curr_time-start_time)/conv_factor>=S*t:
#             runtime_chk=1

#     return OptCost,iter_no

# def ILS(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm):

#     start_time=time.time()

#     tm=0

#     totserv=0
#     totcost=0
#     prev_class=4
#     queue_complete=0
#     weather_state=0

#     Anneal_Seq=[0]*NoA
#     for i in range(NoA):
#         Anneal_Seq[i]=ArrTime_Sorted[i][1]

#     NewCost=Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)
#     OldCost=NewCost
#     OptCost=NewCost
#     Opt_Seq=Anneal_Seq[:]

#     iter_no=0
#     runtime_chk=0
#     AC_remaining=NoA

#     c=0

#     while runtime_chk==0:

#         if c>=100:

#             #Perturb the optimal sequence
#             Anneal_Seq=Opt_Seq[:]

#             perm_size=min(4,AC_remaining) #no. of ACs to shuffle around
#             no_start_pos=AC_remaining-perm_size+1 #no. of possible start positions

#             triangle_dist_size=no_start_pos*(no_start_pos+1)/2

#             z=random.random()*triangle_dist_size
#             totp=0
#             for ii in range(no_start_pos):
#                 if z<(totp+no_start_pos-ii):
#                     pos=ii
#                     break
#                 totp+=no_start_pos-ii

#             remove_perm=[]
#             for ii in range(perm_size):
#                 AC=Anneal_Seq[pos]
#                 Anneal_Seq.remove(AC)
#                 remove_perm.append(AC)

#             old_perm=remove_perm[:]
#             random.shuffle(remove_perm)

#             for ii in range(perm_size):
#                 AC=remove_perm[ii]
#                 Anneal_Seq.insert(pos,AC)

#             # if len(Anneal_Seq)<NoA:
#             # 	print('** len(Anneal_Seq): '+str(len(Anneal_Seq)))

#             NewCost=Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)
#             if iter_no==0 or NewCost<OptCost:
#                 OptCost=NewCost
#                 Opt_Seq=Anneal_Seq[:]
#             OldCost=NewCost

#             c=0

#         #Move one plane in the sequence
#         no_ACs=AC_remaining
#         triangle_dist_size=no_ACs*(no_ACs+1)/2

#         z=random.random()*triangle_dist_size
#         totp=0
#         for ii in range(no_ACs):
#             if z<(totp+no_ACs-ii):
#                 pos=ii
#                 break
#             totp+=no_ACs-ii

#         AC=Anneal_Seq[pos]
#         Anneal_Seq.remove(AC)
#         z1=int(random.random()*3)+1 #no. of places to move
#         z2=random.random() #determine whether to move up or down
#         if z2<0.5:
#             z1=min(z1,pos)
#             pos2=pos-z1
#         else:
#             z1=min(z1,AC_remaining-1-pos)
#             pos2=pos+z1

#         Anneal_Seq.insert(pos2,AC)

#         #End of mutation part

#         NewCost=Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)

#         if NewCost<OldCost:
#             if iter_no==0 or NewCost<OptCost:
#                 OptCost=NewCost
#                 Opt_Seq=Anneal_Seq[:]
#             OldCost=NewCost
#         else:
#             AC=Anneal_Seq[pos2]
#             Anneal_Seq.remove(AC)
#             Anneal_Seq.insert(pos,AC)
#             c+=1

#         iter_no+=1
#         curr_time=time.time()
#         if (curr_time-start_time)/conv_factor>=S*t:
#             runtime_chk=1

#     return OptCost,iter_no

# def Tabu(Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm):

#     start_time=time.time()

#     tm=0

#     totserv=0
#     totcost=0
#     prev_class=4
#     queue_complete=0
#     weather_state=0

#     Anneal_Seq=[0]*NoA
#     for i in range(NoA):
#         Anneal_Seq[i]=i #ArrTime_Sorted[i][1]

#     NewCost=Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)
#     OldCost=NewCost
#     OptCost=NewCost
#     Opt_Seq=Anneal_Seq[:]

#     Tabu_List=[]
#     Tabu_Size=25 #50
#     Tabu_Moves=10 #no. of moves to consider in one step

#     iter_no=0
#     runtime_chk=0
#     AC_remaining=NoA

#     while runtime_chk==0:

#         move_no=0
#         BestCost=-1
#         Old_Seq=Anneal_Seq[:]
#         Best_Seq=Old_Seq[:]

#         no_ACs=AC_remaining

#         while move_no<Tabu_Moves:

#             #Move one plane in the sequence
#             triangle_dist_size=no_ACs*(no_ACs+1)/2

#             z=random.random()*triangle_dist_size
#             totp=0
#             for ii in range(no_ACs):
#                 if z<(totp+no_ACs-ii):
#                     pos=ii
#                     break
#                 totp+=no_ACs-ii

#             AC=Anneal_Seq[pos]
#             Anneal_Seq.remove(AC)
#             z1=int(random.random()*3)+1 #no. of places to move
#             z2=random.random() #determine whether to move up or down
#             if z2<0.5:
#                 z1=min(z1,pos)
#                 pos2=pos-z1
#             else:
#                 z1=min(z1,AC_remaining-1-pos)
#                 pos2=pos+z1

#             Anneal_Seq.insert(pos2,AC)

#             if Anneal_Seq not in Tabu_List:
#                 NewCost=Annealing_Cost(Anneal_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)
#                 if BestCost==-1 or NewCost<BestCost:
#                     BestCost=NewCost
#                     Best_Seq=Anneal_Seq[:]

#             Anneal_Seq=Old_Seq[:]
#             move_no+=1

#         if BestCost>=0 and BestCost<OptCost:
#             OptCost=BestCost
#             Opt_Seq=Best_Seq[:]

#         Tabu_List.append(Best_Seq)
#         #print('Added sequence with cost '+str(BestCost)+' to tabu list, no. of iters is '+str(iter_no))

#         if len(Tabu_List)>Tabu_Size:
#             Tabu_List.remove(Tabu_List[0])
#             #print('Removed sequence with cost ??? from tabu list')

#         Anneal_Seq=Best_Seq[:]

#         iter_no+=1
#         curr_time=time.time()
#         if (curr_time-start_time)/conv_factor>=S*t:
#             runtime_chk=1

#     NewCost=Annealing_Cost(Opt_Seq,Ac_Info,ArrTime,ServTime,ArrTime_Sorted,wlb_tm,wub_tm,0)

#     queue_complete=0
#     prev_class=4
#     reltime=0
#     for i in range(NoA):
#         AC=Opt_Seq[i]
#         Ac_Infoi=Ac_Info[AC]
#         current_class=Ac_Infoi[1]
#         reltime=max(reltime,Ac_Infoi[9])
#         begin_serv=max(reltime,queue_complete)
#         t1=reltime+Ac_Infoi[6]
#         weather_state=weather(reltime,wlb_tm,wub_tm) #weather(begin_serv,wlb_tm,wub_tm)
#         t2=begin_serv+Get_Actual_Serv(AC,prev_class,current_class,weather_state,k,Time_Sep)
#         finish_time=max(t1,t2)

#         f.write('Tabu'+','+str(rep)+','+str(Ac_Infoi[19])+','+str(AC)+','+str(prev_class)+','+str(current_class)+','+str(Time_Sep[prev_class][current_class]/60)+','+str(Ac_Infoi[2])+','+str(Ac_Infoi[9])+','+str(reltime)+','+str(Ac_Infoi[6])+','+str(weather_state)+','+str(begin_serv)+','+str(Get_Actual_Serv(AC,prev_class,current_class,weather_state,k,Time_Sep))+','+str(finish_time)+','+str(max(0,finish_time-(Ac_Infoi[2]+thres1)))+','+str(finish_time-(Ac_Infoi[9]+Ac_Infoi[6]))+','+str(Ac_Infoi[10])+','+str(getcost(Ac_Infoi[2],Ac_Infoi[9],Ac_Infoi[6],finish_time,Ac_Infoi[10],thres1, thres2, lam1, lam2))+',')
#         f.write(str(Ac_Infoi[13])+','+str(Ac_Infoi[14])+',')

#         queue_complete=finish_time
#         prev_class=current_class

#         f.write('\n')

#     return OptCost,iter_no

# def FCFS_rule(Ac_Info,Arr_Pool,Dep_Pool): #Not using this anymore

#     start_time=time.time()

#     Ac_added=[]

#     while len(Arr_Pool)+len(Dep_Pool)>len(Ac_added):

#         if len(Arr_Pool)>0 and len(Dep_Pool)>0:
#             arr_ac=Arr_Pool[0]
#             dep_ac=Dep_Pool[0]
#             if Ac_Info[arr_ac][2]<=Ac_Info[dep_ac][2]:
#                 Ac_added.append(arr_ac)
#                 #print('* Added aircraft '+str(arr_ac)+' arrival to the queue.')
#             else:
#                 Ac_added.append(dep_ac)
#                 #print('* Added aircraft '+str(dep_ac)+' departure to the queue.')
#         elif len(Arr_Pool)>0:
#             arr_ac=Arr_Pool[0]
#             Ac_added.append(arr_ac)
#             #print('* Added aircraft '+str(arr_ac)+' arrival to the queue.')
#         else:
#             dep_ac=Dep_Pool[0]
#             Ac_added.append(dep_ac)
#             #print('* Added aircraft '+str(dep_ac)+' departure to the queue.')

#     end_time=time.time()
#     elap=(end_time-start_time)/conv_factor #convert into minutes
#     if elap<0.001:
#         elap=0.001

#     #elap=fixed_elap/conv_factor

#     return Ac_added,elap