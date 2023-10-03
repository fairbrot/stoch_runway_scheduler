from typing import List
import time
from .utils import weather, getcost, truncexp
from .gamma import gamma_cond_exp

def Genetic_determ(Ac_Info,Arr_Pool,Arr_NotReady,Ac_queue,Left_queue,tm,NoA,k,prev_class,GA_PopList,GA_Info,wiener_cdf,wlb,wub,Opt_List, norm_approx_min: float, tau: int, Max_LookAhead: int, Time_Sep: List[List[int]], thres1:int, thres2: int, lam1: float, lam2: float, tot_arr_cost: float, tot_dep_cost: float, w_rho: float, stepthrough:int, step_summ:int, step_new: int):

	output=0 #output==1 means we're printing results as we go along; output==2 means we're outputting results to "Detailed" csv file
	ee=0
	if stepthrough==1:
		ee=1
	start_time=time.time()

	if k>=norm_approx_min:
		NormalApprox=1
	else:
		NormalApprox=0

	if stepthrough==1:
		st.write('Now entering Genetic procedure'+'\n')

	ArrTime=[0]*NoA
	ServTime=[0]*NoA
	Trav_Time=[0]*NoA

	#Generate arrival and service time percentiles for AC not yet in queue

	ArrTime_Sorted=[]

	for AC in Arr_Pool:
		ArrTime[AC]=Ac_Info[AC][9]
		ArrTime_Sorted.append([ArrTime[AC],AC])

	for AC in Arr_NotReady:
		ArrTime[AC]=max(0,Ac_Info[AC][3]-tau)
		ArrTime_Sorted.append([ArrTime[AC],AC])

	ArrTime_Sorted.sort(key=lambda x: x[0])

	basecost=tot_arr_cost+tot_dep_cost
	if stepthrough==1:
		st.write('basecost is '+','+str(basecost)+'\n'+'\n')
		st.write('Generated results for ACs already in the queue are as follows:'+'\n')
		st.write('AC'+','+'Class'+','+'Time Sep'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

	if len(Ac_queue)>0:

		#Need to generate service times for AC already in the queue; first consider the customer in position 0

		if NormalApprox==0:

			AC=Ac_queue[0]
			Ac_Infoi=Ac_Info[AC]
			rel_time=Ac_Infoi[4]
			sv_time=Ac_Infoi[5]
			cur_class=Ac_Infoi[1]

		# 	#print('AC: '+str(AC)+' Ac_Info: '+str(Ac_Info[AC])+' tm: '+str(tm)+' sv_time: '+str(sv_time)+' prev_class: '+str(prev_class)+' cur_class: '+str(cur_class))
		# 	if Ac_Infoi[12]==1 and (tm-sv_time)*10>=IFR_last_epoch[prev_class][cur_class]:
		# 		ph_B=k-1
		# 	elif (tm-sv_time)*10>=last_epoch[prev_class][cur_class]:
		# 		ph_B=k-1
		# 	else:
		# 		z2=0.5 #use the median
		# 		TotPr=0
		# 		chk_cond=0
		# 		j=0
		# 		while j<=k and chk_cond==0: #for j in range(k+1):
		# 			if Ac_Infoi[12]==1:
		# 				TotPr+=IFR_Serv_Pr[prev_class][cur_class][int((tm-sv_time)*10)][j]
		# 			else:
		# 				TotPr+=Serv_Pr[prev_class][cur_class][int((tm-sv_time)*10)][j]
		# 			if z2<TotPr:
		# 				ph_B=j
		# 				chk_cond=1
		# 			j+=1

		# 	assert ph_B<k

		# 	if stepthrough==1:
		# 		st.write(str(AC)+','+str(Ac_Infoi[1])+','+str(Time_Sep[prev_class][cur_class]/60)+','+str(Ac_Infoi[4])+','+str(Ac_Infoi[6])+','+str(Ac_Infoi[5])+',')

			t1=Ac_Infoi[3]
		# 	if Ac_Infoi[12]==1:
		# 		t2=tm+(w_rho*Time_Sep[prev_class][cur_class]/60)*((k-ph_B)/k)
		# 		#print('Median value ph_B: '+str(ph_B)+' tm sep: '+str((w_rho*Time_Sep[prev_class][cur_class]/60))+' time remaining: '+str((w_rho*Time_Sep[prev_class][cur_class]/60)*((k-ph_B)/k)))
		# 		if stepthrough==1:
		# 			st.write(str((w_rho*Time_Sep[prev_class][cur_class]/60)*((k-ph_B)/k)))
		# 	else:
		# 		t2=tm+(Time_Sep[prev_class][cur_class]/60)*((k-ph_B)/k)
		# 		#print('Median value ph_B: '+str(ph_B)+' tm sep: '+str((Time_Sep[prev_class][cur_class]/60))+' time remaining: '+str((Time_Sep[prev_class][cur_class]/60)*((k-ph_B)/k)))
		# 		if stepthrough==1:
		# 			st.write(str((Time_Sep[prev_class][cur_class]/60)*((k-ph_B)/k)))

			#Get the conditional expectation of service time based on service time elapsed so far

			if Ac_Infoi[12]==1:
				beta=k/(w_rho*Time_Sep[prev_class][cur_class]/60)
			else:
				beta=k/(Time_Sep[prev_class][cur_class]/60)
			alpha=k
			cond_tm=gamma_cond_exp(tm-sv_time,alpha,beta)

			t2=sv_time+cond_tm #tm+(cond_tm-sv_time)
			#print('Time remaining based on mean value: '+str(sv_time+cond_tm-tm))

			queue_complete=max(t1,t2)

			basecost+=getcost(Ac_Infoi[18],Ac_Infoi[9],tau,queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2)
			if stepthrough==1:
				st.write(str(queue_complete)+','+str(Ac_Infoi[10])+','+str(getcost(Ac_Infoi[2],Ac_Infoi[9],tau,queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2))+'\n')

			perm_prev_class=cur_class\

		else:

			AC=Ac_queue[0]
			Ac_Infoi=Ac_Info[AC]
			t1=Ac_Infoi[3]
			sv_time=Ac_Infoi[5]
			cur_class=Ac_Infoi[1]

			Mn=Time_Sep[prev_class][cur_class]/60
			SD=(Mn**2)/k
			t2=truncexp(sv_time+Mn,SD,tm)

			queue_complete=max(t1,t2)

			basecost+=getcost(Ac_Infoi[18],Ac_Infoi[9],tau,queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2)
			AC=Ac_queue[0]
			perm_prev_class=cur_class

		#Now consider the rest of the customers in the queue
		for j in range(1,len(Ac_queue)):

			AC=Ac_queue[j]
			Ac_Infoi=Ac_Info[AC]
			rel_time=Ac_Infoi[4]
			cur_class=Ac_Infoi[1]

			if stepthrough==1:
				st.write(str(AC)+','+str(Ac_Infoi[1])+','+str(Time_Sep[perm_prev_class][cur_class]/60)+','+str(Ac_Infoi[4])+','+str(Ac_Infoi[6])+','+str(Ac_Infoi[5])+',')

			t1=Ac_Infoi[3]
			weather_state=weather(rel_time,wlb,wub) #weather(queue_complete,wlb,wub)
			if weather_state==1: #Ac_Infoi[12]==1:
				t2=queue_complete+(w_rho*Time_Sep[perm_prev_class][cur_class]/60)
				if stepthrough==1:
					st.write(str(w_rho*Time_Sep[perm_prev_class][cur_class]/60)+',')
			else:
				t2=queue_complete+(Time_Sep[perm_prev_class][cur_class]/60)
				if stepthrough==1:
					st.write(str(Time_Sep[perm_prev_class][cur_class]/60)+',')
			queue_complete=max(t1,t2)

			perm_prev_class=cur_class
			basecost+=getcost(Ac_Infoi[18],Ac_Infoi[9],tau,queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2)

			if stepthrough==1:
				st.write(str(queue_complete)+','+str(Ac_Infoi[10])+','+str(getcost(Ac_Infoi[2],Ac_Infoi[9],tau,queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2))+'\n')

	else:

		queue_complete=tm
		perm_prev_class=prev_class

	stored_prev_class=perm_prev_class
	stored_queue_complete=queue_complete

	#Try all the sequences in the population

	for j in range(len(GA_Info)):

		if stepthrough==1:
			st.write('\n'+'Now trying sequence '+','+str(GA_Info[j][0])+'\n')
			st.write('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

		permcost=basecost
		latest_tm=tm
		perm_prev_class=stored_prev_class
		perm_queue_complete=queue_complete
		#perm_weather_state=weather_state

		perm=GA_Info[j][0]
		GA_Info[j][1]+=1
		#gam=0.01
		gam=1/GA_Info[j][1]

		#print('GA_Info: '+str(GA_Info))

		no_ACs=min(Max_LookAhead,len(perm)) #min(AC_List_Length,len(Arr_Pool)+len(Arr_NotReady))
		for index in range(no_ACs):

			#index=perm[i]
			AC=perm[index]
			Ac_Infoi=Ac_Info[AC]
			perm_class=Ac_Infoi[1]
			reltime=max(latest_tm,ArrTime[AC])

			if stepthrough==1:
				st.write(str(AC)+','+str(Ac_Infoi[1])+','+str(Time_Sep[perm_prev_class][perm_class]/60)+','+str(ArrTime[AC])+','+str(reltime)+','+str(tau)+','+str(perm_queue_complete)+',')

			t1=reltime+tau
			if reltime>=wlb and reltime<=wub:
				exp_serv=w_rho*Time_Sep[perm_prev_class][perm_class]/60
			else:
				exp_serv=Time_Sep[perm_prev_class][perm_class]/60
			t2=perm_queue_complete+exp_serv
			if stepthrough==1:
				st.write(str(exp_serv)+',')

			AC_FinishTime=max(t1,t2)
			if t1>=t2:
				straight_into_service=1
			else:
				straight_into_service=0

			GA_Info[j][3][index]=(1-gam)*GA_Info[j][3][index]+gam*straight_into_service

			permcost+=getcost(Ac_Infoi[18],ArrTime[AC],tau,AC_FinishTime,Ac_Infoi[10],thres1,thres2, lam1, lam2) #Ac_Infoi[10]*(AC_FinishTime-(Ac_Infoi[2]+thres))**2
			latest_tm=reltime

			if stepthrough==1:
				st.write(str(AC_FinishTime)+','+str(Ac_Infoi[10])+','+str(getcost(Ac_Infoi[2],ArrTime[AC],tau,AC_FinishTime,Ac_Infoi[10],thres1,thres2, lam1, lam2))+'\n')

			perm_queue_complete=AC_FinishTime
			perm_prev_class=perm_class

		if stepthrough==1:
			st.write('Final cost: '+','+str(permcost)+','+'Gamma: '+','+str(gam)+','+'Old cost: '+','+str(GA_Info[j][2])+',')

		GA_Info[j][2]=(1-gam)*GA_Info[j][2]+gam*permcost
		if stepthrough==1:
			st.write('Total cost: '+','+str(GA_Info[j][2])+','+'Queue probs: '+','+str(GA_Info[j][3])+'\n'+'\n')

	if step_summ==1:
		GA_Info.sort(key=lambda x: x[0])
		for j in range(len(GA_Info)):
			st2.write(str(GA_Info[j][2])+',')
		st2.write('\n')

	#### UPDATE SEQS IN THE OPT LIST ####

	for j in range(len(Opt_List)):

		if stepthrough==1:
			st.write('\n'+'Now trying sequence '+','+str(Opt_List[j][0])+'\n')
			st.write('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

		permcost=basecost
		latest_tm=tm
		perm_prev_class=stored_prev_class
		perm_queue_complete=queue_complete
		#perm_weather_state=weather_state

		perm=Opt_List[j][0]
		Opt_List[j][1]+=1
		#gam=0.01
		gam=1/Opt_List[j][1]

		#print('GA_Info: '+str(GA_Info))

		no_ACs=min(Max_LookAhead,len(perm)) #min(AC_List_Length,len(Arr_Pool)+len(Arr_NotReady))
		for index in range(no_ACs):

			#index=perm[i]
			AC=perm[index]
			Ac_Infoi=Ac_Info[AC]
			perm_class=Ac_Infoi[1]
			reltime=max(latest_tm,ArrTime[AC])

			if stepthrough==1:
				st.write(str(AC)+','+str(Ac_Infoi[1])+','+str(Time_Sep[perm_prev_class][perm_class]/60)+','+str(ArrTime[AC])+','+str(reltime)+','+str(tau)+','+str(perm_queue_complete)+',')

			t1=reltime+tau
			if reltime>=wlb and reltime<=wub:
				exp_serv=w_rho*Time_Sep[perm_prev_class][perm_class]/60
			else:
				exp_serv=Time_Sep[perm_prev_class][perm_class]/60
			t2=perm_queue_complete+exp_serv
			if stepthrough==1:
				st.write(str(exp_serv)+',')

			AC_FinishTime=max(t1,t2)
			if t1>=t2:
				straight_into_service=1
			else:
				straight_into_service=0

			Opt_List[j][3][index]=(1-gam)*Opt_List[j][3][index]+gam*straight_into_service

			permcost+=getcost(Ac_Infoi[18],ArrTime[AC],tau,AC_FinishTime,Ac_Infoi[10],thres1,thres2, lam1, lam2) #Ac_Infoi[10]*(AC_FinishTime-(Ac_Infoi[2]+thres))**2
			latest_tm=reltime

			if stepthrough==1:
				st.write(str(AC_FinishTime)+','+str(Ac_Infoi[10])+','+str(getcost(Ac_Infoi[2],ArrTime[AC],tau,AC_FinishTime,Ac_Infoi[10],thres1,thres2, lam1, lam2))+'\n')

			perm_queue_complete=AC_FinishTime
			perm_prev_class=perm_class

		if stepthrough==1:
			st.write('Final cost: '+','+str(permcost)+','+'Gamma: '+','+str(gam)+','+'Old cost: '+','+str(Opt_List[j][2])+',')

		Opt_List[j][2]=(1-gam)*Opt_List[j][2]+gam*permcost
		if stepthrough==1:
			st.write('Total cost: '+','+str(Opt_List[j][2])+','+'Queue probs: '+','+str(Opt_List[j][3])+'\n'+'\n')

	if step_summ==1:
		Opt_List.sort(key=lambda x: x[0])
		for j in range(len(Opt_List)):
			st2.write(str(Opt_List[j][2])+',')
		st2.write('\n')

	# minperm=0
	# mincost=0
	# for j in range(len(GA_Info)):
	# 	if j==0 or GA_Info[j][2]<mincost:
	# 		mincost=GA_Info[j][2]
	# 		minperm=j

	Opt_List.sort(key=lambda x: x[2])
	GA_Info.sort(key=lambda x: x[2])

	Ac_added=[]
	counter=0
	qp=0

	if len(Arr_Pool)>0:

		if len(GA_Info)>0 and len(Opt_List)>0:
			if Opt_List[0][2]<GA_Info[0][2]:
				perm=Opt_List[0]
			else:
				perm=GA_Info[0]
		elif len(Opt_List)>0:
			perm=Opt_List[0]
		elif len(GA_Info)>0:
			perm=GA_Info[0]
		else:
			assert 1==2

		if perm[0][0] in Arr_Pool:

			if 1==1: #perm[1]>=100 and 1==1: #perm[3][0]>0: #0.05:
				j=0
				counter=perm[1]
				qp=perm[3][0]
				while 1==1: #perm[3][j]>0: #0.05:
					AC=perm[0][j]
					Ac_added.append(AC)
					j+=1
					if j==len(perm[0]):
						break

	# end_time=time.time()
	# elap=(end_time-start_time)/conv_factor
	# #elap=0.01

	#elap=fixed_elap_vnsd/conv_factor

	return Ac_added,counter,qp,stored_queue_complete

