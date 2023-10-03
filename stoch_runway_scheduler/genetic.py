from typing import List
import random
import math
import time
from .utils import weather, getcost, Normal_GetServ, Normal_GetServ_Future, Normal_Conditional_GetServ
from .gamma import Gamma_GetServ, Gamma_GetServ_Future, Gamma_Conditional_GetServ

def Genetic(Ac_Info,Arr_Pool,Arr_NotReady,Ac_queue,Left_queue,tm,NoA,k,prev_class,GA_PopList,GA_Info,wiener_cdf,GA_LoopSize,GA_CheckSize,GA_counter,basecost,wlb,wub,weather_cdf,Opt_List,max_d,soln_evals_tot,soln_evals_num,gamma_cdf,normcdf, norm_approx_min: float, tau: int, Max_LookAhead: int, Time_Sep: List[List[int]], thres1: int, thres2: int, lam1: float, lam2: float, GA_Check_Increment: int, Opt_Size: int, w_rho: float, stepthrough: int, step_summ: int, step_new: int):

	# if len(Arr_Pool)+len(Ac_queue)>0:
	# 	print('Genetic')

	output=0 #output=1 means we're printing results as we go along; output=2 means we're outputting results to "Detailed" csv file
	ee=0
	pruned=0 #indicator of whether or not the number of sequences has gone below the minimum number

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

	wiener_cdf_tau=wiener_cdf[int(10*tau)]

	#Generate arrival and service time percentiles for AC not yet in queue

	for AC in Arr_Pool:

		ArrTime[AC]=max(0,Ac_Info[AC][3]-tau)

		z=int(random.randrange(1,999))
		Trav_Time[AC]=wiener_cdf_tau[z]

		if NormalApprox==0: 
			z=int(random.random()*1000)
			ServTime[AC]=gamma_cdf[z]
			# ServTime[AC]=[0]*k
			# for j in range(k):
			# 	ServTime[AC][j]=random.random()
		else:
			z=int(random.random()*10000)
			ServTime[AC]=normcdf[z]

	for AC in Arr_NotReady:

		z=int(random.randrange(1,999))
		sched=int(10*round(Ac_Info[AC][3]-(tm+tau),1))
		#print('tm: '+str(tm)+' AC: '+str(AC)+' sched: '+str(sched))
		ArrTime[AC]=wiener_cdf[sched][z]+tm
		#st4.write(str(tm)+','+str(Ac_Info[AC][3]-(tm+tau))+','+str(ArrTime[AC])+',')

		z=int(random.randrange(1,999))
		Trav_Time[AC]=wiener_cdf_tau[z]
		#st4.write(str(Trav_Time[AC])+'\n')

		#assert 1==2 

		if NormalApprox==0:
			z=int(random.random()*1000)
			ServTime[AC]=gamma_cdf[z]
			# ServTime[AC]=[0]*k
			# for j in range(k):
			# 	ServTime[AC][j]=random.random()
		else:
			z=int(random.random()*10000)
			ServTime[AC]=normcdf[z]

		# if NormalApprox==0:
		# 	ServTime[AC]=[0]*k
		# 	for j in range(k):
		# 		ServTime[AC][j]=random.random()

	#Before proceeding, randomly generate wlb_gen and wub_gen
	chk=0
	while chk==0:
		if tm>=wlb:
			wlb_gen=wlb
		else:
			#Do wlb_gen
			z=int(random.randrange(1,999))
			sched=int(10*round(wlb-tm,1))
			wlb_gen=weather_cdf[sched][z]+tm
		if tm>=wub:
			wub_gen=wub
		else:
			#Do wub_gen
			z=int(random.randrange(1,999))
			sched=int(10*round(wub-tm,1))
			wub_gen=weather_cdf[sched][z]+tm
		if wlb_gen<=wub_gen:
			chk=1
		else:
			chk=1
			#print('tm: '+str(tm)+' wlb: '+str(wlb)+' wub: '+str(wub)+' wlb_gen: '+str(wlb_gen)+' wub_gen: '+str(wub_gen))

	if stepthrough==1:
		st.write('basecost is '+','+str(basecost)+'\n'+'\n')
		st.write('Generated results for ACs already in the queue are as follows:'+'\n')
		st.write('AC'+','+'Class'+','+'Time Sep'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

	if len(Ac_queue)>0:

		#Need to generate service times for AC already in the queue; first consider the customer in position 0

		if NormalApprox==0:

			# AC=Ac_queue[0]
			# Ac_Infoi=Ac_Info[AC]
			# rel_time=Ac_Infoi[4]
			# sv_time=Ac_Infoi[5]
			# cur_class=Ac_Infoi[1]
			# weather_state=Ac_Infoi[12]

			# if tm>=Ac_Infoi[3]:
			# 	trav_time=Ac_Infoi[6] #travel time has already finished
			# else:
			# 	z=int(random.randrange(1,999))
			# 	sched=int(10*round(Ac_Infoi[3]-tm,1))
			# 	trav_time=wiener_cdf[sched][z]

			# # #print('AC: '+str(AC)+' Ac_Info: '+str(Ac_Info[AC])+' tm: '+str(tm)+' sv_time: '+str(sv_time)+' prev_class: '+str(prev_class)+' cur_class: '+str(cur_class))
			# # if weather_state==1 and (tm-sv_time)*10>=IFR_last_epoch[prev_class][cur_class]:
			# # 	ph_B=k-1
			# # elif (tm-sv_time)*10>=last_epoch[prev_class][cur_class]:
			# # 	ph_B=k-1
			# # else:
			# # 	z2=random.random()
			# # 	TotPr=0
			# # 	chk_cond=0
			# # 	j=0
			# # 	if weather_state==1:
			# # 		Serv_Pr_Row=IFR_Serv_Pr[prev_class][cur_class][int((tm-sv_time)*10)]
			# # 	else:
			# # 		Serv_Pr_Row=Serv_Pr[prev_class][cur_class][int((tm-sv_time)*10)]
			# # 	while j<=k and chk_cond==0: #for j in range(k+1):
			# # 		if weather_state==1:
			# # 			TotPr+=Serv_Pr_Row[j]
			# # 		else:
			# # 			TotPr+=Serv_Pr_Row[j]
			# # 		if z2<TotPr:
			# # 			ph_B=j #ph_B is randomly-sampled number of service phases completed (or remaining??) for AC currently in service
			# # 			chk_cond=1
			# # 		j+=1

			# # assert ph_B<k
			# # serv_time=[0]*(k-ph_B) #k-ph_B is number of phases remaining (??)
			# # for m in range(k-ph_B):
			# # 	serv_time[m]=random.random()

			# if stepthrough==1:
			# 	st.write(str(AC)+','+str(Ac_Infoi[1])+','+str(Time_Sep[prev_class][cur_class]/60)+','+str(Ac_Infoi[4])+','+str(Ac_Infoi[6])+','+str(Ac_Infoi[5])+',')

			# queue_complete,straight_into_service=GetServTime(trav_time,rel_time,prev_class,cur_class,tm,sv_time,ee,weather_state,gamma_cdf)
			# basecost+=getcost(Ac_Infoi[18],Ac_Infoi[9],trav_time,queue_complete,Ac_Infoi[10],thres1,thres2)

			# if stepthrough==1:
			# 	st.write(str(queue_complete)+','+str(Ac_Infoi[10])+','+str(getcost(Ac_Infoi[2],Ac_Infoi[9],trav_time,queue_complete,Ac_Infoi[10],thres1,thres2))+'\n')

			# perm_prev_class=cur_class

			####

			AC=Ac_queue[0]
			Ac_Infoi=Ac_Info[AC]
			rel_time=Ac_Infoi[4]
			sv_time=Ac_Infoi[5]
			cur_class=Ac_Infoi[1]
			weather_state=Ac_Infoi[12]

			if tm>=Ac_Infoi[3]:
				trav_time=Ac_Infoi[6] #travel time has already finished
			else:
				z=int(random.randrange(1,999))
				sched=int(10*round(Ac_Infoi[3]-tm,1))
				trav_time=wiener_cdf[sched][z]

			queue_complete,straight_into_service=Gamma_Conditional_GetServ(k, Time_Sep, trav_time,rel_time,sv_time,prev_class,cur_class,tm,weather_state,gamma_cdf, w_rho)
			basecost+=getcost(Ac_Infoi[18],Ac_Infoi[9],trav_time,queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2)
			perm_prev_class=cur_class

		else:

			AC=Ac_queue[0]
			Ac_Infoi=Ac_Info[AC]
			rel_time=Ac_Infoi[4]
			sv_time=Ac_Infoi[5]
			cur_class=Ac_Infoi[1]
			weather_state=Ac_Infoi[12]

			if tm>=Ac_Infoi[3]:
				trav_time=Ac_Infoi[6] #travel time has already finished
			else:
				z=int(random.randrange(1,999))
				sched=int(10*round(Ac_Infoi[3]-tm,1))
				trav_time=wiener_cdf[sched][z]

			queue_complete,straight_into_service=Normal_Conditional_GetServ(trav_time,rel_time,sv_time,prev_class,cur_class,tm,weather_state, Time_Sep, normcdf, w_rho, k)
			basecost+=getcost(Ac_Infoi[18],Ac_Infoi[9],trav_time,queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2)
			perm_prev_class=cur_class

		#Now consider the rest of the customers in the queue
		for j in range(1,len(Ac_queue)):

			AC=Ac_queue[j]
			Ac_Infoi=Ac_Info[AC]
			rel_time=Ac_Infoi[4]
			cur_class=Ac_Infoi[1]
			weather_state=weather(rel_time,wlb_gen,wub_gen) #weather(queue_complete,wlb_gen,wub_gen)

			z=int(random.randrange(1,999))
			trav_time=wiener_cdf_tau[z] #???????????????

			# if tm>=Ac_Infoi[3]: #this block of code is probably needed but wasn't included in the 5000 experiments for the paper
			# 	trav_time=Ac_Infoi[6] #travel time has already finished
			# else:
			# 	z=int(random.randrange(1,999))
			# 	sched=int(10*round(Ac_Infoi[3]-tm,1))
			# 	trav_time=wiener_cdf[sched][z]

			if NormalApprox==0:

				# # serv_time=[0]*k
				# # for m in range(k):
				# # 	serv_time[m]=random.random()

				# #sv_time=tm

				# if stepthrough==1:
				# 	st.write(str(AC)+','+str(Ac_Infoi[1])+','+str(Time_Sep[perm_prev_class][cur_class]/60)+','+str(Ac_Infoi[4])+','+str(Ac_Infoi[6])+','+str(Ac_Infoi[5])+',')

				# queue_complete,straight_into_service=GetServTime(trav_time,rel_time,perm_prev_class,cur_class,queue_complete,queue_complete,ee,weather_state,gamma_cdf)
				# perm_prev_class=cur_class
				# basecost+=getcost(Ac_Infoi[18],Ac_Infoi[9],trav_time,queue_complete,Ac_Infoi[10],thres1,thres2,lam1,lam2)


				queue_complete,straight_into_service=Gamma_GetServ(k, Time_Sep, rel_time,trav_time,perm_prev_class,cur_class,queue_complete,weather_state,gamma_cdf, w_rho)
				perm_prev_class=cur_class
				basecost+=getcost(Ac_Infoi[18],Ac_Infoi[9],trav_time,queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2)

			else:

				queue_complete,straight_into_service=Normal_GetServ(rel_time,trav_time,perm_prev_class,cur_class,queue_complete,weather_state, Time_Sep, normcdf, w_rho, k)
				perm_prev_class=cur_class
				basecost+=getcost(Ac_Infoi[18],Ac_Infoi[9],trav_time,queue_complete,Ac_Infoi[10],thres1,thres2, lam1, lam2)

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

		GA_Infoj=GA_Info[j]

		perm=GA_Infoj[0]
		GA_Infoj[1]+=1
		#gam=0.01
		#gam=2/(GA_Info[j][1]+1)
		gam=1/GA_Infoj[1]#
		#gam=0.1

		#print('GA_Info: '+str(GA_Info))

		no_ACs=min(Max_LookAhead,len(perm)) #min(AC_List_Length,len(Arr_Pool)+len(Arr_NotReady))
		for index in range(no_ACs):

			#index=perm[i]
			AC=perm[index]
			Ac_Infoi=Ac_Info[AC]
			perm_class=Ac_Infoi[1]
			reltime=max(latest_tm,ArrTime[AC])
			begin_serv=max(reltime,perm_queue_complete)
			weather_state=weather(reltime,wlb_gen,wub_gen) #weather(begin_serv,wlb_gen,wub_gen)

			if output==1:
				ee=1

			if stepthrough==1:
				st.write(str(AC)+','+str(perm_class)+','+str(Time_Sep[perm_prev_class][perm_class]/60)+','+str(ArrTime[AC])+','+str(reltime)+','+str(Trav_Time[AC])+','+str(perm_queue_complete)+',')

			if NormalApprox==0:
				#perm_weather_state=weather(reltime,wlb,wub)
				#AC_FinishTime,straight_into_service=GetServTime_Future(Trav_Time[AC],ServTime[AC],reltime,perm_prev_class,perm_class,perm_queue_complete,ee,weather_state)
				AC_FinishTime,straight_into_service=Gamma_GetServ_Future(k, Time_Sep, reltime,ServTime[AC],Trav_Time[AC],perm_prev_class,perm_class,perm_queue_complete,weather_state, w_rho)
			else:
				AC_FinishTime,straight_into_service=Normal_GetServ_Future(reltime,ServTime[AC],Trav_Time[AC],perm_prev_class,perm_class,perm_queue_complete,weather_state, Time_Sep, w_rho, k)
			GA_Infoj[3][index]=(1-gam)*GA_Infoj[3][index]+gam*straight_into_service

			permcost+=getcost(Ac_Infoi[18],ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi[10],thres1,thres2, lam1, lam2) #Ac_Infoi[10]*(AC_FinishTime-(Ac_Infoi[2]+thres))**2
			latest_tm=reltime

			if stepthrough==1:
				st.write(str(AC_FinishTime)+','+str(Ac_Infoi[10])+','+str(getcost(Ac_Infoi[2],ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi[10],thres1,thres2, lam1, lam2))+'\n')

			perm_queue_complete=AC_FinishTime
			perm_prev_class=perm_class

		if stepthrough==1:
			st.write('Final cost: '+','+str(permcost)+','+'Gamma: '+','+str(gam)+','+'Old cost: '+','+str(GA_Info[j][2])+',')

		# if GA_Infoj[0]==[0,1,3,6,8,7,10,12,2,4,9,11,5,13,14,15,20,16,17,21,22,24,28,23,25,26,27,29,18,19]:
		# 	print('Evaluated cost for sequence '+str(GA_Infoj[0])+' as '+str(permcost))

		# if GA_Infoj[2]>0:
		# 	pct_diff=abs(permcost/GA_Infoj[2]-1)
		# 	if pct_diff>iter_max_d:
		# 		iter_max_d=pct_diff
		# else:
		# 	iter_max_d=1

		GA_Infoj[2]=(1-gam)*GA_Infoj[2]+gam*permcost
		GA_Infoj[4]+=permcost**2
		if stepthrough==1:
			st.write('Total cost: '+','+str(GA_Info[j][2])+','+'Queue probs: '+','+str(GA_Info[j][3])+'\n'+'\n')

	if step_summ==1:
		GA_Info.sort(key=lambda x: x[0])
		for j in range(len(GA_Info)):
			st2.write(str(GA_Info[j][2])+',')
		st2.write('\n')

	######### NOW UPDATE OPT LIST SEQS ###################

	for j in range(len(Opt_List)):

		if stepthrough==1:
			st.write('\n'+'Now trying sequence '+','+str(Opt_List[j][0])+'\n')
			st.write('AC'+','+'Class'+','+'Time Sep'+','+'Arrives in pool'+','+'Release time'+','+'Travel time'+','+'Enters serv'+','+'Actual serv'+','+'Finish time'+','+'Pax weight'+','+'Cost'+'\n')

		permcost=basecost
		latest_tm=tm
		perm_prev_class=stored_prev_class
		perm_queue_complete=queue_complete
		#perm_weather_state=weather_state

		Opt_Listj=Opt_List[j]

		#print('Opt_Listj[0]: '+str(Opt_Listj[0]))
		perm=Opt_Listj[0]
		Opt_Listj[1]+=1
		#gam=0.01
		#gam=2/(GA_Info[j][1]+1)
		gam=1/Opt_Listj[1]
		#gam=0.1

		#print('GA_Info: '+str(GA_Info))

		#print('Opt_List: '+str(Opt_List))
		no_ACs=min(Max_LookAhead,len(perm)) #min(AC_List_Length,len(Arr_Pool)+len(Arr_NotReady))
		for index in range(no_ACs):

			#index=perm[i]
			AC=perm[index]
			Ac_Infoi=Ac_Info[AC]
			perm_class=Ac_Infoi[1]
			reltime=max(latest_tm,ArrTime[AC])
			begin_serv=max(reltime,perm_queue_complete)
			weather_state=weather(reltime,wlb_gen,wub_gen) #weather(begin_serv,wlb_gen,wub_gen)

			if output==1:
				ee=1

			if stepthrough==1:
				st.write(str(AC)+','+str(perm_class)+','+str(Time_Sep[perm_prev_class][perm_class]/60)+','+str(ArrTime[AC])+','+str(reltime)+','+str(Trav_Time[AC])+','+str(perm_queue_complete)+',')

			if NormalApprox==0:
				#perm_weather_state=weather(reltime,wlb,wub)
				#AC_FinishTime,straight_into_service=GetServTime_Future(Trav_Time[AC],ServTime[AC],reltime,perm_prev_class,perm_class,perm_queue_complete,ee,weather_state)
				AC_FinishTime,straight_into_service=Gamma_GetServ_Future(k, Time_Sep, reltime,ServTime[AC],Trav_Time[AC],perm_prev_class,perm_class,perm_queue_complete,weather_state, w_rho)
			else:
				AC_FinishTime,straight_into_service=Normal_GetServ_Future(reltime,ServTime[AC],Trav_Time[AC],perm_prev_class,perm_class,perm_queue_complete,weather_state, Time_Sep, w_rho, k)
			Opt_Listj[3][index]=(1-gam)*Opt_Listj[3][index]+gam*straight_into_service

			permcost+=getcost(Ac_Infoi[18],ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi[10],thres1,thres2, lam1, lam2) #Ac_Infoi[10]*(AC_FinishTime-(Ac_Infoi[2]+thres))**2
			latest_tm=reltime

			if stepthrough==1:
				st.write(str(AC_FinishTime)+','+str(Ac_Infoi[10])+','+str(getcost(Ac_Infoi[2],ArrTime[AC],Trav_Time[AC],AC_FinishTime,Ac_Infoi[10],thres1,thres2, lam1, lam2))+'\n')

			perm_queue_complete=AC_FinishTime
			perm_prev_class=perm_class

		if stepthrough==1:
			st.write('Final cost: '+','+str(permcost)+','+'Gamma: '+','+str(gam)+','+'Old cost: '+','+str(GA_Info[j][2])+',')

		# if Opt_Listj[2]>0:
		# 	pct_diff=abs(permcost/Opt_Listj[2]-1)
		# 	if pct_diff>iter_max_d:
		# 		iter_max_d=pct_diff
		# else:
		# 	iter_max_d=1

		Opt_Listj[2]=(1-gam)*Opt_Listj[2]+gam*permcost
		Opt_Listj[4]+=permcost**2 #Sum of squares

		if stepthrough==1:
			st.write('Total cost: '+','+str(Opt_Listj[2])+','+'Queue probs: '+','+str(Opt_Listj[3])+'\n'+'\n')

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

	GA_Info.sort(key=lambda x: x[2])
	Opt_List.sort(key=lambda x: x[2])

	GA_counter+=1

	if GA_counter>=GA_CheckSize:

		if step_new==1:
			st3.write('GA_counter: '+','+str(GA_counter)+'\n')
			st3.write('Arr_Pool: '+','+str(Arr_Pool)+'\n')
			st3.write('Ac_queue: '+','+str(Ac_queue)+'\n'+'\n')
			st3.write('GA_PopList:'+'\n')
			for i in range(len(GA_Info)):
				st3.write(str(GA_Info[i])+'\n')
			st3.write('\n')
			st3.write('Opt_List:'+'\n')
			for i in range(len(Opt_List)):
				st3.write(str(Opt_List[i])+'\n')
			st3.write('\n')

		# if totserv>315:
		# 	print('tm: '+str(tm)+' GA_counter: '+str(GA_counter)+' len(GA_Info): '+str(len(GA_Info))+' len(Opt_List): '+str(Opt_List))

		t_val=1.96 #97.5th percentile of normal dist

		for j in range(len(GA_Info)):

			GA_Infoj=GA_Info[j]

			if GA_Infoj[2]>0:

				mn1=GA_Infoj[2]
				n1=GA_Infoj[1]
				var1=(GA_Infoj[4]-(n1*mn1**2))/(n1-1)

				for m in range(len(GA_Info)):

					GA_Infom=GA_Info[m]

					if GA_Infom[2]>0:

						mn2=GA_Infom[2]
						n2=GA_Infom[1]
						var2=(GA_Infom[4]-(n2*mn2**2))/(n2-1)

						w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

						if mn1>mn2+w_val:

							# if step_new==1:
							# 	st3.write(str(GA_Infoj[0])+','+' loses to '+','+str(GA_Infom[0])+'\n')

							GA_Infoj[2]=-1
							break

						elif mn2>mn1+w_val:

							# if step_new==1:
							# 	st3.write(str(GA_Infom[0])+','+' loses to '+','+str(GA_Infoj[0])+'\n')

							GA_Infom[2]=-1

			if GA_Infoj[2]>0:

				for m in range(len(Opt_List)):

					Opt_Listm=Opt_List[m]

					if Opt_Listm[2]>0:

						mn2=Opt_Listm[2]
						n2=Opt_Listm[1]
						var2=(Opt_Listm[4]-(n2*mn2**2))/(n2-1)

						w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

						if mn1>mn2+w_val:

							# if step_new==1:
							# 	st3.write(str(GA_Infoj[0])+','+' loses to '+','+str(Opt_Listm[0])+'\n')

							GA_Infoj[2]=-1
							break

						elif mn2>mn1+w_val:

							# if step_new==1:
							# 	st3.write(str(Opt_Listm[0])+','+' loses to '+','+str(GA_Infoj[0])+'\n')

							Opt_Listm[2]=-1

		for j in range(len(Opt_List)):

			Opt_Listj=Opt_List[j]

			if Opt_Listj[2]>0:

				mn1=Opt_Listj[2]
				n1=Opt_Listj[1]
				var1=(Opt_Listj[4]-(n1*mn1**2))/(n1-1)

				for m in range(len(GA_Info)):

					GA_Infom=GA_Info[m]

					if GA_Infom[2]>0:

						mn2=GA_Infom[2]
						n2=GA_Infom[1]
						var2=(GA_Infom[4]-(n2*mn2**2))/(n2-1)

						w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

						if mn1>mn2+w_val:

							# if step_new==1:
							# 	st3.write(str(Opt_Listj[0])+','+' loses to '+','+str(GA_Infom[0])+'\n')

							Opt_Listj[2]=-1
							break

						elif mn2>mn1+w_val:

							# if step_new==1:
							# 	st3.write(str(GA_Infom[0])+','+' loses to '+','+str(Opt_Listj[0])+'\n')

							GA_Infom[2]=-1

			if Opt_Listj[2]>0:

				for m in range(len(Opt_List)):

					Opt_Listm=Opt_List[m]

					if Opt_Listm[2]>0:

						mn2=Opt_Listm[2]
						n2=Opt_Listm[1]
						var2=(Opt_Listm[4]-(n2*mn2**2))/(n2-1)

						w_val=math.sqrt(((t_val**2)*var1/n1)+((t_val**2)*var2/n2))

						if mn1>mn2+w_val:

							# if step_new==1:
							# 	st3.write(str(Opt_Listj[0])+','+' loses to '+','+str(Opt_Listm[0])+'\n')

							Opt_Listj[2]=-1
							break

						elif mn2>mn1+w_val:

							# if step_new==1:
							# 	st3.write(str(Opt_Listm[0])+','+' loses to '+','+str(Opt_Listj[0])+'\n')

							Opt_Listm[2]=-1

		j=0
		while j<len(GA_Info):
			GA_Infoj=GA_Info[j]
			if GA_Infoj[2]<0:
				soln_evals_tot+=GA_Infoj[1]
				soln_evals_num+=1
				# if step_new==1:
				# 	st3.write('soln_evals_tot: '+str(soln_evals_tot)+'\n')
				# 	st3.write('soln_evals_num: '+str(soln_evals_num)+'\n')
				# 	st3.write('soln_evals_avg: '+str(soln_evals_tot/soln_evals_num)+'\n')
				GA_Info.remove([GA_Infoj[0],GA_Infoj[1],GA_Infoj[2],GA_Infoj[3],GA_Infoj[4]])
				# if step_new==1:
				# 	st3.write('Removed entry '+str(j)+' from GA_PopList'+'\n')
			else:
				if step_new==1:
					st3.write('Retained sequence '+','+str(GA_Infoj)+','+' in GA_PopList'+'\n')
				j+=1

		j=0
		while j<len(Opt_List):
			Opt_Listj=Opt_List[j]
			if Opt_Listj[2]<0:
				soln_evals_tot+=Opt_Listj[1]
				soln_evals_num+=1
				# if step_new==1:
				# 	st3.write('soln_evals_tot: '+str(soln_evals_tot)+'\n')
				# 	st3.write('soln_evals_num: '+str(soln_evals_num)+'\n')
				# 	st3.write('soln_evals_avg: '+str(soln_evals_tot/soln_evals_num)+'\n')
				Opt_List.remove([Opt_Listj[0],Opt_Listj[1],Opt_Listj[2],Opt_Listj[3],Opt_Listj[4]])
				# if step_new==1:
				# 	st3.write('Removed entry '+str(j)+' from Opt_List'+'\n')
			else:
				if step_new==1:
					st3.write('Retained sequence '+','+str(Opt_Listj)+','+' in Opt_List'+'\n')
				j+=1

		if len(GA_Info)+len(Opt_List)<=Opt_Size:

			solns_left=len(GA_Info)+len(Opt_List)

			soln_evals_tot+=(solns_left*GA_LoopSize)
			soln_evals_num+=solns_left
			pruned=1

		# #Old R&S code is underneath this part

		# if step_new==1:
		# 	st3.write('GA_counter: '+','+str(GA_counter)+'\n'+'\n')
		# 	st3.write('GA_PopList:'+'\n')

		# #Find the lowest upper CI value
		# min_uppb=0
		# for j in range(len(GA_Info)):
		# 	GA_Infoj=GA_Info[j]
		# 	test_var=(GA_Infoj[4]-(GA_Infoj[1]*(GA_Infoj[2])**2))/(GA_Infoj[1]-1)
		# 	test_lowb=GA_Infoj[2]-1.96*math.sqrt(test_var/GA_Infoj[1])
		# 	test_uppb=GA_Infoj[2]+1.96*math.sqrt(test_var/GA_Infoj[1])
		# 	if step_new==1:
		# 		st3.write('Seq: '+str(GA_Infoj[0])+','+'Evals: '+str(GA_Infoj[1])+' Mn: '+str(GA_Infoj[2])+' qp[0]: '+str(GA_Infoj[3][0])+' LowB: '+str(test_lowb)+' UppB: '+str(test_uppb)+'\n')
		# 	if j==0 or test_uppb<min_uppb:
		# 		min_uppb=test_uppb

		# if step_new==1:
		# 	st3.write('\n'+'Opt_List:'+'\n')

		# for j in range(len(Opt_List)):
		# 	Opt_Listj=Opt_List[j]
		# 	test_var=(Opt_Listj[4]-(Opt_Listj[1]*(Opt_Listj[2])**2))/(Opt_Listj[1]-1)
		# 	test_lowb=Opt_Listj[2]-1.96*math.sqrt(test_var/Opt_Listj[1])
		# 	test_uppb=Opt_Listj[2]+1.96*math.sqrt(test_var/Opt_Listj[1])
		# 	if step_new==1:
		# 		st3.write('Seq: '+str(Opt_Listj[0])+','+'Evals: '+str(Opt_Listj[1])+' Mn: '+str(Opt_Listj[2])+' qp[0]: '+str(Opt_Listj[3][0])+' LowB: '+str(test_lowb)+' UppB: '+str(test_uppb)+'\n')
		# 	if test_uppb<min_uppb:
		# 		min_uppb=test_uppb

		# if step_new==1:
		# 	st3.write('\n'+' Min_uppb is '+str(min_uppb)+'\n'+'\n')

		# j=len(GA_Info)-1
		# while j>=0:
		# 	GA_Infoj=GA_Info[j]
		# 	test_var=(GA_Infoj[4]-(GA_Infoj[1]*(GA_Infoj[2])**2))/(GA_Infoj[1]-1)
		# 	test_lowb=GA_Infoj[2]-1.96*math.sqrt(test_var/GA_Infoj[1])
		# 	if test_lowb>min_uppb:
		# 		GA_Info.remove([GA_Infoj[0],GA_Infoj[1],GA_Infoj[2],GA_Infoj[3],GA_Infoj[4]])
		# 		soln_evals_tot+=GA_Infoj[1]
		# 		soln_evals_num+=1
		# 		if step_new==1:
		# 			st3.write('Removed entry '+str(j)+' from GA_PopList'+'\n')
		# 		if len(GA_Info)+len(Opt_List)==Opt_Size:
		# 			soln_evals_tot+=(Opt_Size*GA_LoopSize)
		# 			soln_evals_num+=Opt_Size
		# 			pruned=1
		# 			break
		# 	j+=-1

		# if len(GA_Info)+len(Opt_List)>Opt_Size:

		# 	j=len(Opt_List)-1
		# 	while j>=0:
		# 		Opt_Listj=Opt_List[j]
		# 		test_var=(Opt_Listj[4]-(Opt_Listj[1]*(Opt_Listj[2])**2))/(Opt_Listj[1]-1)
		# 		test_lowb=Opt_Listj[2]-1.96*math.sqrt(test_var/Opt_Listj[1])
		# 		if test_lowb>min_uppb:
		# 			Opt_List.remove([Opt_Listj[0],Opt_Listj[1],Opt_Listj[2],Opt_Listj[3],Opt_Listj[4]])
		# 			soln_evals_tot+=Opt_Listj[1]
		# 			soln_evals_num+=1
		# 			if step_new==1:
		# 				st3.write('Removed entry '+str(j)+' from Opt_List'+'\n')
		# 			if len(GA_Info)+len(Opt_List)==Opt_Size:
		# 				soln_evals_tot+=(Opt_Size*GA_LoopSize)
		# 				soln_evals_num+=Opt_Size
		# 				pruned=1
		# 				break
		# 		j+=-1

		if GA_counter>=GA_LoopSize:
			GA_CheckSize=GA_Check_Increment
			solns_left=len(GA_Info)+len(Opt_List)
			soln_evals_tot+=(solns_left*GA_LoopSize)
			soln_evals_num+=solns_left
		else:
			GA_CheckSize+=GA_Check_Increment
		if step_new==1:
			st3.write('New GA_CheckSize: '+str(GA_CheckSize)+'\n'+'\n')

	Ac_added=[]

	# if len(Opt_List)>0:
	if len(Arr_Pool)>0:

		if len(Opt_List)>0 and len(GA_Info)>0:
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

			counter=perm[1]
			qp=perm[3][0]

			#perm=Opt_List[0]
			if counter>=GA_Check_Increment or pruned==1:
				if qp>0: #0.05:
					j=0
					while perm[3][j]>0: #0.05:
						AC=perm[0][j]
						Ac_added.append(AC)
						if step_new==1:
							st3.write('Counter is '+','+str(counter)+', ss_prob is '+','+str(perm[3][j])+', Adding AC '+','+str(AC)+'\n')
						j+=1
						if j==len(perm[0]):
							break

		else:
			counter=0
			qp=0

	else:
		counter=0
		qp=0

	# end_time=time.time()
	# elap=(end_time-start_time)/conv_factor
	# #elap=0.001

	# if iter_max_d<max_d:
	# 	max_d=iter_max_d
	#print('max_d: '+str(max_d))

	#elap=fixed_elap/conv_factor

	# if len(Arr_Pool)+len(Ac_queue)>0:
	# 	print('Out of Genetic')

	return Ac_added,counter,qp,max_d,pruned,GA_CheckSize,GA_counter,soln_evals_tot,soln_evals_num