from typing import List
import random
import math

def sample_cond_gamma(t,gamma_cdf):

	#print('sample_cond_gamma')

	if t==0:
		mm=int(1000*random.random())

	else:
		#We want to sample a gamma RV which exceeds t

		lb=0
		ub=1000

		chk=0
		while chk==0:

			v=lb+int((ub-lb)/2) #bisection search

			#print('t: '+str(t)+' v: '+str(v)+' lb: '+str(lb)+' ub: '+str(ub))

			if gamma_cdf[v]<t and gamma_cdf[v+1]>=t:
				chk=1
			elif v==999:
				chk=1
			elif gamma_cdf[v]<t:
				lb=v 
			else:
				ub=v

		v+=1

		nn=(999-v)+1 #no. of random possibilities

		z=int(nn*random.random())
		mm=v+z #value to look up

	y=gamma_cdf[mm]

	#print('Out of sample_cond_gamma')

	return y


def Gamma_GetServ(k: int, Time_Sep: List[List[int]], rel_time,trav_time,prev_class,cur_class,tm,weather_state,gamma_cdf, w_rho: float):

	#This is for ACs that are already in the queue but not yet in service

	t1=rel_time+trav_time

	rate=k/(Time_Sep[prev_class][cur_class]/60)
	if weather_state==1:
		rate*=1/w_rho #=0.5

	z=int(random.random()*1000)
	getserv=gamma_cdf[z]
	getserv*=1/rate #convert to correct scale

	t2=tm+getserv

	if t1<t2:
		straight_into_service=0
		t_out=t2
	else:
		straight_into_service=1
		t_out=t1

	return t_out,straight_into_service

def Gamma_GetServ_Future(k: int, Time_Sep: List[List[int]], rel_time,serv_time,trav_time,prev_class,cur_class,tm,weather_state, w_rho: float):

	#This is for ACs that have not yet been added to the queue

	t1=rel_time+trav_time

	rate=k/(Time_Sep[prev_class][cur_class]/60)
	if weather_state==1:
		rate*=1/w_rho #=0.5

	getserv=serv_time #Look up the stored Gamma(k,1) value, serv_time
	getserv*=1/rate #Convert it to the correct scale

	t2=tm+getserv

	if t1<t2:
		straight_into_service=0
		t_out=t2
	else:
		straight_into_service=1
		t_out=t1

	return t_out,straight_into_service

def Gamma_Conditional_GetServ(k: int, Time_Sep: List[List[int]], trav_time,rel_time,sv_time,prev_class,cur_class,tm,weather_state,gamma_cdf, w_rho: float):

	#This is for the AC currently in service

	t1=rel_time+trav_time

	rate=k/(Time_Sep[prev_class][cur_class]/60)
	if weather_state==1:
		rate*=1/w_rho #=0.5

	cond_serv=sample_cond_gamma(rate*(tm-sv_time),gamma_cdf)
	cond_serv*=1/rate #Convert back to the correct scale

	t2=sv_time+cond_serv

	#max_t=max(t1,t2)

	if t1<t2:
		straight_into_service=0
		t_out=t2
	else:
		straight_into_service=1
		t_out=t1

	return t_out,straight_into_service

def gamma_cond_exp(t,alpha,beta):

	#Site: https://stats.stackexchange.com/questions/338378/closed-form-conditional-expectation-of-gamma-distributed-variable

	#This is assuming mean is alpha/beta, variance is alpha/(beta^2)
	#t is the amount of time that the service has already been in progress

	denom=0
	for j in range(alpha):
		denom+=((beta*t)**j)/math.factorial(j)
	num=((beta*t)**alpha)/math.factorial(alpha)
	x=(alpha/beta)*(1+num/denom)

	#print('t: '+str(t)+' alpha: '+str(alpha)+' beta: '+str(beta)+' x: '+str(x))

	return x