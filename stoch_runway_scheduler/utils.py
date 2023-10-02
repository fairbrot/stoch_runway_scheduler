def weather(tm,wlb,wub): #wlb is starting time for bad weather period, wub is ending time

	if wlb<wub:

		if tm<wlb:
			get_weather_state=0
		elif tm<wub:
			get_weather_state=1
		else:
			get_weather_state=2

	else:

		get_weather_state=0

	return get_weather_state

def getcost(ps_time,pool_time,trav_time,landing_time,pax_weight,thres1,thres2, lam1: float, lam2: float):

	cost=0
	#lam1=0.5 #weight for punctuality
	#lam2=0.5 #weight for queueing HMMM

	if landing_time>ps_time+thres1:
		cost+=lam1*pax_weight*(landing_time-(ps_time+thres1))**2

	if landing_time>pool_time+trav_time+thres2:
		cost+=lam2*pax_weight*(landing_time-(pool_time+trav_time+thres2))**2

	return cost



def Normal_GetServ(rel_time,trav_time,prev_class,cur_class,tm,weather_state):

	#This is for ACs that are already in the queue but not yet in service

	t1=rel_time+trav_time

	sep=Time_Sep[prev_class][cur_class]/60

	if weather_state==1:
		Mn=w_rho*sep
	else:
		Mn=sep

	SD=math.sqrt(Mn**2/k)

	z=int(random.random()*10000)

	t2=tm+normcdf[z]*SD+Mn

	if t1<t2:
		straight_into_service=0
		t_out=t2
	else:
		straight_into_service=1
		t_out=t1

	return t_out,straight_into_service


def Normal_GetServ_Future(rel_time,serv_time,trav_time,prev_class,cur_class,tm,weather_state):

	#This is for ACs that have not yet been added to the queue

	t1=rel_time+trav_time

	sep=Time_Sep[prev_class][cur_class]/60

	if weather_state==1:
		Mn=w_rho*sep
	else:
		Mn=sep

	SD=math.sqrt(Mn**2/k)

	#z=int(random.random()*10000)

	t2=tm+serv_time*SD+Mn #Look up the stored N(0,1) value, serv_time, and convert it to the correct scale

	if t1<t2:
		straight_into_service=0
		t_out=t2
	else:
		straight_into_service=1
		t_out=t1

	return t_out,straight_into_service



def Normal_Conditional_GetServ(trav_time,rel_time,sv_time,prev_class,cur_class,tm,weather_state):

	#This is for the AC currently in service

	t1=rel_time+trav_time

	#Now do the service time
	sep=Time_Sep[prev_class][cur_class]/60

	if weather_state==1:
		Mn=w_rho*sep
	else:
		Mn=sep

	SD=math.sqrt(Mn**2/k)

	t2=truncnorm(sv_time+Mn,SD,tm)

	#max_t=max(t1,t2)

	if t1<t2:
		straight_into_service=0
		t_out=t2
	else:
		straight_into_service=1
		t_out=t1

	return t_out,straight_into_service