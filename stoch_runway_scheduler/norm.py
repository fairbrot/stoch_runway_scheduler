def norm_create_cdf(k):

	t=0
	dt=0.001
	tot_x=0 #total prob so far

	p=0 #thousandth

	norm_cdf=[0]*1001
	norm_cdf[0]=0

	p=1

	while p<=999:

		#tot_x+=(t**(k-1)*math.exp(-t)/math.factorial(k-1))*dt
		tot_x+=math.exp(-(x**2/2))/(math.sqrt(2*math.pi))

		if tot_x>=p/1000:
			norm_cdf[p]=t
			p+=1
			#print(str(p))

		t+=dt

	norm_cdf[1000]=norm_cdf[999]

	return norm_cdf
