from scipy import constants
import pylab as py
import numpy as np
import time as tmm

def decode(data):
	#data must be arranged: (channels,heights,profiles) (C-style, profs change faster)
	data=np.fft.fft(data,None,1) #fft along the heights
	NSA=self.Hinfo.num_hei
	num_codes=self.Hinfo.subcode.shape[0]
	num_bauds=self.Hinfo.subcode.shape[1]
	for code_ind in range(num_codes):
		for hei_ind in range(NSA):
			data[:,hei_ind,code_ind::num_codes]=\
			data[:,hei_ind,code_ind::num_codes]*self.fft_code[code_ind,hei_ind]
	data=np.fft.ifft(data,None,1) #fft along the heightsm
	return data[:,:-num_bauds+1,:]

def enoise(data,navg):
	"""[noise,stdv,i_max,index]=enoise(data,navg)
	Routine to estimate noise level based on Hildebrand and Sekhon
	Modified by Pablo Reyes 15apr2009
	"""
	index=data.argsort()
	data=data[index]
	npts=data.size
	
#   count=max(min(10.,npts),int(npts/50.))
	count=min(10.,npts/2)
	s=data[0:count].sum()
	s2=(data[0:count]**2).sum()
	p=s/(count+1) #added for the case of count==npts
	dp2=(s2-(count+1)*p**2)/count #added for the case of count==npts
	for id in np.arange(count,npts): #last id=npts-1
		s=s+data[id]
		s2=s2+data[id]**2
		p=s/(id+1)
		dp2=(s2-(id+1)*p**2)/id
		if p**2/navg < dp2:
			noise=data[0:id].sum()/id
			stdv=dp2**0.5
			i_max=id
			index=index[0:i_max]
			return noise,stdv,i_max,index
	noise=data.sum()/npts
	stdv=dp2**0.5
	i_max=npts
	return noise,stdv,i_max,index
	

