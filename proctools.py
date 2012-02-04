def decodefft(finfo,data, chans=None):
    if chans==None: chans = finfo.num_chan
    from numpy.fft import fft,ifft
    #output: decoded data with the number of heights reduced
    #two variables are added to the finfo class:
    #deco_num_hei, deco_hrange
    #data must be arranged either: 
    #    (channels,heights,profiles) (C-style, profs change faster)
    #or: (profiles,heights,channels) (C-style, profs change faster)
    #fft along the entire(n=None) acquired heights(axis=1), stores in data
    data=fft(data,n=None,axis=1) 
    NSA=finfo.num_hei
    num_codes=finfo.subcode.shape[0]
    num_bauds=finfo.subcode.shape[1]
    fft_code=fft(finfo.subcode,finfo.num_hei,1).conj()
    #if the order is (channels,heights,profiles)
    if data.shape[0]==chans:  
        for code_ind in range(num_codes):
            for hei_ind in range(NSA):
                data[:,hei_ind,code_ind::num_codes]=\
                data[:,hei_ind,code_ind::num_codes]*fft_code[code_ind,hei_ind]
    else:  #if the order is (profiles,heights,channels)
        for code_ind in range(num_codes):
            for hei_ind in range(NSA):
                data[code_ind::num_codes,hei_ind,:]=\
                data[code_ind::num_codes,hei_ind,:]*fft_code[code_ind,hei_ind]
    data=ifft(data,None,1) #fft along the heightsm
    return data[:,:-num_bauds+1,:]

def MSTISR_corruptdata(year,doy,intg_time):
    from numpy import loadtxt
    from os import path
    filetoread=path.expanduser(
        '~/Research/pythonWC/data/corruptdata/%.4d%.3d_%dsec.txt'%(year, 
            doy, intg_time))
    if path.exists(filetoread):
        return loadtxt(filetoread)
    else:
        print 'File not found ... %s'%filetoread
        return None
        
def enoise3(data,navg):
    """    [noise,stdv,i_max,index]=enoise3(data,navg)
    Routine to estimate noise level based on Hildebrand and Sekhon
    Modified by Pablo Reyes 15apr2009
    """
    import numpy as np
    # pdb.set_trace()
    index = data.argsort()
    data = data[index]
    npts = data.size    
    count=min(npts,max(int(npts/50.),10))
    s = data[0:count].sum()
    s2 = (data[0:count] ** 2.).sum()
    for id in np.arange(count, npts): #last id=npts-1
        s = s + data[id]
        s2 = s2 + data[id] ** 2.
        p = s / (id + 1)
        dp2 = (s2 - (id + 1) * p ** 2.) / id
        if (p) ** 2./navg < dp2:
            noise = data[0:id].sum() / id
            stdv = dp2 ** 0.5
            i_max = id
            index = index[0:i_max]
            return noise, stdv, i_max, index
        if False: # Put True if you want to see the progression
            import pylab as py
            pwr=data
            noise=data[0:id].sum() / id
            stdn = dp2 ** 0.5
            py.figure(200);py.clf()
            py.plot(pwr);
            py.plot([0,200],[noise,noise]);
            py.plot([0,200],py.array([noise,noise])+stdn);
            py.plot([0,200],py.array([noise,noise])-stdn);
            py.ylim(noise+4*stdn*py.array([-1,1]))
            py.title('noise:%f, stdn:%f'%(noise,stdn))
            py.draw()
            py.show()
    noise = data.sum() / npts
    stdv = dp2 ** 0.5
    stdv = dp2
    i_max = npts
    return noise, stdv, i_max, index
