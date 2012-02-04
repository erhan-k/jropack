import pylab as py
def decodefft(finf,data, dropheights = False):
    #output: decoded data with the number of heights reduced
    #two variables are added to the finfo class:
    #deco_num_hei, deco_hrange
    #data must be arranged: 
    #    (channels,heights,times) (C-style, profs change faster)
    #fft along the entire(n=None) acquired heights(axis=1), stores in data
    num_chan = data.shape[0]
    num_ipps = data.shape[2]
    num_codes = finf.subcode.shape[0]
    num_bauds = finf.subcode.shape[1]
    NSA = finf.num_hei + num_bauds - 1
    uppower = py.ceil(py.log2(NSA))
    extra = int(2**uppower - finf.num_hei)
    NSA = int(2**uppower)
    fft_code = py.fft(finf.subcode,n = NSA,axis=1).conj()
    data = py.fft(data,n=NSA,axis=1) #n= None: no cropped data or padded zeros
    for ch in range(num_chan):
        for ipp in range(num_ipps):
            code_i = ipp % num_codes
            data[ch,:,ipp] = data[ch,:,ipp] * fft_code[code_i,:]
    data=py.ifft(data,n=NSA,axis=1) #fft along the heightsm
    if dropheights:
        return data[:,:-extra-(num_bauds-1),:]
    else:
        return data[:,:-extra,:]

def decoder_ht(data,code, iflip = 1., dropheights = False):
    """convolves each data profile with the code sequence"""
    times = data.shape[2]
    hts = data.shape[1]
    chs = data.shape[0]
    numcodes = code.shape[0]
    codelength= code.shape[1]
    code_rev=code[:,::-1] #decoding requires using the inverse of the code
    deflip = iflip
    for ch in range(chs):
        for i in range (times):
            if i % numcodes == 0 :
                deflip=deflip*iflip #call with iflip=-1 if tx has flip
            code_i = i % numcodes
            temp=py.convolve(data[ch,:,i],code_rev[code_i,:])
            data[ch,:,i]=deflip*temp[codelength-1:codelength+hts]
    return data

def decoder(data,code,iflip):
    """convolves each data profile with the code sequence"""
    times=py.shape(data)[0]
    hts=py.shape(data)[1]
    codelength=py.shape(code)[0]
    code_rev=code[::-1] #decoding requires using the inverse of the code
    deflip=1
    #pdb.set_trace()
    for i in range (times):
        temp=py.convolve(data[i,:],code_rev)
        data[i,:]=deflip*temp[codelength-1:codelength+hts]
        deflip=deflip*iflip #call with iflip=-1 if tx has flip
        #pdb.set_trace()
    return data

def test_decodefft():
    from time import time
    class finf:
        pass
    finf.num_chan = 1
    finf.num_hei = 300
    numbauds = 13
    finf.subcode = py.array([[1.,1,1,1,1,-1,-1,1,1,-1,1,-1,1],
            [-1.,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1]])
    data1 = py.random([6,finf.num_hei,128])
    factor = 5.
    data1[0,:numbauds,0] += factor *finf.subcode[0]
    data1[0,50:50+numbauds,0] += factor * finf.subcode[0]
    
    data2 = py.random([6,finf.num_hei,128])
    data2[0,50:50+numbauds,0] += factor * finf.subcode[0]
    datas = (data1,data2)
#    prevfs = py.rcParams['font.size']
    py.rcParams['font.size'] = 8 #this is just to make label lines close
    for i in range(2):
        data = datas[i].copy()
        fig = py.figure()
        nrows = 5
        
        ax = fig.add_subplot(nrows,1,1)
        ax.plot(data[0,:,0],'.-');ax.legend(['data'],loc=2)
        ax.set_xlim([0,finf.num_hei])
        
        ax = fig.add_subplot(nrows,1,2)
        decodedata3 = py.empty_like(data)
        stime = time()
#        for ch in range(data.shape[0]):
#            decodedata3[ch,:,:] = decoder_ht(data[ch,:,:],finf.subcode[0],-1)
        decodedata3 = decoder_ht(data,finf.subcode,1)
        dtime = time() - stime
        ax.plot(abs(decodedata3[0,:,0]),'.-');
        ax.legend(['decoder %f ms'%(dtime*1000)],loc=2)
        ax.set_xlim([0,finf.num_hei])
        #data with code in the middle
        
        data = datas[i].copy()
        stime = time()
        decodedata1 = decodefft(finf,data)
        dtime = time() - stime
        ax = fig.add_subplot(nrows,1,3)
        ax.plot(abs(decodedata1[0,:,0]),'.-');
        ax.legend(['decodefft %f ms'%(dtime*1000)],loc=2)
        ax.set_xlim([0,finf.num_hei])
        data = datas[i].copy()
        decodedata2 = decodefft(finf,data,dropheights = True)

        data = datas[i].copy()
        ax = fig.add_subplot(nrows,1,4)
        ax.plot(abs(decodedata2[0,:,0]),'.-');
        ax.legend(['decodefft droping'],loc=2)
        ax.set_xlim([0,finf.num_hei])
        
        ax = fig.add_subplot(nrows,1,5)
        ax.plot(abs(decodedata1[0,:,0])-abs(decodedata3[0,:,0]),'.-');
        ax.legend(['diff. decoder & decodefft'],loc=2)
        ax.set_xlim([0,finf.num_hei])
    py.draw()
    fig.show()
#    py.rcParams['font.size'] = prevfs

    
