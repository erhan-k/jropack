# jroplot.py
#
#
# Created by Pablo M. Reyes Firpo on 9/12/11.
# Copyright 2011. All rights reserved.
#
# This file contains routines to plot spc and cspc arrays with information
# from the finfo class of the jroread.py routines
#
# History:
#   09/12/11:   by P. Reyes
#       - plotspc, and plotcspc routines
#   09/18/11:   by P. Reyes
#       - Adding the plotrti routine aith dBrange 
#       - Adding dBrange to the other plotting routines
#   10/01/11:   by P. Reyes
#       - Modifying plotpsc2 to plot the noise level in the spc graphs
#       - Adding plotmap for RTI, Coherence, and SNR maps
#   10/08/11:   by P. Reyes
#       - plotmap plots SNR together with vectormap (U and V velocities)
#   10/11/11:   by P. Reyes
#       - plotmap: Added some help documentation
#       - plotmap: The upper noise plot y limits is now calculated from the
#           median of the noise points. Before it was from the minimum, and 
#           the data errors where affecting the limits in some occations.
#   10/24/11:   by P. Reyes
#       - plotcoherence routine

from pylab import *
from time import strftime, gmtime
#from pdb import set_trace
#from IPython.Debugger import Tracer; debug_here = Tracer()
from calendar import timegm
import matplotlib.gridspec as gridspec # useful for axes with diff ratios

#From:
#http://matplotlib.sourceforge.net/examples/pylab_examples/custom_cmap.html
"""
Example: suppose you want red to increase from 0 to 1 over the bottom
half, green to do the same over the middle half, and blue over the top
half.  Then you would use:

cdict = {'red':   ((0.0,  0.0, 0.0),
                   (0.5,  1.0, 1.0),
                   (1.0,  1.0, 1.0)),

         'green': ((0.0,  0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.75, 1.0, 1.0),
                   (1.0,  1.0, 1.0)),

         'blue':  ((0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (1.0,  1.0, 1.0))}

If, as in this example, there are no discontinuities in the r, g, and b
components, then it is quite simple: the second and third element of
each tuple, above, is the same--call it "y".  The first element ("x")
defines interpolation intervals over the full range of 0 to 1, and it
must span that whole range.  In other words, the values of x divide the
0-to-1 range into a set of segments, and y gives the end-point color
values for each segment.

Now consider the green. cdict['green'] is saying that for
0 <= x <= 0.25, y is zero; no green.
0.25 < x <= 0.75, y varies linearly from 0 to 1.
x > 0.75, y remains at 1, full green.

If there are discontinuities, then it is a little more complicated.
Label the 3 elements in each row in the cdict entry for a given color as
(x, y0, y1).  Then for values of x between x[i] and x[i+1] the color
value is interpolated between y1[i] and y0[i+1].

Going back to the cookbook example, look at cdict['red']; because y0 !=
y1, it is saying that for x from 0 to 0.5, red increases from 0 to 1,
but then it jumps down, so that for x from 0.5 to 1, red increases from
0.7 to 1.  Green ramps from 0 to 1 as x goes from 0 to 0.5, then jumps
back to 0, and ramps back to 1 as x goes from 0.5 to 1.

row i:   x  y0  y1
                /
               /
row i+1: x  y0  y1

Above is an attempt to show that for x in the range x[i] to x[i+1], the
interpolation is between y1[i] and y0[i+1].  So, y0[0] and y1[-1] are
never used.
"""

cdict3 = {'red':  ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.8, 1.0),
                   (0.75,1.0, 1.0),
                   (1.0, 0.4, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.9, 0.9),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.4),
                   (0.25,1.0, 1.0),
                   (0.5, 1.0, 0.8),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }
cdict3_reversed = {
          'red':  ((0.0 , 0.5, 0.5),
                   (0.25, 1.0, 1.0),
                   (0.5 , 0.9, 0.9),
                   (0.55, 0.0, 0.0),
                   (1.0 , 0.0, 0.0)),

         'green': ((0.0 , 0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.5 , 0.8, 0.8),
                   (0.75, 0.0, 0.0),
                   (1.0 , 0.0, 0.0)),

         'blue':  ((0.0 , 0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.5 , 0.9, 0.9),
                   (0.75, 1.0, 1.0),
                   (1.0 , 0.5, 0.5))
        }
register_cmap(name='BlueRed3', data=cdict3)
register_cmap(name='BlueRed3_reversed', data=cdict3_reversed)

def plotspec(finf, data, fid=None, chorder=None, chname=None,  dBrange=None,
        noise=None, dshift=None, dshift_err=None, dpi = 80., rowHfactor = 1.0,
        chs2show = None):
    #fid : is the figure number, if None then a new figure number is assign
    if chorder==None: chorder = range(int(finf.num_chan))
    if chname==None:
        chname = []
        for q in range(int(finf.num_chan)):
            chname.append(str(q))
    if data.shape[1] == finf.hrange.size:
        hrange = finf.hrange
    else:
        hrange = arange(data.shape[1])
        if hasattr(finf, 'deco_hrange'):
            if data.shape[1] == finf.deco_hrange.size:
                hrange = finf.deco_hrange
    finf.LTarr = gmtime(finf.shortH['ltime']-finf.timezone_sec)
    fontsize = 8. * 80. / dpi
    if chs2show == None:
        nchan=finf.num_chan
        chs2show = range(nchan)
    else:
        nchan = len(chs2show)
    ncols = 1 + (nchan >= 2) #1 or 2. "1" only if nchan = 1
    nrows = ceil(nchan/2.)
    specH =  rowHfactor * 2.25 * 80. / dpi #height of a spec graph
    TB = (1 + rowHfactor * 0.18) * 0.65  * 80. / dpi #Top and Bottom space
    DH = 0.75 * 80. / dpi #space between rows of spec graphs
    WD = ncols * 5. * 80. / dpi #Total width of two channels
    Hn = nrows * specH + 2 * TB + (nrows - 1) * DH # total height of figure
    TB_r = TB / Hn #Top and Bottom ratio with respect to the total height
    DH_r = DH / specH
    fig = figure(fid, figsize=(WD, Hn), dpi=dpi)

    pwrprofile = data.mean(2)
    if dBrange == None:
        if noise == None:
            # smallsig_ind = find((finf.hrange>60) & (finf.hrange<800))
            # vmin = 10*log10(nanmin(pwrprofile[chs2show,smallsig_ind]))
            # vmax = 10*log10(nanmax(pwrprofile[chs2show,smallsig_ind]))
            vminmax = vmin = 10*log10(nanmin(pwrprofile[chs2show,:]))
            for i in range(nchan):
                if 10*log10(nanmin(pwrprofile[chs2show[i],:])) > vminmax:
                    vminmax = 10*log10(nanmin(pwrprofile[chs2show[i],:]))
            vmin = vmin - 0.5
            vmax = vminmax + 7.5
        else:
            vmin = 10*log10(nanmin(noise[chs2show])) - 0.5
            vmax = 10*log10(nanmax(noise[chs2show])) + 7.5
    else:
        vmin = dBrange[0]
        vmax = dBrange[1]
    fig.clf()
    specitems = 4
    gs = gridspec.GridSpec(int(nrows) , int(ncols * specitems),
            width_ratios=[.125,.37,.025,.15] + (ncols-1) * [.125,.37,.025,0])
    gs.update(left=0.1,right=0.95,hspace=DH_r,wspace=0,
            bottom = TB_r, top = (1.-TB_r)  )
    pwrcolors=['b-','r-','g-','y-','c-','m-']
    data[data<=0]=nan #to prevent zeors or negative numbers to be taken log10
    for ch in range(nchan):
        ax = fig.add_subplot(gs[specitems*ch+1])
        pcm = ax.pcolormesh(finf.vel_arr, hrange,
                            10*log10(data[chorder[chs2show[ch]],:,:]),
                vmin=vmin,vmax=vmax)
        if dshift != None and dshift_err != None:
            vm = ax.errorbar(dshift[chs2show[ch],:], hrange,
                             xerr = dshift_err[chs2show[ch],:],
                             fmt = 'k.-', ms = 0.5)
        elif dshift != None and dshift_err == None:
            vm = ax.plot(dshift[chs2show[ch],:],hrange,'k.-',ms = 0.5)
        ax.set_xlim(finf.vel_arr[[0,-1]])
        ax.set_ylim(hrange[[0,-1]])

        cbar = fig.colorbar(pcm,
                cax=axes(gs[specitems*ch+2].get_position(fig)))
#        cbar = fig.colorbar(pcm)
        cbar.ax.tick_params(labelsize = fontsize)
        cbar.ax.set_title('dB',fontsize=fontsize)
        bname = 'Bm '+chname[chorder[chs2show[ch]]]
        ax.set_title(bname+strftime(' %y-%b-%d %H:%M:%S L.T.',finf.LTarr),
                fontsize = fontsize)
        ax.set_xlabel('Doppler Velocity [m/s]',fontsize=fontsize)
        #ax.set_ylabel('Height [km]',fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        ax.set_yticklabels([])
        #ax2 will be the axis of power profile
        ax2 = fig.add_subplot(gs[specitems*ch])
        if noise!=None:
            nleveldB = 10*log10(noise[chorder[chs2show[ch]]])
            ax2.plot(nleveldB*array([1,1.]), [hrange[0],hrange[-1]],'k',
                    dashes=[5,2],#5 points on, 2 points off
                    linewidth=0.8)
            ax2.set_title('Noise %.2f dB'%nleveldB, fontsize=fontsize)
        else:
            ax2.set_title('Power Profile %s'%chname[chorder[chs2show[ch]]],
                    fontsize=fontsize)
        if ch == 0:
            pwrprofile[pwrprofile<=0]=nan
            for chn in range(1,nchan):
                ax2.plot(10*log10(pwrprofile[chorder[chs2show[chn]],:]),hrange,
                        pwrcolors[chn % size(pwrcolors)])
        ax2.plot(10*log10(pwrprofile[chorder[chs2show[ch]],:]),hrange,
                pwrcolors[chs2show[ch] % size(pwrcolors)])
        ax2.set_xlim([vmin,vmax])
        ax2.set_ylim(hrange[[0,-1]])
        ax2.grid(True)
        ax2.set_xlabel('Power [dB]',fontsize=fontsize)
        #ax2.set_yticklabels([])
        ax2.set_ylabel('Height [km]',fontsize=fontsize)
        ax2.tick_params(labelsize=fontsize)
    fig.suptitle('Spectra Plots, %s'%strftime(
            ' %Y-%b-%d %H:%M:%S LT',finf.LTarr),fontsize=fontsize*1.2)
    draw()
    fig.show()
    return fig

def plotcoherence(finf, cspc, spc, fid=None, chname=None, facecolor=None,
        drawinside=True, dpi = 80., pairs2show = None, crange = None,
        rowHfactor = 1.0):
#    coherence = zeros_like(cspc)
    if pairs2show == None:
        pairs = finf.pairs
    else:
        pairs = array(pairs2show)
    if pairs.shape == (2,) : pairs.shape = (1,2)
    coherence = nan * empty((pairs.shape[0],)+cspc.shape[1:],dtype=cspc.dtype)
    for pp in range(pairs.shape[0]):
        denominator = sqrt(prod(spc[pairs[pp],:,:],axis=0))
        for i in range(finf.pairs.shape[0]):
            if (pairs[pp] == finf.pairs[i]).prod() == 1:
                coherence[pp,:,:] = cspc[i,:,:] / denominator
    if chname==None:
        chname = []
        for q in range(finf.num_chan):
            chname.append(str(q))
    finf.LTarr = gmtime(finf.shortH['ltime']-finf.timezone_sec)

    if spc.shape[1] == finf.hrange.size:
        hrange = finf.hrange
    else:
        hrange = arange(spc.shape[1])
        if hasattr(finf, 'deco_hrange'):
            if spc.shape[1] == finf.deco_hrange.size:
                hrange = finf.deco_hrange
    fontsize = 7.5 * 80. / dpi
    ncols = 2 * (1 + (pairs.shape[0] >= 2)) #1 or 2. "1" only if num pairs = 1
    nrows = ceil(pairs.shape[0]/2.)
    spcH = rowHfactor * 1.3 * 80. / dpi #height of a spc graph
    TB = (1 + rowHfactor * 0.18) * 0.65 * 80. / dpi #Top and Bottom space
    DH = 0. * 80. / dpi #space between rows of spc graphs
    WD = ncols * 3.5 * 80. / dpi #Total width of two channels
    Hn = nrows * spcH + 2 * TB + (nrows - 1) * DH # total height of figure
    TB_r = TB / Hn #Top and Bottom ratio with respect to the total height
    DH_r = DH / spcH
    fig = figure(fid, figsize=(WD, Hn), dpi=dpi, facecolor=facecolor)
    crosspwr = absolute(coherence)
    crossang = angle(coherence)
    pwrprofile = crosspwr.mean(2)
    angprofile = crossang.mean(2)
    if crange == None:
        vmin,vmax = (0., 1.)
    else:
        vmin,vmax = crange
    fig.clf()
    cspecitems = 4
    gs = gridspec.GridSpec(int(nrows) , int(ncols * cspecitems),
            width_ratios=[0.125,0.25, 0.025, 0.05, 0.125,0.25, 0.025,0.15]
            + (ncols/2 - 1) * [0.125,0.25, 0.025, 0.05, 0.125,0.25, 0.025,0.0])
    gs.update(left=0.1,right=0.95,hspace=DH_r,wspace=0,
            bottom = TB_r, top = (1.-TB_r)  )
    pwrcolors=['b-','r-','g-','y-','c-','m-']
#    crosspwr[crosspwr<=0]=nan #prevent zeros or neg. numbers to be taken log10
    for pp in range(pairs.shape[0]):
        ax = fig.add_subplot(gs[2 * cspecitems*pp+1])
        pcm = ax.pcolormesh(finf.vel_arr,hrange,
                crosspwr[pp,:,:],vmin=vmin,vmax=vmax)
        ax.set_xlim(finf.vel_arr[[0,-1]])
        ax.set_ylim(hrange[[0,-1]])
        cbar = fig.colorbar(pcm,
                cax=axes(gs[2*cspecitems*pp+2].get_position(fig)))
        cbar.ax.tick_params(labelsize = fontsize)
        bname = 'pair (%s,%s)'%tuple(array(chname)[pairs[pp,:]])
        if pp in [0,1]:
#            cbar.ax.set_title('coh.',fontsize=fontsize)
            ax.set_title('Coherence',fontsize = fontsize)
        if pp in [2*((pairs.shape[0]-1)/2), 2*((pairs.shape[0]-1)/2)+1]:
            ax.set_xlabel('Doppler Velocity [m/s]',fontsize=fontsize)
        else:
            ax.set_xticklabels([])
        #ax.set_ylabel('Height [km]',fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        ax.set_yticklabels([])
        ax2 = fig.add_subplot(gs[2 * cspecitems*pp])
        if pp == 0:
            for ppn in range(1,pairs.shape[0]):
                ax2.plot(pwrprofile[ppn,:],hrange,
                        pwrcolors[ppn % size(pwrcolors)])
        if pp in [0,1]:
            ax2.set_title('Coherence',fontsize=fontsize)
        ax2.plot(pwrprofile[pp,:],hrange,
                pwrcolors[pp % size(pwrcolors)])
        ax2.set_xlim([vmin,vmax])
        ax2.set_ylim(hrange[[0,-1]])
        ax2.grid(True)
        if pp in [2*((pairs.shape[0]-1)/2), 2*((pairs.shape[0]-1)/2)+1]:
            ax2.set_xlabel('Coherence ',fontsize=fontsize)
        else:
            ax2.set_xticklabels([])
        #ax2.set_yticklabels([])
        ax2.set_ylabel('Height [km]\n%s'%bname,fontsize=fontsize)
        ax2.tick_params(labelsize=fontsize)
        ax = fig.add_subplot(gs[2 * cspecitems*pp+5])
        pcm = ax.pcolormesh(finf.vel_arr,hrange,
                crossang[pp,:,:] * 180 / pi,vmin=-180,vmax=180)
        ax.set_xlim(finf.vel_arr[[0,-1]])
        ax.set_ylim(hrange[[0,-1]])
        cbar = fig.colorbar(pcm,
                cax=axes(gs[2*cspecitems*pp+6].get_position(fig)))
        cbar.ax.tick_params(labelsize = fontsize)
        bname = 'pair (%s,%s)'%tuple(array(chname)[pairs[pp,:]])
        if pp in [0,1]:
            cbar.ax.set_title('Phase',fontsize=fontsize)
            ax.set_title('NCS Phase',fontsize = fontsize)
        if pp in [2*((pairs.shape[0]-1)/2), 2*((pairs.shape[0]-1)/2)+1]:
            ax.set_xlabel('Doppler Velocity [m/s]',fontsize=fontsize)
        else:
            ax.set_xticklabels([])
        #ax.set_ylabel('Height [km]',fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        ax.set_yticklabels([])
        ax2 = fig.add_subplot(gs[2 * cspecitems*pp+4])
        if pp == 0:
            for ppn in range(1, pairs.shape[0]):
                ax2.plot(angprofile[ppn,:] * 180. / pi,hrange,
                        pwrcolors[ppn % size(pwrcolors)])
        if pp in [0,1]:
            ax2.set_title('Phase',fontsize=fontsize)
        ax2.plot(angprofile[pp,:] * 180. / pi,hrange,
                pwrcolors[pp % size(pwrcolors)])
        ax2.set_xlim([-180,180])
        ax2.set_ylim(hrange[[0,-1]])
        ax2.grid(True)
        if pp in [2*((pairs.shape[0]-1)/2), 2*((pairs.shape[0]-1)/2)+1]:
            ax2.set_xlabel('Phase [deg]',fontsize=fontsize)
        else:
            ax2.set_xticklabels([])
        #ax2.set_yticklabels([])
#        ax2.set_ylabel('Height [km]\n%s'%bname,fontsize=fontsize)
        ax2.tick_params(labelsize=fontsize)
        ax2.set_yticklabels([])

    fig.suptitle('Normalized Cross-Spectrum (NCS), %s'%strftime(
            ' %Y-%b-%d %H:%M:%S LT',finf.LTarr),fontsize=fontsize*1.2)
    draw()
    fig.show()
    return fig


def mark_large_gaps(rti,acqtime_arr,noise=None):
    timediff = acqtime_arr[1:] - acqtime_arr[:-1]
    meantdiff = timediff[~isnan(timediff)].mean()
    markindex = find(timediff > 2 * meantdiff)
#    markindex = find(timediff > 2 * nanmin(timediff))
    interp_time = acqtime_arr[markindex] + meantdiff
    acqtime_arr = insert(acqtime_arr,markindex+1,interp_time)
    rti = insert(rti,markindex+1,nan,axis = 2)
    if noise != None:
        noise = insert(noise,markindex+1,nan,axis = 1)
        return rti, acqtime_arr, noise
    else:
        return rti, acqtime_arr

def plotmap(finf=None, data=None, fid=None, acqtime_arr=None, chname=None,
        chorder=None, crange=None,facecolor=None,drawinside=True,maptype='RTI',
        mapunits='dB', tlim=None, hlim=None, noise=None, rowHfactor=1.0,
        cmap = None, cmap_reversed = False, vectormap=None, dpi = 80,
        order = 'ch_h_t', colWfactor=1.0):
    """
    Inputs
     finf:
        finfo class with the parameters of the experiment.
        If finfo=None, channels will be obtained from the data shape, hrange
        will be consecutive integers and timezone will be -5 hours(Jicamarca)
     data:
        2D or 3D array containing RTI map.
            2D: heights x times
            3D: channels x heights x times
    Example1 (Simple use):
     jroplot.plotmap(None, data)

    Example2 :
     jroplot.plotmap(finf, data)

    Comments:
        P. Reyes (11oct2011): The order of the input parameters should be:
            data, finf=None, fid=None, ......
    """
    #fid : is the figure number, if None then a new figure number is assign
    # data expected to be: (ch, heights, time)
    # if acqtime_arr is None then a simple array of numbers is created
    if data.ndim == 2: # if dimensions = 2 :only 1 channel data
        data.shape = (1,) + data.shape
    if finf==None: #If not finf given, then generate a fake finf
        from time import timezone
        class finf:
            pass
        finf.num_chan = data.shape[0]
        finf.hrange = arange(data.shape[1])
        finf.timezone_sec = 5 * 60 *60. #At Jicamarca
    if data.shape[0] != finf.num_chan:
        order = 't_ch_h'
    if order == 't_ch_h':
        data = data.transpose([1,2,0])
        if noise != None:
            noise = noise.transpose([1,0])
        if vectormap != None:
            vectormap = vectormap.transpose([1,2,0])
    if data.shape[1] == finf.hrange.size:
        hrange = finf.hrange
    else:
        hrange = arange(data.shape[1])
        if hasattr(finf, 'deco_hrange'):
            if data.shape[1] == finf.deco_hrange.size:
                hrange = finf.deco_hrange

    nchan = int(finf.num_chan)
    if cmap == None:
        if cmap_reversed:
            cmap = cm.jet_r
        else:
            cmap = cm.jet
    if chorder==None: chorder = range(nchan)
    if chname == None:
        chname = []
        for q in range(nchan):
            chname.append('Ch:%s'%str(q))
    if hlim==None:
        hlim = hrange[[0,-1]]

    if acqtime_arr == None:
        if order == 't_ch_h':
            datearr = arange(data.shape[0])
        elif order == 'ch_h_t':
            datearr = arange(data.shape[2])
        sdata = ''
        txrange = array(datearr[[0,-1]])
    else:
        #modify data to insert no numbers (nans) in between large gaps of data.
        if vectormap != None:
            vectormap, junk = mark_large_gaps(vectormap,acqtime_arr)
        if noise != None:
            data, acqtime_arr,noise = mark_large_gaps(data,acqtime_arr,noise)
        else:
            data, acqtime_arr = mark_large_gaps(data,acqtime_arr)
        # time formatter
        #######################################################################
        StarT_t = gmtime(acqtime_arr[0] - finf.timezone_sec)
        LTstartday = timegm([StarT_t.tm_year, 1,StarT_t.tm_yday ,0 ,0 ,0])
        LT1 = gmtime(LTstartday)  #getting the start time array
        dt  = date2num(datetime.datetime(LT1.tm_year, LT1.tm_mon, LT1.tm_mday))
        timevec_hr = (acqtime_arr - LTstartday - finf.timezone_sec) / 3600.
        datearr = dt + timevec_hr / 24. # to display time in the x axis
        dateFmt = DateFormatter('%H:%M')
        #######################################################################
        if tlim == None:
            StopT_t = gmtime(acqtime_arr[-1] - finf.timezone_sec)
            LTstopday = timegm([StopT_t.tm_year, 1,StopT_t.tm_yday ,0 ,0 ,0])
            days2display = 1 + ceil((LTstopday - LTstartday)/(24.*3600.))
            tlim = [0., 24.* days2display]
        if tlim[1]>24.:
            days2display = ceil(tlim[1] / 24.)
            LT2 = gmtime(LTstartday + 3600.*24.*(days2display-1))
            sdate = strftime('%Y-%b-%d to ',LT1) + strftime('%Y-%b-%d',LT2)
        else:
            sdate = strftime('%Y-%b-%d', LT1)
        txrange = dt + array(tlim)/24.
    if crange != None:
        if size(crange) == 2:
            vmin = [crange[0]]*nchan
            vmax = [crange[1]]*nchan
        elif size(crange)/2 == nchan and crange.__len__() == nchan:
            vmin = []
            vmax = []
            for q in range(nchan):
                vmin.append(crange[q][0])
                vmax.append(crange[q][1])
        else:
            crange = None #So that it gets calculated :
    if crange == None: #This redundancy of question is in case crange is bad
        smallsig_ind = find((hrange>60) & (hrange<800))
        if mapunits=='dB':
            valmin = 10*log10(nanmin(data[:,smallsig_ind,:]))
            valmax = 10*log10(nanmax(data[:,smallsig_ind,:]))
        else:
            valmin = (nanmin(data[:,smallsig_ind,:]))
            valmax = (nanmax(data[:,smallsig_ind,:]))
        vmin = [valmin]*int(finf.num_chan)
        vmax = [valmax]*int(finf.num_chan)
    fontsize = 8 * 80. / dpi
    ncols = 1+1 #extra one for the colorbar
    noiseindex = 0
    nrows = finf.num_chan
    chrowindexes = range(int(finf.num_chan))
    if maptype[:3] in ['RTI','SNR']:
        if noise != None:
            nrows = finf.num_chan + 1
            noiseindex = 1
    elif maptype[:3]=='Coh':
        nrows = finf.num_pairs
        if finf.pairs.shape == (2,) : finf.pairs.shape = (1,2)
        chrowindexes = range(int(finf.num_pairs))

    rowH = rowHfactor * 1.4 * 80. / dpi #height of a row
    TB = (1 + rowHfactor * 0.18) * 0.65 * 80. / dpi #Top and Bottom space
    DH = 0. * 80. / dpi #space between rows of spc graphs
    Hn = nrows * rowH + 2 * TB + (nrows - 1) * DH # total height of figure
    #Top and Bottom ratio with respect to the total height
    TB_r = TB / Hn
    DH_r = DH / rowH
    # colorbar
    CW = 0.2 * 80 / dpi #colorbar width
    LW = 0.9 * 80 / dpi #Left space
    RW = 0.6 * 80 / dpi #Right space
    DW = 0.5 * 80. / dpi #space between cols graph and colorbar
    LegW = 0.3 * 80 / dpi #leggend spacing
    PW = colWfactor*(ncols-1)*7 * 80. / dpi #Total plot width of the plot
    GW = PW + DW + CW # Plot width + Colorbar width
    W = LW + RW + GW #total width
    DW_r = DW / ((CW+PW)/2)
    LW_r = LW / W
    RW_r = RW / W

    fig = figure(fid, figsize=(W, Hn), dpi=dpi, facecolor=facecolor)
    fig.clf()
    gs = gridspec.GridSpec(int(nrows) , ncols,
            width_ratios=[PW/(PW+CW), CW/(PW+CW)])
    if mapunits=='dB':
        data[data<=0]=nan #to prevent zeroes or negative numbers for logarithm

    gs.update(left=LW_r,right=(1. - RW_r),hspace=DH_r,wspace=DW_r, 
            bottom = TB_r, top = (1.-TB_r)  )

    pwrcolors=['b-','r-','g-','y-','c-','m-']

    if maptype[:3] in ['RTI','SNR'] and noise != None:
        #Add the noise lines in the first row
        ax = fig.add_subplot(gs[0])
        ylabel = 'Noise'
        if mapunits == 'dB':
            noise = 10*log10(noise)
            ylabel = 'Noise [dB]'
        Nmedianmax = Nmedianmin = median(noise[:int(finf.num_chan),:])
        for ch in range(int(finf.num_chan)):
            temparray = noise[chorder[ch],:].squeeze()
            Nmedian = median(temparray)
            if Nmedian > Nmedianmax:
                Nmedianmax = Nmedian
            if Nmedian < Nmedianmin:
                Nmedianmin = Nmedian
            pt = ax.plot(datearr, temparray,
                        pwrcolors[ch % size(pwrcolors)])
        prevfs = rcParams['font.size']
        rcParams['font.size'] = 4 #this is just to make label lines close
        leg = ax.legend(chname, bbox_to_anchor=(1.+ LegW/PW, 1.),
                loc='upper left')#,
#                borderaxespad=0.)#,
#                labelspacing=-0.3, borderpad = 0.05,
#                handlelength = 1.25, handletextpad = 0.1)
        rcParams['font.size'] = prevfs
        setp(leg.get_texts(), fontsize = fontsize)
        if mapunits == 'dB':
#            Nrange = [floor(Nchmin-3), ceil(Nminmax + 10)] # accounting galaxy
            Nrange = [floor(Nmedianmin-5), ceil(Nmedianmax + 10)] # galaxy
        else:
#            Nrange = [Nchmin*0.9, Nminmax*11] # accounting for galaxy
            Nrange = [Nmedianmin*0.8, Nmedianmax*11.] # accounting for galaxy
        ax.set_ylabel(ylabel, fontsize = fontsize)
        ax.set_ylim(Nrange)
        ax.tick_params(labelsize = fontsize)
        axR = ax.twinx() #ticks on the right
        axR.tick_params(labelsize = fontsize)
        axR.set_ylim(Nrange)
        axT = ax.twiny() #ticks on the right
        axT.tick_params(labelsize = fontsize)
        ax.grid()
        if acqtime_arr != None:
            ax.xaxis_date()
            ax.xaxis.set_major_formatter(dateFmt)
            axR.xaxis_date()
            axR.xaxis.set_major_formatter(dateFmt)
            axT.xaxis_date()
            axT.xaxis.set_major_formatter(dateFmt)
        ax.set_xticklabels([])
        axR.set_xticklabels([])
#        axT.set_xticklabels([])
        ax.set_xlim(txrange)
        axT.set_xlim(txrange)
        axR.set_xlim(txrange)

    for ch in chrowindexes:
        ax = fig.add_subplot(gs[(ch+noiseindex)*ncols])
        if mapunits=='dB':
            temparray = 10*log10(data[chorder[ch],:,:])
        elif maptype[:3]=='Coh':
            temparray = (data[ch,:,:])
        else:
            temparray = (data[chorder[ch],:,:])
        pcm = ax.pcolormesh(datearr,hrange,temparray,
                vmin=vmin[ch],vmax=vmax[ch],
                cmap = cmap)
        cbar = fig.colorbar(pcm,
                cax=axes(gs[(ch+noiseindex)*ncols+1].get_position(fig)))
        cbar.ax.tick_params(labelsize = fontsize)
        if vectormap!= None:
            tfactor = 3
            hfactor = 3
            U = vectormap[0,::hfactor,::tfactor]
            V = vectormap[1,::hfactor,::tfactor]
            U[abs(U)>300]=nan
            V[abs(V)>80]=nan
            U = U / 5
            Q = ax.quiver(datearr[::tfactor],hrange[::hfactor], U, V)
        if maptype[:3]=='Coh':
            bname = 'pair (%s,%s)'%tuple(array(chname)[finf.pairs[ch,:]])
        else:
            bname = chname[chorder[ch]]
        axR = ax.twinx() #ticks on the right
        axR.tick_params(labelsize = fontsize)
        axR.set_ylim(hlim)
        axT = ax.twiny() #ticks on the right
        axT.tick_params(labelsize = fontsize)
        if not acqtime_arr == None:
            axT.xaxis_date()
            axT.xaxis.set_major_formatter(dateFmt)
            ax.xaxis_date()
            ax.xaxis.set_major_formatter(dateFmt)
        if ch == chrowindexes[0]:
            cbar.ax.set_title('%s'%mapunits,fontsize = fontsize)
            if noise != None:
                axT.set_xticklabels([])
        else:
            axT.set_xticklabels([])
        ax.set_ylabel('%s\nHeight [km]'%bname, fontsize = fontsize)
        ax.tick_params(labelsize = fontsize)

        if ch == chrowindexes[-1]:
            ax.set_xlabel('Local Time', fontsize = fontsize)
        else:
            ax.set_xticklabels([])
            axR.set_xticklabels([])
        ax.set_xlim(txrange)
        ax.set_ylim(hlim)
        axT.set_xlim(txrange)
        ax.grid()
    fig.suptitle('%s map, %s'%(maptype,sdate), fontsize=fontsize*1.2)

    if drawinside:
        draw()
        fig.show()
    return fig
