from obspy.core import read
import glob
import matplotlib.pyplot as plt
from scipy.misc import electrocardiogram
from scipy.signal import find_peaks
import numpy as np
import math  
from time import process_time
import ntpath
st = process_time()
from rf import RFStream, read_rf, IterMultipleComponents, get_profile_boxes
plt.style.use('ggplot')


def calc_h_kappa(vp = 6.3,p = 0.06,w1=0.75,w2 = 0.25,outfile = "h-kappa-values.txt",data_dir_loc = "../results/dataRF", outloc="./"):
    # vp = 6.3
    # p = 0.06 #ray parameter
    # w1=0.75 # Weight Ps to 0.75
    # w2 = 0.25 # Weight PpPs to 0.25
    # outfile = 'h-kappa-values.txt'
    f= open(outloc+outfile,'w')
    # data_dir_loc = "../results/dataRF"
    data_files = glob.glob(data_dir_loc+"/*-rf_profile_rfs.h5")
    # print(data_files)
    for data in data_files:

        # print(f"data file is {data}")
        network = ntpath.basename(data).split('-')[0]
        station = ntpath.basename(data).split('-')[1]   
        st = read_rf(data)

        st = st.select(component="L")
    #    st = st.trim2(-25, 75, 'onset')
        len_trace_list=[]
        for tr in st:
            lentr=tr.stats.npts
            len_trace_list.append(lentr)
        if len(set(len_trace_list))> 1:
            # print(f"Unequal trace lengths for {data}")
            continue
        st = st.stack()
        for index,trace in enumerate(st):

            errorphase = False
            nbphase = 0
            [xpeaks, ypeaks] = [],[]
            trace.filter('bandpass', freqmin=0.005, freqmax=2)
            t = trace.stats.starttime
            pps = trace.stats.sampling_rate
    #        print(trace.stats)
            trace.trim(t+24, t+44)
            xpeaks, ypeaks = find_peaks(trace, height=0.02, distance=50)

            if len(xpeaks) > 2:
                if len(xpeaks) < 5:
                    print('nb of peaks =',len(xpeaks))
                    plt.plot(trace)
                    plt.plot(xpeaks, trace[xpeaks], "x")
                    plt.plot(np.zeros_like(trace), "--", color="gray")
    #                plt.show()
                    t0 = xpeaks[0]/pps
                    t1 = xpeaks[1]/pps
                    t2 = xpeaks[2]/pps
                    if len(xpeaks) > 3:
                        t3 = xpeaks[3]/pps
                    if t0 < 2.5 and t0 > 0:
                        t0 = t0
                    else:
                        t0 = np.NaN
                        errorphase = True
                    if t1 < 7.0 and t1 > 2.6:
                        t1 = t1
                    else:
                        t1 = np.NaN
                        errorphase = True
                    if t2 < 14.0 and t2 > 7.5:
                        t2 = t2
                    else:
                        if t3 and t3 < 14.0 and t2 > 7.5:
                            t2 = t3
                        else:
                            t2 = np.NaN
                            errorphase = True
    #                print(t0,t1,t2)
                
                    try:
                        if w1+w2 != 1:
                            raise ValueError('Weights are not properly defined')
                    except ValueError as e:
                        exit(str(e))

    # Measure the difference between theory and data:
                    if not errorphase:
                        numpoints = 1000
                        hs = np.linspace(20,40,numpoints)
                        Kappas = np.linspace(1.5, 2.5, numpoints)
                        H, K = np.meshgrid(hs, Kappas)
                        depth1 = (t1-t0)/(np.sqrt((K/vp)**2-(p)**2)-np.sqrt((1/vp)**2-(p)**2))
                        depth2 = (t2-t0)/(np.sqrt((K/vp)**2-(p)**2)+np.sqrt((1/vp)**2-(p)**2))
                        deltas = np.absolute((w1*depth1 + w2* depth2) - H)

    ## VISUALIZATION
                        # fig, ax = plt.subplots(2,2,figsize=(8,6),gridspec_kw={"width_ratios":[1, 0.05]})
                        fig = plt.figure()
                        
                        delta_lvs = np.linspace(np.amin(deltas),np.amax(deltas),30)


                        fig, axes = plt.subplots(nrows=2, ncols=1)
                        cmap = plt.get_cmap('rainbow_r')
                        CS = axes[0].contourf(H, K, deltas, levels=delta_lvs,cmap=cmap)
                        result = np.where(deltas == np.amin(deltas))
        # print(H[result], K[result], deltas[result])
                        axes[0].plot(H[result],K[result],'ko')
                        f.write(f"{network},{station},{trace.stats.station_latitude:.4f},{trace.stats.station_longitude:.4f},{H[result][0]:.2f},{K[result][0]:.2f}\n")
                        # axes[0].clabel(CS, inline=1, fontsize=10, fmt='%2.1f', colors='w')
                        axes[0].set_title(r'$H$-$\kappa$ grid search')
                        axes[0].set_xlabel('H')
                        axes[0].set_ylabel(r'$\kappa$')


                        # ax2 = plt.subplot2grid((2, 3), (0, 2), colspan=1)
                        # plt.colorbar(CS,ax=ax2)

                        times_data = np.arange(0,len(trace.data))/pps
                        axes[1].plot(times_data,trace.data)
                        axes[1].plot(xpeaks/pps, trace[xpeaks], "x")
                        #print(trace.data[np.where(times_data == t0)])
                        #axes[1].plot(np.zeros_like(trace), "--", color="gray")
                        axes[1].annotate('P', (t0+0.2, trace.data[np.where(times_data == t0)]), textcoords='data', size=10)
                        axes[1].annotate('PS', (t1, trace.data[np.where(times_data == t1)]), textcoords='data', size=10)
                        axes[1].annotate('PpPs', (t2, trace.data[np.where(times_data == t2)]), textcoords='data', size=10)
                        plt.tight_layout()

                        fig.subplots_adjust(right=0.82)
                        cbar_ax = fig.add_axes([0.85, 0.56, 0.03, 0.36])
                        fig.colorbar(CS, cax=cbar_ax)

                        # fig.tight_layout()

                        plt.savefig(outloc+f'H-K/{network}-{station}-h-k_outfile-{index}.png')
                    else:
                        print('bad peaks')
    f.close()


if __name__=="__main__":
    calc_h_kappa()