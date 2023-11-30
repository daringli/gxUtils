#!/usr/bin/env python

from netcdf_util import netcdf_file


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

import sys

from scipy.signal import argrelmax, argrelmin

from eik_tools import EikFile as GeomFile

def remove_shallow_wells(y, imax, imin, well_threshold=0.05):
    filtered_imax = list(imax)
    filtered_imin = list(imin)
    # You need to be this deep to be considered a minimum
    yrange = np.max(y) - np.min(y)
    i = 0 
    while i < len(filtered_imax) - 1:
        Dy1 = y[filtered_imax[i]] - y[filtered_imin[i]]
        Dy2 = y[filtered_imax[i+1]] - y[filtered_imin[i]]
        # smallest Dy sets the depth of well
        if Dy1 > Dy2:
            Dy = Dy2
            to_del = i + 1
        else:
            Dy = Dy1
            to_del = i

        #DBs.append(DB)
        if Dy/yrange > well_threshold:
            #filtered_imax.append(maxmin[i][0])
            i = i + 1
        else:
            del filtered_imax[to_del]
            del filtered_imin[i]
            
    imin = np.array(filtered_imin)
    imax = np.array(filtered_imax)
    #imax, imin = remove_duplicate_extrema(y,imax,filtered_imin)
    return imax, imin



def remove_duplicate_extrema(y, imax, imin):
    imax = np.array(imax)
    imin = np.array(imin)
    
    Nmax = len(imax)
    Nmin = len(imin)
    if Nmax  == Nmin + 1:     
        mask1 = imax[:-1] < imin
        mask2 = imax[1:] > imin
        mask = np.logical_and(mask1,mask2)
        if not(mask.all()):
            do_careful = True
            #print(mask1)
            #print(mask2)
            print("Masks don't agree")
        else:
            filtered_imax = imax
            filtered_imin = imin
            do_careful = False
    else:
        do_careful = True

    #do_careful = True
    if do_careful:
        # this should always work, but is slower
        print("Warning! Something odd with the maxes and mins. Starting careful well construction")

        filtered_imax = list(imax)
        filtered_imin = []
        i = 0
        while i < len(filtered_imax)-1:
            i0 = filtered_imax[i]
            i1 = filtered_imax[i+1]
            #print(i1)
            #print(imin)
            indices_between = np.logical_and(imin<i1, imin>i0)
            imins_between = imin[indices_between]
            Nimins_between = len(imins_between)
            #pdb.set_trace()
            if Nimins_between > 1:
                # take lowest minimum
                # and proceed to next maximum
                _tmp = np.argmin(y[imins_between])
                filtered_imin.append(imins_between[_tmp])
                i = i + 1
            elif Nimins_between == 0:
                # there is no minimum between the maximum
                # remove the smallest of the two maximum
                # and repeat
                if y[filtered_imax[i]] > y[filtered_imax[i+1]]:
                    del filtered_imax[i+1]
                else:
                    del filtered_imax[i]
            else:
                # just 1 minimum between
                filtered_imin.append(imins_between[0])
                i = i + 1
                
        filtered_imax = np.array(filtered_imax)
        filtered_imin = np.array(filtered_imin)
            
    return filtered_imax, filtered_imin


class ClickHandler(object):
    def __init__(self, eik_input = 'gx.eik.out', bg = None):
        #self.imax_pressed = False
        #self.imin_pressed = False
        self.maxlines = []
        self.minlines = []

        self.bg = bg
        self.imin = None
        self.imax = None

        self.gf = GeomFile(eik_input)
        self.outputs = ["bmag", "gradpar", "gbdrift", "cvdrift", "cvdrift0", "gds2", "gds21", "gds22"]
        self.Noutputs = len(self.outputs)
        self.has_extremas = [False] * self.Noutputs
        self.imaxs = [None] * self.Noutputs
        self.imins = [None] * self.Noutputs

        self.true_maxlines = [ [] for _ in range(self.Noutputs) ]
        self.true_minlines = [ [] for _ in range(self.Noutputs) ]
        self.distance_lines = [ [] for _ in range(self.Noutputs) ]
        self.distance_texts  = [ [] for _ in range(self.Noutputs) ]
        
        self.last_x = None
        self.L = np.nan


    def set_background(self):
        image = plt.imread(self.bg)
        aspect = 690/194.
        height = 0.1
        bot = 1 - height
        width = height * aspect
        left = 0.5 - width/2
        #left = 0.25
        #bot = 0.75
        background_ax = plt.axes([left, bot, width, height]) # create a dummy subplot for the background
        background_ax.set_zorder(-1) # set the background subplot behind the others
        background_ax.imshow(image, aspect='auto') # show the backgroud image
        background_ax.get_xaxis().set_visible(False)
        background_ax.get_yaxis().set_visible(False)
        background_ax.set_axis_off()
        
    def find_avg(self):
        for i,y in enumerate(self.ys):
            self.y_avg[i] = np.mean(y)
        
    def find_extrema(self, i):
        if not self.has_extremas[i]:
            y = self.ys[i]

            imax = np.sort(argrelmax(y, order=5,mode='wrap')[0])
            imin = np.sort(argrelmin(y, order=5,mode='wrap')[0])
            #imax, imin = remove_duplicate_extrema(y, imax, imin)
            imax, imin = remove_shallow_wells(y, imax, imin)

            #iextrema = np.concatenate((imax,imin))
            #iextrema.sort()

            imax.sort()
            imin.sort()
            
            self.has_extremas[i] = True
            self.imaxs[i] = imax
            self.imins[i] = imin
        
        
    def plot(self):
        Noutputs = self.Noutputs
        ntheta = self.gf.ntheta + 1
        
        self.ys = np.zeros((Noutputs,ntheta))
        self.y_avg = np.zeros(Noutputs)
        # TODO: not sure if scaled theta or theta is most appropriate. depends on gradpar definition. look up.
        self.x = self.gf.scaled_theta
        self.fig, axes = plt.subplots(4,2, sharex=True)
        self.axes = axes.flatten()
        for i, output in enumerate(self.outputs):
            self.ys[i] = getattr(self.gf,output)
            self.axes[i].plot(self.x,self.ys[i])
            self.axes[i].set_title(output)
        
        self.find_avg()

            
        self.axes[-2].set_xlabel(r"$\theta$")
        self.axes[-2].set_xlabel(r"$\theta$")
        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        if self.bg is not None:
            self.set_background()
        plt.title(self.gf.eikfile)
        plt.show()


    def onclick(self,event):
        if event.inaxes is not None:
            event_ax = event.inaxes
        else:
            print("Must press a subplot (a matplotlib axis)")
            return

        event_i = None
        for i,ax in enumerate(self.axes):
            if event_ax == ax:
                event_i = i
                break
        else:
            print("Axis not recognized??")
            return

        self.find_extrema(i)
        min_pressed = False
        max_pressed = False
        
        event_x = event.xdata
        if str(event.button) == "MouseButton.LEFT":
            color = 'b'
            #minline = self.axes[event_i].axvline(event_x, color=color)
            index_of_minimum  = (np.abs(self.x[self.imins[event_i]]-event_x)).argmin()
            self.imin = self.imins[event_i][index_of_minimum]
            min_pressed = True
            
        elif str(event.button) == "MouseButton.RIGHT":
            color = 'r'
            #minline = self.axes[event_i].axvline(event_x, color=color)
            #print(self.imaxs[event_i])
            index_of_maximum  = (np.abs(self.x[self.imaxs[event_i]]-event_x)).argmin()
            self.imin = self.imaxs[event_i][index_of_maximum]
            max_pressed = True

       
            

        if min_pressed or max_pressed:
            self.this_x = self.x[self.imin]
            if self.last_x is not None:
                # we've got two points
                # can calculate distance
                self.L = np.abs(self.this_x - self.last_x)
                print("Distance: " + str(self.L))
                
            
            # update all subplots
            tmps = [None] * self.Noutputs
            for j in range(self.Noutputs):
                #print(j)
                tmps[j] = self.axes[j].axvline(self.x[self.imin], linestyle='dashed', color=color)
                self.true_minlines[j].append(tmps[j])
                #print(self.true_minlines[j])
                if len(self.true_minlines[j]) > 2:
                    self.true_minlines[j][0].remove()
                    del self.true_minlines[j][0]
                if not np.isnan(self.L):
                    l = self.axes[j].plot([self.this_x, self.last_x], [self.y_avg[j], self.y_avg[j]], color='k')
                    self.distance_lines[j].append(l[0])
                    t_x = 0.5*self.this_x +  0.5*self.last_x
                    t_y = self.y_avg[j] * 1.01
                    t = self.axes[j].text(t_x, t_y, "{:.2f}".format(self.L), )
                    self.distance_texts[j].append(t)
                    if len(self.distance_lines[j]) > 1:
                        self.distance_lines[j][0].remove()
                        del self.distance_lines[j][0]
                        self.distance_texts[j][0].remove()
                        del self.distance_texts[j][0]

            self.last_x = self.this_x
            plt.draw()
        #self.set_background()

        

if __name__=="__main__":
    import sys
    from glob import glob
    
    if len(sys.argv) > 1:
        dir = sys.argv[1]
    else:
        dir = '.'

    geofile = sorted(glob(dir + "/gx_*_geo.nc"))[-1]
    ch = ClickHandler(geofile)
    ch.plot()
