# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 10:17:25 2017

@author: taizo kawano
"""
"""
###########################################################################

Copyright (c) 2018, Taizo kawano


# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

###########################################################################
"""
import os
import sys
#import time
#import datetime
#import csv
import pandas
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#from scipy import stats

########################################
#dot plot function
#median with inter quatile range are shown as bars
#_thedata: list of list: [[],[]...]. not numpy array.
# ylim: (min, max)
# or melted form of pandas.Dataframe. 
#if melted form, label must come first column
#**kwargs could be labels, groupnames ylim, size,thickness,
# binnum, coeff, sort, figsize
#191227 mct measures of central tendency = mean average and sd
# also errorbar = "sd" or "se"
#211118 separete getposfordot getsift function 
#to utilize other type of plots (two way anova)
#231020 spyder5.4.3 inline plot need yticks([]) etc. to avoid doubly ticks
#####################################
def dotplots(_thedata, **kwargs):
    colorlist = ["black","orangered", "deepskyblue", "lime",
                 "blueviolet", "green","deeppink","orange"]
    markerlist = [".","2", "x", "s",
                 "*", "d","1","p"]
    if "col" in kwargs:
        col = kwargs["col"]
    else:
        col="black"
    if "figsize" in kwargs:
        figsize = kwargs["figsize"]        
    else:
        figsize = None        
    fig = plt.figure(figsize = figsize)
    #231020 these lines needed to avoid doubly labeld axis
    #but disappear y-tick on python 3.6
    #plt.yticks([])
    plt.xticks([])
    if sys.version_info.minor > 6 and sys.version_info.major == 3:
        plt.yticks([])        
        #on python 3.6 these delete spines?
        plt.gca().spines['left'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)

    
    #plt.gca().yaxis.set_ticks_position('left')
    #plt.gca().xaxis.set_ticks_position('bottom')
    #plt.tick_params(labelbottom='off')
    datalist = []
    if "labels" in kwargs:
        labels = kwargs["labels"]        
    else:
        labels = None        
    if "groupnames" in kwargs:
        templabels = kwargs["groupnames"]        
    else:
        templabels = None        
    datalist = []
    #incase that _thedata is melted dataframe
    if type(_thedata) == pandas.core.frame.DataFrame:
        if templabels is None:
            templabels = _thedata[_thedata.columns[0]].unique()
        if labels is None:
            labels = templabels
        groupnumber = len(templabels)
        print("groupnumber "+str(groupnumber))
        for x in templabels:            
            sub = list(_thedata[_thedata.columns[1]][_thedata[_thedata.columns[0]]==x])
            datalist.append(sub)
    #if the _thedata is list
    else:
        datalist = _thedata
        groupnumber = len(_thedata)        
        if labels is None:
            labels = [str(x) for x in np.arange(groupnumber)+1]
            
    print("labels " + str(labels))
    #return datalist
    if "ylim" in kwargs:
        ylim = kwargs["ylim"]        
    else:
        #20171205 not work yet. need fix
        ylim = (min([min(a) for a in datalist]), 
                max([max(a) for a in datalist]))        
    ax = fig.add_subplot(1,1,1)
    _therange = ylim[1]-ylim[0]
    ax.set_ylim(ylim[0]-_therange*0.05, ylim[1]+_therange*0.05)
    ax.set_xlim(-0.5, groupnumber-1+0.5)
    
    if "size" in kwargs:
        size = kwargs["size"]        
    else:
        size = 1        
    if "thickness" in kwargs:
        thickness = kwargs["thickness"]        
    else:
        thickness = 1        
    if "sort" in kwargs:
        sort = kwargs["sort"]        
    else:
        sort = False        
    if "outlier" in kwargs:
        outlier = kwargs["outlier"] 
    else:
        outlier = "none" 
    #xticklabels rotation
    if "rotation" in kwargs:
        rotation = kwargs["rotation"] 
        #horizontalaligment = "right"
        horizontalaligment = "center"
    else:
        rotation = 0
        horizontalaligment = "center"
    if "mct" in kwargs:
        mct = kwargs["mct"]
        print("mct specify")
    else:
        mct = "median"

    ax.set_xticks(range(groupnumber))
    ax.set_xticklabels(labels)
    samplesize = [len(x) for x in datalist]
    print("samplesize "+str(samplesize))
    labelswithsize = [labels[x]+"\n"+str(samplesize[x]) for x in range(len(labels))]
    ax.set_xticklabels(labelswithsize, rotation = rotation, ha = horizontalaligment)

    for i in range(groupnumber):
    #for adata in _thedata:
        #ydata = _thedata
        adata = np.array(datalist[i])
        #ydata = analyzeddataframe["meanfoq"]
        #ydata = ltqadf[0]["qduration"]
        themin = adata.min()
        themax = adata.max()
        
        #divide data into 4 so that each contains same number of sample
        q25, q75 = np.percentile(adata, [25, 75])
        
        if mct == "median":        
            #median
            themct = np.median(adata)
            #inter quatile range. 
            upperrange = q75
            lowerrange = q25
        elif mct == "mean":
            #plot mean
            #print("mean")
            themct = np.mean(adata)
            if "errorbar" in kwargs:
                if kwargs["errorbar"] == "sd":
                    upperrange = themct + np.std(adata)
                    lowerrange = themct - np.std(adata)
                elif kwargs["errorbar"] == "se":
                    upperrange = themct + np.std(adata)/np.sqrt(len(adata))
                    lowerrange = themct - np.std(adata)/np.sqrt(len(adata))
        print("mct = " + mct)

        #mct measures of central tendency
        ax.hlines(y = themct, xmin = i-0.3, xmax = i+0.3, 
                  colors = "black", linewidths = thickness)
        #plot inter quatile or sd range. 
        ax.hlines(y = [lowerrange,upperrange], xmin = i-0.1, xmax = i+0.1, 
                  colors = "black", linewidths = thickness)
        ax.vlines(x=i, ymin = lowerrange, ymax=upperrange, 
                  colors = "black", linewidths = thickness)
        
        
        #histdata = np.histogram(adata, bins=np.linspace(themin, themax, num=25))
        if "binnum" in kwargs:
            binnum = kwargs["binnum"]        
        else:
            binnum = 10        
        binslist=np.linspace(themin, themax, num=binnum)
        kwargs["bins"] = binslist
        #outlier handling
        if outlier != "none":
            upperlimit = q75+1.5*(q75-q25)
            lowerlimit = q25-1.5*(q75-q25)
            if outlier == "line":
                ax.vlines(x=i, ymin = lowerlimit, ymax=upperlimit, 
                      colors = "gray", linewidths = thickness/2)
            #outlier are plotted first
            if outlier == "oc":#open circle
                boolarray = np.logical_or(np.array(adata) > upperlimit, np.array(adata) < lowerlimit)
                outdata = adata[boolarray]
                ax.scatter([i for a in range(len(outdata))],outdata,
                           s = size*2, color="black",facecolors="none",
                           marker = "o")
                pdata = adata[~boolarray]
        else:
            pdata = adata
        #20171129 use kwargs to specify coeff
        if "coeff" in kwargs:
            coeff = kwargs["coeff"]        
        else:
            #coeff = 1/histdata[0].max()**2*(np.sqrt(size)/3)
            coeff = 0.5
            
        """
        print("binslist {}".format(binslist))
        histdata = np.histogram(pdata, bins = binslist)
        #print(histdata)
        #numpy.digitize Return the indices of the bins to 
        #which each value in input array belongs. this might better?
        #max number as shuring along x axis
        #if many data is exist use this
        #coeff = histdata[0].max()*1.5
        #only a few, use this
        #coeff = 1/size*5
        #20170914 changed to shiftarrays[j]*coeff, so try this way
        #coeff = 1/histdata[0].max()**2*(np.sqrt(size)/3)
        #20171129 use kwargs to specify coeff
        if "coeff" in kwargs:
            coeff = kwargs["coeff"]        
        else:
            #coeff = 1/histdata[0].max()**2*(np.sqrt(size)/3)
            coeff = 0.5

        #array for how far shift from the middle line
        def getxshiftarray(anarray):
            returnlist = []
            for a in anarray:
                shiftarray = np.array([int((i+1)/2) if i%2 ==0 else -1*int((i+1)/2) for i in range(a)])
                if a%2 ==0:
                    shiftarray=shiftarray+0.5
                returnlist.append(shiftarray)
            return returnlist
        
        shiftarrays = getxshiftarray(histdata[0])
        maxnuminbin = histdata[0].max()
        for j in range(len(histdata[1])-1):
            subsetdata = pdata[(histdata[1][j] <= pdata) & (pdata < histdata[1][j+1])]
            #the right most bin is include max.
            if j == len(histdata[1])-2:
                subsetdata = pdata[(histdata[1][j] <= pdata) & (pdata <= histdata[1][j+1])]
            if len(subsetdata)>0:
                if sort:
                    subsetdata = np.sort(subsetdata)
                #print(str(j) +" "+ str(subsetdata))
                #ax.scatter(shiftarrays[j]/coeff+i,subsetdata,s = size ,color="black")
                #ax.scatter(shiftarrays[j]*coeff+i,subsetdata,s = size ,color="black")
                ax.scatter(shiftarrays[j]/maxnuminbin*coeff+i,subsetdata,
                           s = size ,color="black")
                #ax.scatter(shiftarrays[j]/maxnuminbin*coeff+i,
                #           np.sort(subsetdata),s = size ,color="black")
        """
        #211118 using def getposfordot(_data, **kwargs) get x and y array
        #actuall dot ploting part
        xyarray = getposfordot(pdata, **kwargs)
        ax.scatter(xyarray[0]*coeff+i, xyarray[1], color=col,
                   s=size)
        #"""
    #ax.set_xticks(hrtick*60*60/interval)
    #ax.set_xticklabels(hrtick)    
    #ax.set_ylim(ylim[0]-_therange*0.05, ylim[1]+_therange*0.05)
    #ax.set_xlim(-0.5, groupnumber-1+0.5)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    fig.tight_layout()
    #print("test231020")

    return fig
























#211118 just copy pased from lethargus figure script. need to generelize.
#_df 1st column treatment 2nd genotype
def dotplotfor2wanova(_df, _param, **kwargs):
        
    colorlist = ["black","orangered", "deepskyblue", "lime",
                 "blueviolet", "green","deeppink","orange"]
    markerlist = [".","2", "x", "s",
                 "*", "d","1","p"]
    if "figsize" in kwargs:
        figsz = kwargs["figsize"]        
    else:
        figsz = (3,4)       
    
    #xpos = [1,2]
    levels0 = _df.iloc[:,0].unique()
    levels1 = _df.iloc[:,1].unique()
    
    
    fig = plt.figure(figsize = figsz)
    ax = fig.add_subplot(1,1,1)
    # 
    for i,(l1, col) in enumerate(zip(levels1, colorlist)):
        # 
        for j, (l0, ms) in enumerate(zip(levels0, markerlist)):
            gdf = _df[(_df.iloc[:,0]==l0) & (_df.iloc[:,1]==l1)]
            xyarray = getposfordot(gdf[_param])
            xcen = l0+i*0.2
            #x_ = [l0 for a in range(len(gdf))]
            #ax.scatter(xyarray[0]*0.5+l0, xyarray[1], color=col, marker = ms)
            
            ax.scatter(xyarray[0]*0.2+xcen, xyarray[1], color=col, marker = ".")
    
            #mean
            mean_a = np.mean(gdf[_param])
            ax.plot([xcen-0.15,xcen+0.15],[mean_a,mean_a],
                    color=col, linewidth = 1)
            #std error
            """
            sd_a = np.std(gdf.readout)
            se_a = sd_a/np.sqrt(len(gdf.readout))
            ax.scatter(l0, mean_a, s=4, color=col)
            ax.plot([l0,l0],[mean_a-se_a,mean_a+se_a],
                    color=col, linewidth = 1)
            """
        ax.annotate("{} {}".format(_df.columns[1],l1),
                    xy=(0.1, 0.9 - i*0.1),
                    color = col,
                    xycoords='axes fraction', fontsize=8,
                    horizontalalignment='left', verticalalignment='bottom')
    
    ax.set_xticks((0,1))
    ax.set_xticklabels([])
    ax.set_xlabel( _df.columns[0])
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, np.max(_df[_param])*1.05)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.tight_layout()
    return fig
        



#211118 implemented to dotplot library and used dotplots function 
#211022 for like a 2 factor dotplot eg. drag and genotype
#_data is array or pandas series?
def getposfordot(_data, **kwargs):
    sort = False
    if "sort" in kwargs:
        sort = kwargs["sort"]        
    if "bins" in kwargs:
        bin_ = kwargs["bins"]        
    else:
        bin_ = 10

    pdata = _data
    #print("binslist {}".format(bin_))

    histdata = np.histogram(pdata, bins = bin_)
    #print(histdata)
    #array for how far shift from the middle line
    def getxshiftarray(anarray):
        returnlist = []
        for a in anarray:
            shiftarray = np.array([int((i+1)/2) if i%2 ==0 else -1*int((i+1)/2) for i in range(a)])
            if a%2 ==0:
                shiftarray=shiftarray+0.5
            returnlist.append(shiftarray)
        return returnlist
    
    shiftarrays = getxshiftarray(histdata[0])    
    maxnuminbin = histdata[0].max()
    returnval = np.empty((2, len(pdata)))
    
    sti=0
    endi=0
    for j in range(len(histdata[1])-1):
        subsetdata = pdata[(histdata[1][j] <= pdata) & (pdata < histdata[1][j+1])]
        #the right most bin is include max.
        if j == len(histdata[1])-2:
            subsetdata = pdata[(histdata[1][j] <= pdata) & (pdata <= histdata[1][j+1])]
        endi=sti+len(subsetdata)
        #print(subsetdata)
        #print("s e {} {}".format(sti, endi))
        if len(subsetdata)>0:
            if sort:
                subsetdata = np.sort(subsetdata)
            #print("shiftarrays[j] "+str(shiftarrays[j]))
            #x
            returnval[0,sti:endi]=shiftarrays[j]/maxnuminbin
            #y
            returnval[1,sti:endi]=subsetdata
            sti=endi
    return returnval








#if __name__ == "__main__":
#    return 
    
"""
dotplotfig = dotplot.dotplots(tempdf,
                              groupnames = uniquegroupnames,
                              #mct = "mean",
                              #errorbar = "se",
                              mct = "median",
                              #errorbar = "se",
                              ylim = (0,max(tempdf[param])),
                              binnum = 10, size = 5, thickness = 1,
                              figsize = dotfigsz)
                              #figsize = (3,4))
"""





















