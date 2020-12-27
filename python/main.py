#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
reads ../data/olenina_Table_ICES_id_volpcell.csv,
calculates parameters like Volume, Area, Equivalent Spherical Area/Volume, Surface extension, Aspect ratio, etc..
and writes the extended data with additional parameters in ../data/olenina_Table_ICES_id_volpcell_geom.csv
! The script does not yet work for the datasets from Roselli et al., which employs a different format to describe dimensions !
Tested for: python 3.6
"""

#Built in / Generic imports
import os

#Library imports
import csv

#Local imports
#from calc_VA_funcs import *
from calc_geom_funcs import *

__author__ = "Onur Kerimoglu"
__email__  = "onur.kerimoglu@uol.de"
__credits__ = ["Onur Kerimoglu", "Alexey Ryabov"]
__version__ = "1.0.0" # December 2020
#__version__ = "0.1.0" October 2018 (named calculate_VA.py back then)

def main():
    cwd = os.getcwd()
    rootpath = os.path.dirname(cwd)

    dataset = 'Olenina' #options: Olenina, Roselli

    fin = os.path.join(rootpath, 'data', 'olenina_Table_ICES_id_volpcell.csv')
    fout = os.path.join(rootpath, 'data', 'olenina_Table_ICES_id_volpcell_geom.csv')

    #open the input file
    with open(fin, 'rU') as f:
        reader = csv.reader(f)
        Din = list(reader)
        headersin = Din[0]
        del Din[0:2]

    error_tolerance = 0.5
    headersnew = [u'Volume [um3]',u'Area [um2]',u'Eqv_Rad_v [um]',u'Eqv_Rad_s [um]',u'Surf_ext [-]',
                  u'Lmin [um]',u'Lmid [um]',u'Lmax [um]',u'DmMid2Min [-]',u'DmMax2Mid [-]',u'AspectRatio [-]',
                  u'Prolate [Y/N]',u'Elongation']
    #append headers
    headersout = headersin+headersnew

    #to compare the calculated volume with the provided volume
    Vcol = headersin.index('Average_Cellular_Volume um3')
    Dout = [None]*len(Din)
    errtol_crossed = False
    maxerr = 0
    shapedict = {'known':{}, 'unknown':{}}
    for r in range(len(Din)):
        V,A,gshape,EqRadV,EqRadS,SurfExt,LMin,LMid,LMax,DmMid2Min,DmMax2Mid,AspectR,Prolate,Elongation = get_geom_props(Din[r],headersin,dataset)

        #replace the possibly tuned gshape name
        Din[r] = correct_gshape(Din[r],gshape,headersin)
        #register the gshape either to known or unknown category
        cat = 'unknown' if np.isnan(V) else 'known'
        shapedict = register2shapedict(shapedict,cat, gshape)
        #check the calculated volume against the V provided in the data
        if (not np.isnan(V)) and (Din[r][Vcol] != '-'):
            Vin = np.float(Din[r][Vcol])
            err = abs(V-Vin)/Vin
            if err > error_tolerance:
                print ('! species id %s (%s): calculated V (%s) differs from the provided (%s) by more than %s%%'%(Din[r][0],gshape,V,Vin,100*error_tolerance))
                errtol_crossed = True
            maxerr=err if err>maxerr else maxerr
        Vstr = str(V) if not np.isnan(V) else 'NA'
        Astr = str(A) if not np.isnan(A) else 'NA'
        EqRadVs = str(EqRadV) if not np.isnan(EqRadV) else 'NA'
        EqRadSs = str(EqRadS) if not np.isnan(EqRadS) else 'NA'
        SurfExts = str(SurfExt) if not np.isnan(SurfExt) else 'NA'
        LMins = str(LMin) if not np.isnan(LMin) else 'NA'
        LMids = str(LMid) if not np.isnan(LMid) else 'NA'
        LMaxs = str(LMax) if not np.isnan(LMax) else 'NA'
        DmMid2Mins = str(DmMid2Min) if not np.isnan(DmMid2Min) else 'NA'
        DmMax2Mids = str(DmMax2Mid) if not np.isnan(DmMax2Mid) else 'NA'
        AspectRs = str(AspectR) if not np.isnan(AspectR) else 'NA'
        Prolates = str(Prolate) if not np.isnan(Prolate) else ''
        #Elongation = Elongation if not np.isnan() else 'NA'
        Dout[r]=Din[r]+[Vstr,Astr,EqRadVs,EqRadSs,SurfExts,LMins,LMids,LMaxs,DmMid2Mins,DmMax2Mids,AspectRs,Prolates,Elongation]
    with open(fout, 'w') as fp:
        writer = csv.writer(fp, delimiter=',')
        writer.writerows([headersout])
        writer.writerows(Dout)
    print ('data created:' + fout)
    print ('maximum error between calculated and provided volume: %5.2f%%'%(maxerr*100))
    if errtol_crossed:
        print ('see above for the instances where the error tolerance (%s%%) for the volume calculation was exceeded'%(100*error_tolerance))
    else:
        print ('error tolerance threshold (%s%%) was never exceeded'%(100*error_tolerance))
    plotshapedict(shapedict, fout.replace('.csv', '_shape_register'))

def get_geom_props(fields,headersin,dataset):
    gshape, V, A, Dsorted = calc_geom(fields, headersin, dataset)

    pi = np.pi

    #check if the V-A relationship is plausible, eg., A can't be smaller than the area of an equivalent sphere
    if (dataset=='Olenina') and (not np.isnan(A)):
        r=(V/ (4./3*pi))**(1./3) #radius of the eq. sphere
        As=4. * pi * r**2.
        reldif=(As-A)/As
        if reldif>0.01:
            print ('! For %s: A is implausible (lower than the eq. sphere). Nullifying Area'%gshape)
            A=np.nan

    # equvivalent spherical radius based on volume
    EqRadV = (V / (4 / 3 * pi)) ** (1 / 3.)
    # and equvivalent spherical radius based on surface area:
    EqRadS = (A / (4 * pi)) ** (1 / 2.)
    # surface extension
    SurfExt = A / (4 * pi * EqRadV ** 2)
    # minimal, middle and maximal linear cell dimensions
    # Dsorted=[1.,2.,3.]
    LMin = Dsorted[0]
    LMid = Dsorted[1]
    LMax = Dsorted[2]
    # aspect ratio
    DmMid2Min = LMid / LMin
    DmMax2Mid = LMax / LMid
    AspectR = LMax / LMin
    # Prolate or not
    Prolate = 0
    if DmMax2Mid < DmMid2Min:
        AspectR = 1. / AspectR  # for the non-prolates, AspectR=Lmin/Lmax
        Prolate = 0  # i.e, ind
    else:
        Prolate = 1  # i.e., ~ind
    # Classify: Prolate / Oblate / Compact / Unknown
    CompactCellRatio = 1.5
    if (DmMax2Mid < CompactCellRatio) & (DmMid2Min < CompactCellRatio):
        Elongation = 'C'
    elif DmMax2Mid >= DmMid2Min:
        Elongation = 'P'
    elif DmMax2Mid < DmMid2Min:
        Elongation = 'O'
    else:
        Elongation = 'U'

    return (V, A, gshape, EqRadV, EqRadS, SurfExt, LMin, LMid, LMax, DmMid2Min, DmMax2Mid, AspectR, Prolate, Elongation)

def correct_gshape(fields,gshape, headersin):
    shapecol = headersin.index('Geometric_shape')
    fields[shapecol] = gshape
    return (fields)

def register2shapedict(shapedict,cat,gshape):
    if gshape not in shapedict[cat].keys():
        shapedict[cat][gshape] = 1
    else:
        shapedict[cat][gshape]=shapedict[cat][gshape]+1
    return shapedict

def plotshapedict(shapedict,fname):
    import matplotlib.pyplot as plt
    knownshapes0=list(shapedict['known'].keys())
    knownvalues0=[shapedict['known'][shape] for shape in knownshapes0]
    knowncols=['b' for shape in knownshapes0]
    #sort, descending
    sortinds=np.flipud(np.argsort(knownvalues0))
    knownvalues=[knownvalues0[i] for i in sortinds]
    knownshapes=[knownshapes0[i] for i in sortinds]

    unknownshapes0 = list(shapedict['unknown'].keys())
    unknownvalues0 = [shapedict['unknown'][shape] for shape in unknownshapes0]
    unknowncols = ['r' for shape in unknownshapes0]
    # sort, descending
    sortinds = np.flipud(np.argsort(unknownvalues0))
    unknownvalues = [unknownvalues0[i] for i in sortinds]
    unknownshapes = [unknownshapes0[i] for i in sortinds]

    shapes=unknownshapes+knownshapes
    values=unknownvalues+knownvalues
    cols=unknowncols+knowncols
    ypos=np.arange(len(shapes))

    #fig
    fig=plt.figure(figsize=(8,8))
    ax=plt.subplot(111)
    plt.subplots_adjust(left=0.37,bottom=0.1,right=0.98,top=0.98)
    hall=plt.barh(ypos,values,align='center',color=cols)
    plt.yticks(ypos,shapes)
    plt.ylim(-0.5,ypos[-1]+0.5)
    plt.xlabel('# of species')
    plt.ylabel('Cellular Shape')
    h=[hall[0],hall[len(unknownshapes)]]
    plt.legend(np.flipud(h),['known shapes (%s)'%sum(knownvalues),'unknown shapes (%s)'%sum(unknownvalues)])
    #plt.show()
    plt.savefig(fname+'.png')
    print ('figure saved:'+fname+'.png\n')


if __name__=='__main__':
    main()