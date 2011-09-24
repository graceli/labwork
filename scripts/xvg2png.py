#from scitools.std import *
import numpy as np

def xvg2array(xvgf):
    f1 = open(xvgf, 'r')
    xlist = []
    ylist = []
    for line in f1:
        if line.startswith('#') or line.startswith('@'):
            pass
        else:
            linelist = line.split()
            if len(linelist) >= 2:
                xlist.append(eval(linelist[0]))
                ylist.append(eval(linelist[1]))
    xarray = np.array(xlist,dtype=float)
    yarray = np.array(ylist,dtype=float)
    f1.close()
    return xarray, yarray

def xvg_average(xvgf):
    yarray = xvg2array(xvgf)[1]
    ave = average(yarray)
    return ave    
    

# def xvg2png(xvgf, outputfile):
#     arraylist = xvg2array(xvgf)
#     xarray = arraylist[0]
#     yarray = arraylist[1]
#     plot(xarray,yarray,legend=xvgf,show=False,hardcopy=outputfile)

# def twoxvg2png(xvgf1,xvgf2,outputfile):
#     """ process the xvg file to obtain array data"""

#     arraylist1 = xvg2array(xvgf1)
#     xarray1 = arraylist1[0]
#     yarray1 = arraylist1[1]

#     arraylist2 = xvg2array(xvgf2)
#     xarray2 = arraylist2[0]
#     yarray2 = arraylist2[1]
    
#     plot(xarray1,yarray1,'b',legend=xvgf1,show=False)
#     hold('on')
#     plot(xarray2,yarray2,'r',legend=xvgf2,show=False)
#     hold('off')
#     hardcopy(outputfile)
#     print 'done'

# def threexvg2png(xvgf1,xvgf2,xvgf3,outputfile):
    
#     arraylist1 = xvg2array(xvgf1)
#     xarray1 = arraylist1[0]
#     yarray1 = arraylist1[1]

#     arraylist2 = xvg2array(xvgf2)
#     xarray2 = arraylist2[0]
#     yarray2 = arraylist2[1]

#     arraylist3 = xvg2array(xvgf3)
#     xarray3 = arraylist3[0]
#     yarray3 = arraylist3[1]

#     plot(xarray1,yarray1,'g',legend=xvgf1,show=False)
#     hold('on')
#     plot(xarray2,yarray2,'b',legend=xvgf2,show=False)
#     hold('on')
#     plot(xarray3,yarray3,'r',legend=xvgf3,show=False)
#     hold('off')
#     hardcopy(outputfile)
#     print 'done'

# if __name__ == '__main__':
#     twoxvg2png('e2et100w300p0.xvg', '../m300/e2et100m300p0.xvg', 'e2e300.png')
#     twoxvg2png('gyrt100w300p0.xvg', '../m300/gyrt100m300p0.xvg', 'gyr300.png')
#     twoxvg2png('hbppt100w300p0.xvg', '../m300/hbppt100m300p0.xvg', 'hbpp300.png')
#     twoxvg2png('hbpst100w300p0.xvg', '../m300/hbpst100m300p0.xvg', 'hbps300.png')

#     threexvg2png('./v100/e2ev100.xvg', './w100/e2et100w100p0.xvg', './m100/e2et100m100p0.xvg', 'e2e100.png')
#     threexvg2png('./v100/gyrv100.xvg', './w100/gyrt100w100p0.xvg', './m100/gyrt100m100p0.xvg', 'gyr100.png')
#     threexvg2png('./v100/hbnumv100.xvg', './w100/hbppt100w100p0.xvg', './m100/hbppt100m100p0.xvg', 'hbpp100.png')
#     twoxvg2png('./w100/hbpst100w100p0.xvg', './m100/hbpst100m100p0.xvg', 'hbps100.png')
