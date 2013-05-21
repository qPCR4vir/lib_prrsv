from numpy import *
from scipy.stats import pearsonr

def pearsonMat(a,names=[]):
    ax=1
    takeAx=0

    dim = a.shape[ax]
    rc=zeros((dim,dim),dtype=float)
    pairs=[]

    for i in range(dim):
        for j in range(i,dim):
            x=a[:,i]
            y=a[:,j]
            cond = (x>0)&(y>0)
            x=x[cond]
            y=y[cond]
            c,p = pearsonr(x,y)
            print i,j,c,p
            rc[i,j]=c
            rc[j,i]=c
            
            try:
                iName=names[i]
            except IndexError:
                iName=i
            try:
                jName=names[j]
            except IndexError:
                jName=j
            pairs.append((iName,jName,c))
    return rc,pairs


def dscatter(x,y,axis=None):
    """
    """
    from pylab import scatter,colorbar,axis
    from scipy.stats import pearsonr

    nBins=50
    maxValue = max((x.max(),y.max()))

    bins=logspace(0.0,ceil(log10(maxValue)),num=nBins)
    
    # remove zeros
    c = (x>0)&(y>0)
    c1 = (x>1)&(y>1)

    R=pearsonr(x[c],y[c])[0]
    R1=pearsonr(x[c1],y[c1])[0]
    
    r=column_stack((x[c],y[c]))
    d,dBins=histogramdd(r,bins=(bins,bins))

    pairs=set((tuple(t) for t in r.tolist()))
    pArray=array([t for t in pairs])
    ux=pArray[:,0]
    uy=pArray[:,1]
    d.shape
    def lookupDensity(x0,y0):
        i=dBins[0][dBins[0]<x0].shape[0]
        if i>=d.shape[0]:
            i=d.shape[0]-1
        j=dBins[1][dBins[1]<y0].shape[0]
        if j>=d.shape[1]:
            j=d.shape[1]-1
        #print i,j
        return d[i,j]

    dz=[]
    batches={}
    for i in xrange(len(ux)):
        z=lookupDensity(ux[i],uy[i]) +1 
        dz.append(z)
        if z not in batches:
            batches[z]=[]
        batches[z].append((ux[i],uy[i]))

    
    scatter(log10(ux),log10(uy),s=10,c=log10(dz),edgecolors='none')
    axis('scaled')
    l = floor(10*log10(min((ux.min(),uy.min()))))/10.0 - 0.1
    u =  ceil(10*log10(max((ux.max(),uy.max()))))/10.0 + 0.1
    axis([l,u,l,u])

    return R,R1

    #colorbar()
    
    #return ux,uy,dz,[l,u,l,u]

    
def semiLogHistLostFound(i,o):
    from pylab import histogram,semilogx

    f=i[(o>0) &(i>0)]
    l=i[(o<1) &(i>0)]

    d=i.sum()
    fNorm=f/d
    lNorm=l/d
    
    fHist=histogram(fNorm,logspace(-6,-2.0,30))
    lHist=histogram(lNorm,logspace(-6,-2.0,30))
    
    semilogx(fHist[1][1:],fHist[0],label='out>0')
    semilogx(lHist[1][1:],lHist[0],label='out=0')

def logLogRatioFoundLost(i,o,**pltKwds):
    
    from pylab import histogram,semilogx,loglog

    f=i[(o>0) &(i>0)]
    l=i[(o<1) &(i>0)]

    d=i.sum()
    fNorm=f/d
    lNorm=l/d
    
    fHist=histogram(fNorm,logspace(-6,-2.0,30))
    lHist=histogram(lNorm,logspace(-6,-2.0,30))
    
    #semilogx(fHist[1][1:],fHist[0].astype(float)/lHist[0])
    loglog(fHist[1][1:],fHist[0].astype(float)/lHist[0],**pltKwds)
                 
                    
def semiLogFracFound(i,o,**pltKwds):
    from pylab import histogram,semilogx,loglog
    f=i[(o>0) &(i>0)]
    l=i[(o<1) &(i>0)]

    d=i.sum()
    fNorm=f/d
    lNorm=l/d
    
    fHist=histogram(fNorm,logspace(-6,-2.0,30))
    lHist=histogram(lNorm,logspace(-6,-2.0,30))
    semilogx(fHist[1][1:],fHist[0].astype(float)/(fHist[0]+lHist[0]),**pltKwds)
