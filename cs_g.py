#%% Finding and Correcting Cycle Slips
# Greg Starr, Sebastijan Mrak

#%% Import
import numpy as np
from scipy.interpolate import interp1d

#%% Old Functions

def quickCheck(data):
    """
    same as Sebastijan's algorithm, just checks for the 1 or 2 slips-in-a-row
    pattern within nonmathematical arbitrary bounds
    """
    l = np.empty(data.shape)
    l[:] = data
    #Get interaval of observation
    receiving_interval = np.where(np.isfinite(data))[0]
    start = receiving_interval[0]
    stop = receiving_interval[-1]
    l = l[start : stop]
    # Find NaNs in the time of obseravtions, caused by loss of locks and stick 
    idx = np.isfinite(l)
    l_corr = l[idx]
    # Calculate the Higher order differences on the dataset
    N_diff = 3
    data_diff = np.hstack(([np.nan]*N_diff,np.diff(l_corr,N_diff)))
    data_round = np.round(data_diff)
    
    cs_ix=[]
    cs_value = []
    for i in range(data_round.shape[0]-6):
        t = data_round[i:i+4]
        if(abs(t[1]-t[3])<=.2*abs(t[1]+t[3]) and 
           abs(t[2]+t[1]+t[3])<abs(.2*t[2]) and abs(t[1])>5 and abs(t[3])>5):
            cs_ix.append(i+1)
            cs_value.append(data_diff[i+1])
        else:
            test = np.array([0,-2*t[0]+t[3]-t[1],t[0]-2*t[3]-t[2],0])
            if(t@test<1 and abs(t[0])>5 and abs(t[3])>5):
                cs_ix.append(i)
                cs_value.append(data_diff[i])
                cs_ix.append(i+1)
                cs_value.append(data_diff[i+3])

    if len(cs_ix) > 0:
        cs_ix = np.hstack((cs_ix, 0))
    # Create new array for corrected values
    data_corrected = np.nan*np.zeros(len(l_corr))
    count = 0
    CS = len(cs_ix)  
    i = 0
    # cycle slip repair
    if CS > 0:
        # until first cycle slip, rewrite the values
        while i < cs_ix[0]:
            data_corrected[i] = l_corr[i]
            i += 1
        # After first cycle slip correct the values with sum of previous 
        # cycle slips
        while i < len(l_corr):
            if i == cs_ix[count]:
                count = count + 1
            data_corrected[i] = l_corr[i] - sum(np.round(cs_value[0 : count]))
            i += 1
    else:
        data_corrected = l_corr
    # Replace original input data with at non-NaN indexes 
    l[idx] = data_corrected

    return l

#%% Simulation
def slips(slipv):
    dim = slipv.shape[0]
    a = np.zeros((dim,dim))
    a[0]=np.arange(dim)**2
    a[0] += slipv
    for i in np.arange(1,dim):
        a[i] = np.hstack((np.nan,np.diff(a[i-1])))
    return a

def randomSlips(n):
    ix = np.random.choice(np.arange(20,80),n,replace=False)
    s = np.zeros(100)
    vs = []
    for i in ix:
        v = np.random.randint(-20,20)
        s[i:] += v
        vs.append(v)
        
    info = np.vstack((ix,vs))
    a = slips(s)
    return a,info
    
#%% New detect/correct ideally based on probability/detection theory

"""
fix slips in sequence or start from largest slips and work your way down?
"""

def newplan(sig, window):
    """
    passes through input signal "sig" checks if a value belongs to the random
    noise of the past "window" of samples by using chebyshev probability bound,
    could be better if noise is assumed Gaussian? frame as more of a detection 
    problem?
    """
    wave = np.array([1,-2,1])
    idx = np.isfinite(sig)
    out = np.empty(sig[idx].shape)
    out = sig[idx].ravel()
    val = []
    pos = []
    last_std = 10 #dummy value
    last_mean = 0 #dummy value
    for i in np.arange(window,out.shape[0]-3):
        if abs(out[i]-last_mean)>10*last_std: # Chebyshev prob <= .01
            val.append(out[i])
            pos.append(i)
            out[i:i+3] -= wave*int(out[i])
        last_std = np.std(out[i-window:i])
        last_mean = np.mean(out[i-window:i])
            
    return val,pos,out
    

def removeworst(sig,pos=None,val=None):
    """
    assumes the highest absolute value in a data set is a cycle slip, removes
    it according to the integer scalar A*[1,-2,1] structure, this is meant to 
    be used repetitively until the original signal satisfies some criteria
    
    alternatively if position and value are specified, this just removes the
    slip according to the same structure
    """
    out = sig.ravel()
    wave = np.array([1,-2,1])    
    if pos!=None and val!=None:
        out[pos-1:pos+2]-=wave*val
        return out,pos,val
    out = sig.ravel()
    wave = np.array([1,-2,1])
    idx = np.argmax(abs(sig))
    val = int(out[idx]/2)
    out[idx-1:idx+2] += wave*val
            
    return out,idx,val
    
    
def fix(data):
    """
    first full D/C, 
    
    1) finds NaNs in data, performs linear interpolation, then 
    treats those spots as 2 consecutive cycle slips of the same value
    
    2) fix overall data until the largest absolute value is less than
    (arbitrarily) 10
    """
    receiving_interval = np.where(np.isfinite(data))[0]
    start = receiving_interval[0]
    stop = receiving_interval[-1]
    
    data_copy = data.ravel()
    data_copy = data_copy[start:stop]
    not_nan = np.isfinite(data_copy)
    ixs = np.arange(data_copy.shape[0])
    interp = interp1d(ixs[not_nan],data_copy[not_nan])
    
    nan_pos = np.where(~not_nan)[0]
        
    filled_data = interp(ixs)
    
    third_diff = np.hstack(np.diff(filled_data,3))
    
    fixes = {}
    fixed = third_diff
    for pos in nan_pos:
        a = int(np.sqrt(np.mean(third_diff[pos-3:pos]**2)))*np.sign(third_diff[pos])
        fixes[pos] = -a
        fixes[pos+1] = -a
        fixed,idx,fix_value = removeworst(fixed,pos-2,a)
        fixed,idx,fix_value = removeworst(fixed,pos-1,a) 
    
    fixed,idx,fix_value = removeworst(fixed)
    fixes[idx+2] = fix_value
    
    while(abs(fix_value)>10):
        fixed,idx,fix_value = removeworst(fixed)
        idx+=2
        if(idx in fixes): 
            fixes[idx]+=fix_value
        else:
            fixes[idx] = fix_value
    for slip in fixes:
        filled_data[slip:] += fixes[slip]

    return filled_data,fixes
    
#%% New Plan: Between NaNs First
"""
first fix data between NaNs, then connect each interval and correct for slips
involving NaNs
"""
    
def betweenNan(sig,l):
    """
    First between NaN checker, uses change in variance from one window to the
    next to detect cycle slips, currently not grounded in math, if the variance
    goes up by 40 then it is considered a slip, I still need to work out what
    change in variance probabilistically constitutes a cycle slip
    l is window length
    """
    if(len(sig)<l): return [],[]
    out = np.empty(len(sig))
    out[:] = sig
    wave = np.array([1,-2,1])
    i=0
    v = np.empty(len(out)-l)
    idx = []
    val = []
    while i<len(out)-l-2:
        v[i] = np.var(out[i:i+l])
        if i>0 and v[i]-v[i-1]>40:
            a = int(out[i+l-1]/2)
            val.append(a)
            idx.append(i+l)
            out[i+l-2:i+l+1] += a*wave
            i-=l
        i+=1
    return idx,val


def fix2(data):
    """
    half D/C, uses the first between NaNs checker, does not fix NaNs yet,
    this one needs work because if the variance spikes from one window to the
    next, has it been caused by the "a" jump or the "-2a" jump? this needs 
    to check further one window to see which part of the structure the jump
    in variance is caused by. What of several slips in a row with no NaNs?
    
    I'm pretty sure this works best on smallish slips, or multiple small slips 
    in a row. Single large jumps fuck with it like in Mah8. Wierd though, because
    it works on Mah2
    """
    
    receiving_interval = np.where(np.isfinite(data))[0]
    start = receiving_interval[0]
    stop = receiving_interval[-1]
    
    data_copy = data.ravel()
    data_copy = data_copy[start:stop]
    not_nan = np.isfinite(data_copy) # a relic from interpolation
    
    nan_pos = np.where(~not_nan)[0]
    fixes = {}
        
    for i in range(len(nan_pos)+1):
        if i==0:
            sig = np.diff(data_copy[:nan_pos[0]],3)
            offset = 0
        elif i==len(nan_pos):
            sig = np.diff(data_copy[nan_pos[-1]:],3)
            offset = nan_pos[-1]+1
        else:
            sig = np.diff(data_copy[nan_pos[i-1]:nan_pos[i]],3)
            offset = nan_pos[i-1]+1
        idx,val = betweenNan(sig,10)
        for slip in range(len(idx)):
            if(slip in fixes): 
                fixes[idx[slip]+offset] = val[slip]
            else:
                fixes[idx[slip]+offset] = val[slip]
    
    for slip in fixes:
        data_copy[slip:] += fixes[slip]

    return data_copy,fixes
    
 
#%% Best Version

def betweenNan2(L1,L2,l):
    """
    so far has proved to be the best between NaN checker, it takes the
    normalized cross covariance between L1 and L2 derivatives which should
    go to zero or below during a cycle slip, I need to work out the math behind
    this though which should be kind of tough
    """
    idx,val = [],[]
    nan_pos = np.where(np.isnan(L2))[0]
    for i in range(len(nan_pos)+1):

        if i==0:
            td1 = np.diff(L1[:nan_pos[0]],3)
            td2 = np.diff(L2[:nan_pos[0]],3)
            offset = 0
        elif i==len(nan_pos):
            td1 = np.diff(L1[nan_pos[-1]+1:],3)
            td2 = np.diff(L2[nan_pos[-1]+1:],3)
            offset = nan_pos[-1]+1
        else:
            td1 = np.diff(L1[nan_pos[i-1]+1:nan_pos[i]],3)
            td2 = np.diff(L2[nan_pos[i-1]+1:nan_pos[i]],3)
            offset = nan_pos[i-1]+1
        if len(td1)-l<0: continue
        v = np.empty(len(td1)-l)
        for i in range(len(td1)-l):
            v[i] = np.correlate((td1[i:i+l]-np.mean(td1[i:i+l]))/np.std(td1[i:i+l]),
                                (td2[i:i+l]-np.mean(td2[i:i+l]))/np.std(td2[i:i+l]))
            if(i>0 and v[i]-v[i-1]>.9*l ):
                idx.append(i+offset)
                val.append(int(td2[i-2]/2))
    return idx,val
        
    
def fix3(L1,L2):
    """
    so far the best full D/C uses the cross covariance technique for between
    NaNs and interpolation for NaNs. I still need to implement the floor and 
    ceiling checker to determine which integer best corrects the data,
    furthermore the math could use a little work.
    """

    receiving_interval = np.where(np.isfinite(L2))[0]
    start = receiving_interval[0]
    stop = receiving_interval[-1]
    
    l1 = L1.ravel()
    l1 = l1[start:stop]
    l2 = L2.ravel()
    l2 = l2[start:stop]
    not_nan = np.isfinite(l2)
    ixs = np.arange(l2.shape[0])
    interp = interp1d(ixs[not_nan],l2[not_nan])
    filled_data = interp(ixs)
    
    #FIND BETWEEN NANS
    fixes = {}
    idx,val = betweenNan2(l1,l2,20)
    for slip in range(len(idx)):
        if(idx[slip] in fixes): 
            fixes[idx[slip]] = val[slip]
        else:
            fixes[idx[slip]] = val[slip]

    #FIND NANS
    third_diff = np.diff(filled_data,3)
    fixed = third_diff
    nan_pos = np.where(~not_nan)[0]
    for pos in nan_pos:
        a = int(np.sqrt(np.mean(third_diff[pos-3:pos]**2)))*np.sign(third_diff[pos])
        fixes[pos] = -a
        fixes[pos+1] = -a
        fixed,idx,fix_value = removeworst(fixed,pos-2,a)
        fixed,idx,fix_value = removeworst(fixed,pos-1,a) 
    
    #FIX SLIPS
    for slip in fixes:
        filled_data[slip:] += fixes[slip]

    return filled_data,fixes
    