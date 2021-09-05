import numpy as np

# ----------------------------- #
# ----- Fourier transform ----- #
# ----------------------------- #

############################################################################################
##### Definition ###########################################################################
############################################################################################
####                                                                                    ####
#### Discrete Fourier transform                                                         ####
####   FT[x_{k}]  = a_{n}=\sum _{k=0} ^{N-1} x_{k} \exp{-i2 \pi nk/N}                   ####
#### Inverse discrete Fourier transform                                                 ####
####   IFT[a_{n}] = x_{k}=(1/N) \sum _{n=-(N-1)/2} ^{+(N-1)/2} a_{n} \exp{+i2 \pi nk/N} ####
####                                                                                    ####
############################################################################################

def ft_mod2_unfold(fs, b2s, mu):
    fs_neg=-fs[::-1]
    fs_unfold=fs_neg
    fs_unfold=np.append(fs_unfold, 0)
    fs_unfold=np.append(fs_unfold, fs)

    ln=2*len(b2s)+1
    b2_0=(mu*ln)**2
    b2s_flip=b2s[::-1]
    b2s_unfold=b2s_flip
    b2s_unfold=np.append(b2s_unfold, b2_0)
    b2s_unfold=np.append(b2s_unfold, b2s)

    return fs_unfold, b2s_unfold

def ft_mod2_fold(fs, bs):
    n=len(fs)
    # ----- Debug ----- #
    if len(fs)%2==0:
        print('Error')
        sys.exit()

    fs_fold=fs[int((n+1)/2):] #Only positive frequency
    bs_fold=bs[int((n+1)/2):]

    return fs_fold, bs_fold

# Convolution (direct)
def conv_calc_dir(bs, cs):
    # ----- Debug ----- #
    if not len(bs)==len(cs):
        print('Error')
        sys.exit()

    n=len(bs)
    ds=np.zeros(n)
    for i in range(n):
        bs_adj=bs[::-1]
        bs_adj=np.roll(bs_adj, -(int((n-1)/2)-i))
        d=np.sum(bs_adj*cs)
        ds[i]=d

    return ds

# Convolution (efficient)
def conv_calc_eff(bs, cs):
    # ----- Debug ----- #
    if not len(bs)==len(cs):
        print('Error')
        sys.exit()

    # ----- Adjust array for the use of numpy.fft.ifft ----- #
    n=len(bs)
    bs_adj=np.roll(bs, -int((n-1)/2))
    cs_adj=np.roll(cs, -int((n-1)/2))

    # ----- numpy.fft.ifft ----- #
    xs=np.fft.ifft(bs_adj)
    ys=np.fft.ifft(cs_adj)

    # ----- numpy.fft.fft ----- #
    ####################################################################
    #### Convolution theorem ###########################################
    ## a_{n} \odot b_{n} := \sum_{l} a_{n-l}b_{l} = N FT[x_{k}*y_{k}] ##
    ####################################################################
    ds=n*np.fft.fft(xs*ys)
    ds=ds.real
    ds=np.roll(ds, int((n-1)/2))

    return ds

