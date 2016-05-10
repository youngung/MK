from for_lib import vm
import numpy as np
import matplotlib.pyplot as plt
pi=np.pi
sin=np.sin
cos=np.cos

def main():
    fig=plt.figure()
    ax=fig.add_subplot(111)
    th = np.linspace(-pi,pi)
    x=cos(th)
    y=sin(th)
    z=np.zeros(len(th))
    s=np.array([x,y,z,z,z,z]).T
    print s.shape

    X=[]
    Y=[]
    for i in xrange(len(s)):
        phi, dphi, d2phi = vm(s[i])
        print s[i]
        ys=s[i]/phi
        X.append(ys[0])
        Y.append(ys[1])
    ax.plot(X,Y)

if __name__=='__main__':
    main()
