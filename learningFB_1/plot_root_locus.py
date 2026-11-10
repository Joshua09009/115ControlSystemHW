"""
產生 RootLocus_Ex510_Feedback.tex 用到的兩張圖：
  figs/three_cases.png : omega = 2, 1, 0.5 三種情況的根軌跡
  figs/w2_vs_K.png     : 兩個 w^2 隨 K 的變化（omega = 2, 0.5）

系統：L(s) = 1 / ( s (s+2) [ (s+1)^2 + omega^2 ] )
方法：對一串 K 值用 numpy.roots 解特徵方程式 a(s) + K = 0，把所有根畫在 s 平面上。
（figs/sketch.jpg 是課堂投影片截圖，不是由本程式產生。）

執行：python plot_root_locus.py   （需要 numpy、matplotlib）
"""
import os
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

FIG=os.path.join(os.path.dirname(os.path.abspath(__file__)),'figs')
os.makedirs(FIG,exist_ok=True)

def locus(w, ax, title, marks):
    a=np.poly1d([1,2,0]); b=np.poly1d([1,2,1+w*w]); A=a*b
    Ks=np.concatenate([np.linspace(0,30,20000),np.linspace(30,2000,4000)])
    R=np.array([np.sort_complex(np.roots(A.coeffs+np.r_[0,0,0,0,K])) for K in Ks])
    ax.plot(R.real.ravel(),R.imag.ravel(),'.',ms=0.8,color='#1f77b4')
    P=A.r; ax.plot(P.real,P.imag,'x',ms=9,mew=2,color='k')
    for K,lab in marks:
        r=np.roots(A.coeffs+np.r_[0,0,0,0,K]); ax.plot(r.real,r.imag,'o',ms=5,mfc='none',color='r')
    t=np.linspace(-4,4,2); 
    for s in [1,-1]: ax.plot(-1+t, s*t,'--',lw=0.7,color='gray')
    ax.axhline(0,color='k',lw=0.5); ax.axvline(0,color='k',lw=0.5)
    ax.set_xlim(-4,2); ax.set_ylim(-3,3); ax.set_aspect('equal'); ax.grid(alpha=.3)
    ax.set_title(title,fontsize=11); ax.set_xlabel('Real axis')
fig,axs=plt.subplots(1,3,figsize=(12,4.4))
locus(2,axs[0],r'$\omega=2$ (Example 5.10)',[(4,''),(6.25,'')])
locus(1,axs[1],r'$\omega=1$ (boundary)',[(1,'')])
locus(0.5,axs[2],r'$\omega=0.5$',[(0.25,''),(1.25**2/4,'')])
axs[0].set_ylabel('Imaginary axis')
plt.tight_layout(); plt.savefig(os.path.join(FIG,'three_cases.png'),dpi=200)
# w^2 trajectory plot for omega=0.5 and 2
fig,axs=plt.subplots(1,2,figsize=(10,3.6))
for ax,w in zip(axs,[2,0.5]):
    K=np.linspace(0,(1+w*w)**2/4,400); d=np.sqrt((1+w*w)**2-4*K)
    ax.plot(K,(1-w*w+d)/2,label=r'$w_+^2$'); ax.plot(K,(1-w*w-d)/2,label=r'$w_-^2$')
    ax.axhline(0,color='k',lw=0.6); ax.axvline(w*w,color='r',ls=':',lw=1)
    ax.set_title(rf'$\omega={w}$: real roots of $w^2$ vs $K$'); ax.set_xlabel('K'); ax.legend(); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'w2_vs_K.png'),dpi=200)
