
# coding: utf-8

# In[34]:

from matplotlib import pyplot as plt
import numpy as np
import scipy.sparse as spar
from numpy import linalg as la
from scipy import linalg as spla
from scipy.sparse import linalg as sla
import math
from operator import itemgetter
from scipy import optimize as op
import random as rd
import networkx as nx
import os
import sys


# In[35]:

def randomGraph(n=10,fill=.3,diag=True):
    rank=True
    while rank:
        A=(np.random.rand(n,n)<fill)*1
        if diag:
            for j in xrange(n):
                A[j,j]=0
            if rankL(A)==n-1 and rankL(A.T)==n-1 and sum(A.sum(0)==0)==0 and sum(A.sum(1)==0)==0:
                rank=False
        else:
            rank=False
    return A


# In[36]:

def rankL(Q):
    D=np.diag(np.sum(Q,axis=1))
    L=D-Q
    return la.matrix_rank(L)

# In[38]:

def makeclusters(n,g,alp=.8,beta=.3):
    if (n%g)!=0:
        raise ValueError('number of nodes not divisible by number of groups')
    B=np.zeros((n,n))
    rank=True
    while rank:
        for i in xrange(g):
            B[i*(n/g):(i+1)*(n/g),i*(n/g):(i+1)*(n/g)]=randomGraph(n/g,alp)

        for i in xrange(g-1):
            for j in xrange(i+1,g):
                #j=i+1
                B[i*(n/g):(i+1)*(n/g),j*(n/g):(j+1)*(n/g)]=randomGraph(n/g,beta,False)
                B[j*(n/g):(j+1)*(n/g),i*(n/g):(i+1)*(n/g)]=randomGraph(n/g,beta,False)
        if rankL(B)==n-1 and rankL(B.T)==n-1 and sum(B.sum(0)==0)==0 and sum(B.sum(1)==0)==0:
            
            rank=False
    
    return B


# In[39]:

def pLaw(n,exp):
    rank=True
    while rank:
        z=nx.utils.create_degree_sequence(n,nx.utils.powerlaw_sequence,exponent=exp)
        #y=nx.utils.create_degree_sequence(10,nx.utils.powerlaw_sequence,exponent=1.5)
        G=nx.directed_configuration_model(z,z)
        Q=np.array(nx.adj_matrix(G.to_directed()).todense())
        Q[Q>0]=1
        for i in xrange(n):
            Q[i,i]=0
        if rankL(Q)==n-1 and rankL(Q.T)==n-1 and sum(Q.sum(0)==0)==0 and sum(Q.sum(1)==0)==0:
                rank=False
    return Q


# In[42]:

def adj_mat(datafile, n):
    """ Parse the data stored in 'datafile' and form the
    adjacencem matrix of the corresponding graph.
    'n' is the number of rows of the nxn sparse matrix formed. """
    W=np.zeros((n,n))
    with open(datafile,"r") as f:
        for L in f:
            words=L.strip().split()
            try:
                m=int(words[0])
                z=int(words[1])
                if m<n and z<n:
                    W[m,z]=1
            except:
                pass
    return W


# In[43]:

def f(x):
    return la.norm(x,ord=1)/la.norm(x,ord=2)
def g(x,fid,esp):
    return esp-la.norm(x-fid,ord=2)
def h(x,L,n,esp,lam):
    return esp-la.norm(np.dot(np.eye(n)*lam-L,x),ord=2)/la.norm(x,ord=2)
def h2(x,L,esp,y):
    return esp-la.norm(np.dot(L,x)-y,ord=2)/la.norm(x,ord=2)


# In[44]:

def areasphere(fun,args,n,num=1000):
    count=0
    for i in xrange(num):
        x=np.random.randn(n,n)
        y=np.sqrt((x**2).sum(1))
        if fun(y,*args)>0:
            count+=1
    return count*1./num
    


# In[45]:

def rightespvalue(fun,args,n,num=1000):
    count=0
    for i in xrange(num):
        x=np.random.randn(n,n)
        y=np.sqrt((x**2).sum(1))
        count+=fun(y,*args)
    return count*1./num
    


# In[46]:

def espchoser(A,n,lam):
    esp=0
    args=[A.T,n,esp,lam]
    esp=-1*rightespvalue(h,args,n,num=1000)
    args=[A.T,n,esp,lam]
    while areasphere(h,args,n,num=1000):
        esp*=.995
        args=[A.T,n,esp,lam]
    return esp


# In[49]:

def limiterPer(A,tol=10**(-6),which="normal"):
    n=len(A)
    Dinv=np.diag(1./A.sum(1))
    if which=="normal":
        Lt=A
    else:
        Lt=np.dot(Dinv,A)
    lams,vec=la.eig(Lt.T)
    index=np.argsort(lams)[-1]
    lam=np.real(lams[index])
    x0=np.real(vec[:,index])
    x=x0.copy()
    esp=espchoser(Lt.T,n,lam)
    const2={'type':'ineq','fun':h,'args':(Lt.T,n,esp,lam)}
    value=op.minimize(f,x0,constraints=(const2),options={'maxiter':2000})
    x=value.x/la.norm(value.x,2)
    return x,x0,esp


# In[50]:

def limiterPerChange(A,esps=np.arange(0,7,.1),tol=10**(-6),num=10,which="normal"):
    n=len(A)
    z=len(esps)
    change=np.zeros((num,z))
    ms=np.zeros(z)
    Dinv=np.diag(1./A.sum(1))
    if which=="normal":
        Lt=A
    else:
        Lt=np.dot(Dinv,A)
    lams,vec=la.eig(Lt.T)
    index=np.argsort(lams)[-1]
    lam=np.real(lams[index])
    x0=np.real(vec[:,index])
    x=x0.copy()
    for i,esp in enumerate(esps):
        const2={'type':'ineq','fun':h,'args':(Lt.T,n,esp,lam)}
        value=op.minimize(f,x0,constraints=(const2),options={'maxiter':2000})
        x=value.x/la.norm(value.x,2)
        change[:,i]=np.argsort(np.abs(x))[::-1][:num]
        ms[i]=np.sum(np.abs(x)>tol)
    return change,ms



# In[158]:

def strbuild2(A,eigvec,vec,titles):
    n=len(A)
    def getNeighbors(A,index):
        return list(np.flatnonzero(A[index]))
    sor=np.argsort(np.abs(eigvec))[::-1]
    eigvec*=np.sqrt(n)
    
    f=[]
    f.append('{\n')
    f.append('\t"titles":[')
    for t in titles:
        f.append('"'+t+'"')
        if t!=titles[-1]:
            f.append(',')
    f.append('],\n')
    f.append('\t"nodes":[\n')
    for i in xrange(n):
        f.append('\t\t{"name":"'+str(i)+'","eigvalue":'+("%.2f" % eigvec[i])+',"values":'+str(list(vec[i].astype(int)))+'}')
        if i!=n-1:
            f.append(',')
        f.append('\n')
    f.append('\t],\n')
    f.append('\t"links":[\n')
    for i in xrange(n):
        t=getNeighbors(A,i)
        for j in t:
            f.append('\t\t{"source":'+str(i)+',"target":'+str(j)+'}')
            if i!=n-1 or j!=t[-1]:
                f.append(',')
            f.append('\n')
    f.append('\t]\n')
    f.append('}')
    return ''.join(f)


# In[145]:

def centrailtyM(A,num=5):
    G=nx.DiGraph(A)
    ranks=np.zeros((num,8))
    ranks[:,0]=np.argsort(nx.in_degree_centrality(G).values())[::-1][:num]
    ranks[:,1]=np.argsort(nx.closeness_centrality(G).values())[::-1][:num]
    ranks[:,2]=np.argsort(nx.betweenness_centrality(G).values())[::-1][:num]
    ranks[:,3]=np.argsort(nx.eigenvector_centrality_numpy(G).values())[::-1][:num]
    ranks[:,4]=np.argsort(nx.katz_centrality_numpy(G,weight=None).values())[::-1][:num]
    ranks[:,5]=np.argsort(nx.pagerank_numpy(G,weight=None).values())[::-1][:num]
    return ranks


# In[146]:

def makevec(n,ranks):
    m=ranks.shape[1]
    vec=np.zeros((n,m))
    for i in xrange(m):
        for z in ranks[:,i]:
            vec[z,i]=1
    return vec


# In[159]:

def CodeBase(A,num=10,name="name"):
    if not os.path.exists(name):
        os.makedirs(name)
    n=len(A)
    titles=["Indegree","Closeness","Betweenness","Eigenvector","Katz","Pagerank","Sparse Adjacency","Sparse Random Walk"]
    ranks=centrailtyM(A,num=num)
    val,x0,esp1=limiterPer(A,tol=10**(-6),which="normal")
    sparrank=np.argsort(np.abs(val))[::-1][:num]
    ranks[:,-2]=sparrank
    val,x02,esp2=limiterPer(A,tol=10**(-6),which="rw")
    sparrank=np.argsort(np.abs(val))[::-1][:num]
    ranks[:,-1]=sparrank
    np.save(name+"/"+name,A)
    np.savetxt(name+"/"+name+"ranks.csv", ranks)
    vec=makevec(n,ranks)
    string=strbuild2(A,x0,vec,titles)
    with open(name+"/"+name+".json",'w') as f:
        f.write(string)
    change,ms=limiterPerChange(A,esps=np.linspace(0,esp1,11),tol=10**(-6),num=6,which="normal")
    np.savetxt(name+"/"+name+"norChange.csv", change.T)
    change,ms=limiterPerChange(A,esps=np.linspace(0,esp2,11),tol=10**(-6),num=6,which="rw")
    np.savetxt(name+"/"+name+"rwChange.csv", change.T)

def main(args):
    if args[2]=="-c":
        A=makeclusters(int(args[3]),int(args[4]),float(args[5]),float(args[6]))
    elif args[2]=="-p":
        A=A=pLaw(int(args[3]),float(args[4]))
    elif args[2]=="-r":
        A=randomGraph(int(args[3]),float(args[4]))
    else:
        A=adj_mat(args[2], int(args[3]))

    CodeBase(A,num=int(args[1]),name=args[0])

if __name__ == "__main__":
    args=sys.argv[1:]
    main(args)

