# -*- coding: utf-8 -*-
"""
Created on Fri Mar 02 22:00:26 2018

@author: d0v0b
"""
import numpy as np

sqrt_eps = 1.4901e-08
eps = 2.2204e-16
ex = np.array([[1], [0], [0]])
ey = np.array([[0], [1], [0]])
ez = np.array([[0], [0], [1]])
p01 = np.array([[0],[0],[0]])
p12 = np.array([[0.32],[0],[0.78]])
p23 = np.array([[0],[0],[1.075]])
p34 = np.array([[1.142],[0],[0.2]])
p45 = np.array([[0],[0],[0]])
p56 = np.array([[0],[0],[0]])
p6T = np.array([[0.2],[0],[0]])


def rotx(q):
    R = np.matrix([[1,0,0],[0,np.cos(q),-np.sin(q)],[0,np.sin(q),np.cos(q)]])
    return R

def roty(q):
    R = np.matrix([[np.cos(q),0,np.sin(q)],[0,1,0],[-np.sin(q),0,np.cos(q)]])
    return R

def rotz(q):
    R = np.matrix([[np.cos(q),-np.sin(q),0],[np.sin(q),np.cos(q),0],[0,0,1]])
    return R

def hat(k):
    khat = np.matrix([[0,-1*k[2],k[1]],[k[2],0,-1*k[0]],[-1*k[1],k[0],0]])
    return khat

def subproblem0(p,q,k):
    if ((k.transpose().dot(p)>sqrt_eps) or (k.transpose().dot(q)>sqrt_eps)):
        print 'k must be perpendicular to p and q'
    
    ep = p/np.linalg.norm(p)
    eq = q/np.linalg.norm(q)
    
    theta = 2*np.arctan2(np.linalg.norm(ep-eq),np.linalg.norm(ep+eq))
    
    if k.transpose().dot(np.cross(p,q,axis=0))<0:
      theta=-theta
    
    return theta

def subproblem1(k,p,q):
    
    if np.linalg.norm(p-q)< sqrt_eps:
        theta=0
        return theta
              
    k=k/np.linalg.norm(k)
    pp = p-(p.transpose().dot(k))*k
    qp = q-(q.transpose().dot(k))*k

    epp = pp/np.linalg.norm(pp)
    eqp = qp/np.linalg.norm(qp)
    
    theta = subproblem0(epp,eqp,k)
#    theta = np.arctan2(k.transpose()*(np.cross(epp,eqp)),epp.transpose()*eqp)
    
    tol=1e-2
    if (np.abs(np.linalg.norm(p)-np.linalg.norm(q))>tol):
      print '*** Warning *** ||p|| and ||q|| must be the same!!!'
    
    return theta

def subproblem2(k1,k2,p,q):
    
    k12 = np.dot(k1.transpose(),k2)[0][0]
    pk = np.dot(p.transpose(),k2)[0][0]
    qk = np.dot(q.transpose(),k1)[0][0]
    
    
    # check if solution exists
    denominator = np.square(k12)-1
    if np.abs(denominator) < eps:
        theta1=[]
        theta2=[]
        print 'subproblem2 no solution (***1***)'
        return [theta1, theta2] 
    
    a = np.matmul(np.array([[k12,-1],[-1,k12]]),np.array([pk,qk]).T)/denominator  
    bb = (np.square(np.linalg.norm(p)) - np.square(np.linalg.norm(a)) - 2*a[0]*a[1]*k12)
    
    if np.abs(bb.all())<eps:
        bb=0
    
    if bb.all()<0:
        theta1=[]
        theta2=[]
        print 'subproblem2 no solution (***2***)'
        return [theta1, theta2]
    
    # check if there is only 1 solution
    gamma = np.sqrt(bb)/np.linalg.norm(np.cross(k1,k2,axis=0))
#    if np.abs(gamma)<eps:
#      c1 = np.array([k1, k2, cross(k1,k2)]).dot(a.append(gamma))
#      theta2 = [subproblem1(k2,p,c1)]
#      theta1 = [-1*subproblem1(k1,q,c1)]  

    # general case: 2 solutions    
    theta1 = np.zeros([2,1])
    theta2 = np.zeros([2,1])
    c1 = np.matmul(np.concatenate((k1, k2, np.cross(k1,k2,axis=0)),axis=1),np.array([np.append(a,gamma)]).T)
    c2 = np.matmul(np.concatenate((k1, k2, np.cross(k1,k2,axis=0)),axis=1),np.array([np.append(a,-1*gamma)]).T)
    theta2[0] = subproblem1(k2,p,c1)
    theta2[1] = subproblem1(k2,p,c2)
    
    theta1[0] = -1*subproblem1(k1,q,c1)
    theta1[1] = -1*subproblem1(k1,q,c2)
    
    return [theta1, theta2]


def subproblem3(k,p,q,d):

    pp = p-(k.transpose().dot(p))*k
    qp = q-(k.transpose().dot(q))*k

    dpsq = d**2 - (k.transpose().dot(p-q))**2
    
    if dpsq<0:
        theta=[]
        return theta
    
    if dpsq==0:
        theta = subproblem1(k,pp/np.linalg.norm(pp),qp/np.linalg.norm(qp))
        return theta
      
    bb = (np.linalg.norm(pp)**2+np.linalg.norm(qp)**2-dpsq)/(2*np.linalg.norm(pp)*np.linalg.norm(qp))
    
    if np.abs(bb)>1.0:
        theta=[]
        return theta
    
    phi = np.arccos(bb)
    
    theta0 = subproblem1(k,pp/np.linalg.norm(pp),qp/np.linalg.norm(qp))
    theta = np.zeros([2,1])
    
    theta[0] = theta0+phi
    theta[1] = theta0-phi
    return theta

def subproblem4(p,q,k,d):
    
    a = p.transpose()*hat(k)*q
    b = -p.transpose()*hat(k)*hat(k)*q
    c = d - (p.transpose().dot(q) - b)
    
    phi = np.arctan2(b,a)
    
    if np.abs(c/np.sqrt(a**2 + b**2))>1.0:
        theta = []
        return theta
    
    theta = np.zeros([2,1])
    psi = np.arcsin(c/np.sqrt(a**2+b**2))
    
    theta[0] = -phi+psi
    theta[1] = -phi-psi+np.pi
    return theta

def abb_irb_6640_IK(R06, p0T):

    Q_all = np.array([])
    ## Q1
    ## d=p'*rot(k,theta)*q
    d1 = ey.T.dot(p12+p23+p34)
    v1 = p0T-R06*p6T
    k1 = -ez
    p1 = ey
    Q1 = subproblem4(p1,v1,k1,d1)
    
    for q1 in Q1:
        
        ## Q3
        ## norm(q-exp(k x theta) p) = d
        d3 = np.linalg.norm(rotz(-q1)*(p0T-R06*p6T)-p12)
        v3 = p23
        k3 = ey 
        p3 = -p34
        Q3 = subproblem3(k3,p3,v3,d3)
        
        for q3 in Q3:
            
            ## Q2
            ## exp(k x theta) p = q
            v2 = np.array(rotz(-q1)*(p0T-R06*p6T)-p12)
            k2 = ey
            p2 = np.array(p23 + roty(q3)*p34)
            Q2 = subproblem1(k2,p2,v2)
    
            ## Q5 Q6
            ## exp(k1 x theta1) * exp(k2 x theta2) p = q
            v4 = np.array(roty(-q3)*roty(-Q2)*rotz(-q1)*R06*ex)
            k41 = ex
            k42 = ey
            p4 = ex
            Q4,Q5 = subproblem2(k41,k42,p4,v4)
            
            for q4,q5 in zip(Q4,Q5):
                
                ## Q4
                ## exp(k x theta) p = q
                k6 = ex
                p6 = ey 
                v6 = np.array(roty(-q5)*rotx(-q4)*roty(-q3)*roty(-Q2)*rotz(-q1)*R06*ey)
                Q6 = subproblem1(k6,p6,v6)
    
    #            Q_all.append([q1,Q2,q3,q5,q6,Q4])
                Q_all = np.hstack((Q_all, np.array([q1,Q2,q3,q4,q5,Q6]).reshape(6)))
    
    return np.array(Q_all.reshape(-1,6))




#q = [np.pi/3,np.pi/4,0,np.pi/5,np.pi/2,np.pi]
#R01 = rotz(q[0])
#R12 = roty(q[1])
#R23 = roty(q[2])
#R34 = rotx(q[3])
#R45 = roty(q[4])
#R56 = rotx(q[5])
#R06 = R01*R12*R23*R34*R45*R56
#p0T = np.array(R01*p12+ R01*R12*p23 + R01*R12*R23*p34 + R06*p6T)
#
#Q_all = abb_irb_6640_IK(R06, p0T)
#print np.array(Q_all[0])
    


#S1 = rotx(0.5)
#S2 = roty(0.5)
#S3 = rotz(0.5)
#H = hat(p12)
    
#print subproblem1(ez,p12/np.linalg.norm(p12),p34/np.linalg.norm(p34))
#print subproblem3(ey,np.array([[-1.142],[0],[-0.2]]),p23,1.7117)
#print subproblem2(ex,ey,p12/np.linalg.norm(p12),ez)
#print np.matmul(S1,np.matmul(S2,S3))

#print subproblem4(ey,np.array([[1.0145],[1.7572],[0.8740]]),-ez,0)
#print S1*S2*S3