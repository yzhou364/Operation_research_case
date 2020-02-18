# -*- coding: utf-8 -*-
from gurobipy import *
import time
import numpy as np

#Setting the parameter for the case study
J = 1
T_BAR = 12
N_BAR = 1
I_t = 5.32
PAI = [50.17, 35.17, 35.15, 57.17, 90.00, 146.94, 95.00, 95.00, 90.61, 60.39,
95.62, 60.25]
TAO = 2
C_1 = 75
D_1 = 75
Q_1_LOW = 8.5
Q_1_HIGH = 42
Q_MINUS = 70
Q_PLUS = 70
V_LOW = 1500
V_HIGH = 3300
V_0 = V_t_HIGH = 2107.85
S_HIGH = GRB.INFINITY
W_1 = 0.0583333
Y_1 = 0
E_1 = 2
Q_10 = 0
G_10 = 0
U_10 = 0
Q_1_MINUS = -27
P_1_MINUS = -21.5
THETA_LOW = 0

def func(q1t, v1t):
    R0 = 0.01
    L = [4.09863600116008, -1.25535942295343, 0.160530264942775, -9.76201903589132e-3, 0.000309429429972963,
         -492928898248035e-6, 3.11519548768e-8]
    L_low = 385
    K = [307.395, 3.88e-1, -4.37e-4, 2.65e-7, -8.87e-11, 1.55e-14, -1.11e-18]
    summation = 0
    temp1 = 0
    for hh in range(7):
        temp1 = L[hh] * q1t ** hh
        temp = 0
        for kk in range(7):
            temp = temp + K[kk] * (v1t ** kk) *(10000**kk)

        temp = temp - L_low - R0 * q1t * q1t
        temp1 = temp1 * temp
        summation += temp1
    p1t = 9.81 / 1000 * q1t * summation
    return p1t
#GRB model
m1 = Model("MILP1")

v={}
s={}
q={}
p={}
w={}
w_bar={}
g={}
y={}
y_bar={}
u={}

#Adding variables
for t in range(T_BAR):
    v[t] = m1.addVar(vtype = 'C', name = 'water volume in t')
    s[t] = m1.addVar(vtype = 'C', name = 'spillage in t' )
    for j in range(J):
        q[j,t] = m1.addVar(vtype = 'C',name = 'water flow')
        p[j,t] = m1.addVar(vtype = 'C',name = 'power generated')
        w[j,t] = m1.addVar(vtype = 'B',name = 'shot down phrase')
        w_bar[j,t] = m1.addVar(vtype = 'B',name = 'start-up phrase')
        g[j,t] = m1.addVar(vtype = 'B',name = 'status')
        y[j,t] = m1.addVar(vtype = 'B',name = 'shutdown phrase')
        y_bar[j,t] = m1.addVar(vtype = 'B',name = 'start-up phrase')
        u[j,t] = m1.addVar(vtype = 'B',name = 'status of pump')
        
#Adding connstraint
for j in range(J):
    for t in range(T_BAR):
        q[j,t].lb = Q_1_MINUS
        q[j,t].ub = Q_1_HIGH   
        p[j,t].lb = P_1_MINUS
        p[j,t].ub = 50#
for t in range(T_BAR):
    v[t].lb = V_LOW
    v[t].ub = V_HIGH
    s[t].lb = 0
    s[t].ub = S_HIGH
#Fix the initial variable 
q[0,0].lb=q[0,0].ub = 0
v[0].lb=v[0].ub = V_0
g[0,0].lb=g[0,0].ub = 0
u[0,0].lb=u[0,0].ub = 0


m1.setObjective(quicksum(TAO*PAI[t]*p[0,t]-C_1*w_bar[0,t]-(D_1+PAI[t]*E_1)*y_bar[0,t] for t in range(T_BAR)))
m1.modelSense = GRB.MAXIMIZE

#Constraint 8.2-8.15
m1.addConstr(v[11] - V_0 == 0,name = "8.2")     
m1.addConstrs((v[t] - v[t-1] -0.36*TAO*(I_t-q[0,t]-s[t]) == 0 for t in range(1,T_BAR)),name= "8.3")
m1.addConstrs((q[0,t]-(Q_1_MINUS*u[0,t]+Q_1_LOW*g[0,t]) >= 0 for t in range(T_BAR)),name="8.4")
m1.addConstrs((q[0,t]-(Q_1_MINUS*u[0,t]+Q_1_HIGH*g[0,t]) <= 0 for t in range(T_BAR)),name="8.5")
m1.addConstrs(((q[0,t]+q[0,t-1])+TAO*Q_MINUS >=0 for t in range(1,T_BAR)),name="8.6")
m1.addConstrs(((q[0,t]+q[0,t-1])-TAO*Q_PLUS <=0 for t in range(1,T_BAR)),name="8.7")
m1.addConstrs((s[t]-(W_1*w_bar[0,t]+Y_1*y_bar[0,t]) >= 0 for t in range(T_BAR)),name="8.8")
m1.addConstrs((q[0,t]+s[t]-THETA_LOW >=0 for t in range(T_BAR)), name ="8.9")
m1.addConstrs((g[0,t]-g[0,t-1]-(w_bar[0,t]-w[0,t]) ==0 for t in range(1,T_BAR)), name="8.10")
m1.addConstrs((w_bar[0,t]+w[0,t] <= 1 for t in range(T_BAR)), name ="8.11")
m1.addConstrs((u[0,t]-u[0,t-1]-(y_bar[0,t]-y[0,t]) == 0 for t in range(1,T_BAR)),name="8.12")
m1.addConstrs((y_bar[0,t]+y[0,t] <= 1 for t in range(T_BAR)),name="8.13")
m1.addConstrs((g[0,t]+u[0,t] <= 1 for t in range(T_BAR)),name="8.14")
m1.addConstrs((u[0,t] <= N_BAR-1 for t in range(T_BAR)), name="8.15")


#Add more variables
Z = 5
R = 5
d ={}
z ={}
lamda = {}
for t in range(T_BAR):
    for r in range(R):
        d[t,r] = m1.addVar(vtype = 'B',name = 'membership')
for j in range(J):
    for t in range(T_BAR):
        for i in range(Z):
            z[j,t,i] = m1.addVar(vtype = 'B',name = 'contiguity')
            lamda[j,t,i] =m1.addVar(vtype = 'C',name = 'weight')
            lamda[j,t,i].lb = 0
            lamda[j,t,i].ub = 1
#Add more parameter
H_INT = []
num = 0
for i in range(R+1):
    H_INT.append(V_LOW + num*(V_HIGH-V_LOW)/R)
    num+=1

Q1_JR=[]
P1 = np.zeros((Z, R))
P1_delta = []
for i in range(Z):
    Q1_JR.append((Q_1_LOW + i * (Q_1_HIGH-Q_1_LOW) / (Z)))
for i in range(Z):
    for r in range(R):
        P1[i][r] = func(Q1_JR[i], (H_INT[r] + H_INT[r + 1]) / 2)
        

for r in range(R):
    P1_delta.append(max((P1[i][R- 1] - P1[i][r]) for i in range(Z)))

            
####
#Constraint 8.16-8.25
m1.addConstrs((q[0,t]-quicksum(Q1_JR[i]*lamda[0,t,i] for i in range(Z))-Q_1_MINUS*u[0,t] == 0 for t in range(T_BAR)),name ="8.18")
m1.addConstrs((quicksum(lamda[0,t,i] for i in range(Z))-g[0,t] == 0 for t in range(T_BAR)),name="8.19") ###make infeasible
m1.addConstrs((lamda[0,t,i] - z[0,t,i] <= 0 for t in range(T_BAR) for i in range(Z)),name="8.20")
m1.addConstrs((z[0,t,i]+z[0,j,k] <= 1 for t in range(T_BAR) for k in range(Z) for i in range(k-1)),name ="8.21") ###
m1.addConstrs((quicksum(d[t,r] for r in range(R)) == 1 for t in range(T_BAR)),name = "8.22")
m1.addConstrs((p[0,t] - quicksum(P1[i,r]*lamda[0,t,i] for i in range(Z)) - P_1_MINUS*u[0,t] - P1_delta[r]*(1-d[t,r]) <= 0 for t in range(T_BAR) for r in range(R)),name = "8.23")
m1.addConstrs((v[t] - quicksum(H_INT[r-1]*d[t,r] for r in range(1,R)) >= 0 for t in range(T_BAR)),name ="8.24")
m1.addConstrs((v[t] - quicksum(H_INT[r]*d[t,r] for r in range(R)) <= 0 for t in range(T_BAR)),name ="8.25")



m1.update()
m1.optimize()
print(v[11])