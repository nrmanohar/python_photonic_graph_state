from stabilizer import *
import math
import matplotlib.pyplot as plt

def rref(state):
    N=state.size
    K=N
    KU=0
    NL=0
    while NL<N-1 and KU< K-1:
        stabs = state.stabilizers()
        stabs = remove_sign(stabs)
        column_stabs = []
        for k in range(KU,K):
            column_stabs.append(stabs[k][NL])
        identity = False
        oneitem = False
        twoitem = False
        distinct_pauli = set(column_stabs)
        if len(distinct_pauli)==1:
            if 'I' in distinct_pauli:
                identity = True
        elif len(distinct_pauli)<=2:
            oneitem = True
        else:
            twoitem = True
        if identity:
            NL+=1
        elif oneitem:
            for k in range(KU,K):
                if stabs[k][NL]!='I':
                    state.swap(KU,k)
                    break
            stabs = state.stabilizers()
            stabs = remove_sign(stabs)
            for k in range(KU+1,K):
                if stabs[k][NL]!='I':
                    state.row_add(KU,k)
            NL+=1
            KU+=1
        elif twoitem:
            for k in range(KU,K):
                if stabs[k][NL]!='I':
                    reference1 = stabs[k][NL]
                    state.swap(KU,k)
                    break
            for k in range(KU,K):
                if stabs[k][NL]!='I' and stabs[k][NL]!= reference1:
                    reference2 = stabs[k][NL]
                    state.swap(KU+1,k)
                    break
            stabs = state.stabilizers()
            stabs = remove_sign(stabs)
            for k in range(KU+2,K):
                if stabs[k][NL]=='I':
                    continue
                elif stabs[k][NL]==reference1:
                    state.row_add(KU,k)
                elif stabs[k][NL]==reference2:
                    state.row_add(KU+1,k)
                elif stabs[k][NL]!='I':
                    state.row_add(KU,k)
                    state.row_add(KU+1,k)
            NL+=1
            KU+=2

def heightfunction(state):
    rref(state)
    state.gaussian()
    gauss = state.gauss
    leftmost = []
    for i in range(state.size):
        for j in range(state.size):
            if gauss[i,j]!=0:
                leftmost.append(j+1)
                break
    height = []
    for i in range(state.size+1):
        count = sum(j > i for j in leftmost)
        height.append(state.size-i-count)
    return height

def plot_height(state):
    height = heightfunction(state)
    x_val = []
    for i in range(state.size+1):
        x_val.append(i)
    tickers = range(math.floor(min(height)), math.ceil(max(height))+1)
    plt.grid(color = 'blue', linewidth = 0.5)
    plt.plot(x_val,height,color='blue')
    plt.scatter(x_val,height,color='blue')
    plt.yticks(tickers)
    plt.title('Target Height Function')
    plt.xlabel('x')
    plt.ylabel('h(x)')
    plt.show()

def num_emitters(state):
    height = heightfunction(state)
    emitters = max(height)
    return emitters

def remove_sign(stabs):
    for i in range(len(stabs)):
        stabs[i] = stabs[i].lstrip('-')
    return stabs



