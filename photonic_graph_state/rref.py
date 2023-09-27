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
            else:
                oneitem = True
        elif len(distinct_pauli)==2:
            if 'I' in distinct_pauli:
                oneitem = True
            else:
                twoitem = True
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

def photonic_circuit_solver(state):
    stabs = state.stabilizers()
    n_e = num_emitters(state)
    n_p = state.size
    N = n_p+n_e
    for i in range(len(stabs)):
        stabs[i]+=n_e*'I'
    for i in range(n_e):
        stabs.append(n_p*'I'+i*'I'+'Z'+(n_e-i-1)*'I')
    target_state = Stabilizer(n_e+n_p,stabs)
    protocol = []
    for j in range(n_p,0,-1):
        height = heightfunction(target_state)
        photonindex = j-1
        d = height[j]-height[j-1]
        if d<0:
            'Time Reverse Measurement'
            gauss = target_state.gauss.tolist()
            indexfinder = [0 for i in range(N)]
            index=-1
            for i in range(n_p,N):
                indexfinder[i] = 1
                if indexfinder in gauss:
                    index = i
                    for k in range(N):
                        if indexfinder == gauss[k]:
                            a = k 
                indexfinder[i] = 0
            if index != -1:
                stab = [target_state.tab[a,index],target_state.tab[a,index+N]]
                if stab == [0,1]:
                    if target_state.signvector[a]==1:
                        target_state.clifford('X',index)
                        protocol.append(['X',index])
                elif stab == [1,0]:
                    if target_state.signvector[a]==1:
                        target_state.clifford('Z',index)
                        protocol.append(['Z',index])
                    target_state.clifford('H',index)
                    protocol.append(['H',index])
                elif stab == [1,1]:
                    if target_state.signvector[a]==0:
                        target_state.clifford('Z',index)
                        protocol.append(['Z',index])
                    target_state.clifford('S',index)
                    protocol.append(['S',index])
                    target_state.clifford('H',index)
                    protocol.append(['H',index])
                protocol.append(['Measure',index,photonindex])
                target_state.clifford('H',index)
                target_state.clifford('CNOT',index,photonindex)

            else:
                'More thorough rotation required'

        'Photon Absoprtion'

        'Identify Stabilizer'
        for i in range(N):
            toggler = True
            for k in range(photonindex):
                if target_state.tab[i,k]!=0 or target_state.tab[i,k+N]!=0:
                    toggler = False
                    break
            if target_state.tab[i,photonindex]==0 and target_state.tab[i,photonindex+N]==0:
                toggler = False
            if toggler:
                a = i
                break
        
        'Bring into Z'

        for i in range(photonindex,N):
            stab = [target_state.tab[a,i],target_state.tab[a,N+i]]
            if stab == [1,0]:
                protocol.append(['H',i])
                target_state.clifford('H',i)
            elif stab == [1,1]:
                protocol.append(['S',i])
                target_state.clifford('S',i)
                protocol.append(['H',i])
                target_state.clifford('H',i)
        
        'Disentangle all but one emitter'

        for i in range(n_p,N):
            stab = [target_state.tab[a,i],target_state.tab[a,N+i]]
            if stab == [0,1]:
                emitter = i
                break

        for i in range(emitter+1,N):
            if [target_state.tab[a,i],target_state.tab[a,N+i]]== [0,1]:
                protocol.append(['CNOT',i,emitter])
                target_state.clifford('CNOT',i,emitter)
        protocol.append(['CNOT',emitter,photonindex])
        target_state.clifford('CNOT',emitter,photonindex)

        'Clear out stabilizers'

        for i in range(N):
            if i!=a and [target_state.tab[i,photonindex],target_state.tab[i,N+photonindex]]!=[0,0]:
                target_state.row_add(a,i)

    rref(target_state)

    for i in range(n_p):
        if target_state.signvector[i]==1:
            target_state.clifford('X',i)
            protocol.append([['X'],i])
    return protocol.reverse()

def remove_sign(stabs):
    for i in range(len(stabs)):
        stabs[i] = stabs[i].lstrip('-')
    return stabs



