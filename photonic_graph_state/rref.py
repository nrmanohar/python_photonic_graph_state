from stabilizer import *
import math
import matplotlib.pyplot as plt

def rref(state):
    N=state.size
    K=N
    KU=0
    NL=0
    while NL<N-1 and KU<K-1:
        zeroitem = True
        oneitem = False
        twoitem = False
        r1=N
        r2=N
        for k in range(KU,K):
            if state.tab[k,NL]!=0 or state.tab[k,NL+N]!=0:
                r1 = k
                zeroitem = False
                oneitem = True
                break
        for k in range(r1,K):
            if state.tab[k,NL]!=0 or state.tab[k,NL+N]!=0:
                if state.tab[k,NL]!=state.tab[r1,NL] or state.tab[k,NL+N]!=state.tab[r1,NL+N]:
                    r2 = k
                    oneitem = False
                    twoitem = True
                    break
        if zeroitem:
            NL+=1
        elif oneitem:
            state.swap(KU,r1)
            for i in range(KU+1,K):
                if state.tab[i,NL]!=0 or state.tab[i,NL+N]!=0:
                    state.row_add(KU,i)
            KU+=1
            NL+=1
        elif twoitem:
            state.swap(KU,r1)
            state.swap(KU+1,r2)
            for i in range(KU+2,K):
                if state.tab[i,NL]!=0 or state.tab[i,NL+N]!=0:
                    if state.tab[i,NL]==state.tab[KU,NL] and state.tab[i,NL+N]==state.tab[KU,NL+N]:
                        state.row_add(KU,i)
                    elif state.tab[i,NL]==state.tab[KU+1,NL] and state.tab[i,NL+N]==state.tab[KU+1,NL+N]:
                        state.row_add(KU+1,i)
                    else:
                        state.row_add(KU,i)
                        state.row_add(KU+1,i)
            NL+=1
            KU+=2

def heightfunction(state):
    rref(state)
    N=state.size
    leftmost = []
    for i in range(state.size):
        for j in range(state.size):
            if state.tab[i,j]!=0 or state.tab[i,j+N]:
                leftmost.append(j+1)
                break
    height = []
    for i in range(state.size+1):
        count = sum(j > i for j in leftmost)
        height.append(state.size-i-count)
    return height

def plot_height(state):
    height = heightfunction(state)
    x_val = [i for i in range(state.size+1)]
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
    n_e = num_emitters(state)
    n_p = state.size
    N = n_e+n_p
    target_state = Stabilizer(N)
    for i in range(n_p):
        target_state.signvector[i] = state.signvector[i]
        for j in range(n_p):
            target_state.tab[i,j] = state.tab[i,j]
            target_state.tab[i,j+N] = state.tab[i,j+n_p]
    protocol = []
    for h in range(n_p,0,-1):
        height = heightfunction(target_state)
        photonindex = h-1
        d = height[h]-height[h-1]
        if d<0:
            rows = []
            for i in range(N):
                toggler = True
                for j in range(n_p):
                    if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                        toggler = False
                        break
                if toggler:
                    rows.append(i)
            sums=[]
            for row in rows:
                sum=0
                for j in range(n_p,N):
                    if target_state.tab[row,j]!=0 or target_state.tab[row,j+N]!=0:
                        sum+=1
                sums.append(sum)
            row = rows[sums.index(min(sums))]
            for i in range(n_p,N):
                if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                    emit = i
                    if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                        protocol.append(['S',i])
                        target_state.clifford('S',i)
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    break
            for i in range(emit+1,N):
                if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                    if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                        protocol.append(['S',i])
                        target_state.clifford('S',i)
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    protocol.append(['CNOT',i,emit])
                    target_state.clifford('CNOT',i,emit)
            if target_state.signvector[row]==1:
                target_state.clifford('X',emit)
                protocol.append(['X',emit])
            target_state.clifford('H',emit)
            target_state.clifford('CNOT',emit,photonindex)
            protocol.append(['Measure',emit,photonindex])
        rref(target_state)
        for i in range(N):
            toggler = True
            if target_state.tab[i,photonindex]==0 and target_state.tab[i,photonindex+N]==0:
                toggler = False
            if toggler:
                for j in range(photonindex):
                    if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                        toggler = False
                        break
            if toggler:
                row = i
                break
        emit = -1
        if target_state.tab[row,photonindex]==1 and target_state.tab[row,photonindex+N]==0:
            protocol.append(['H',photonindex])
            target_state.clifford('H',photonindex)
        elif target_state.tab[row,photonindex]==1 and target_state.tab[row,photonindex+N]==1:
            protocol.append(['S',photonindex])
            target_state.clifford('S',photonindex)
            protocol.append(['H',photonindex])
            target_state.clifford('H',photonindex)

        for i in range(n_p,N):
            if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                emit = i
                if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                    protocol.append(['H',i])
                    target_state.clifford('H',i)
                elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                    protocol.append(['S',i])
                    target_state.clifford('S',i)
                    protocol.append(['H',i])
                    target_state.clifford('H',i)
                break
        if emit!= -1:
            for i in range(emit+1,N):
                if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                    if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                        protocol.append(['S',i])
                        target_state.clifford('S',i)
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    protocol.append(['CNOT',i,emit])
                    target_state.clifford('CNOT',i,emit)
            target_state.clifford('CNOT',emit,photonindex)
            protocol.append(['Emission',emit,photonindex])
        for i in range(N):
            if target_state.tab[i,photonindex+N]!=0 or target_state.tab[i,photonindex]!=0:
                if i!=row:
                    target_state.row_add(row,i)
    rref(target_state)
    sums = []
    for i in range(n_p,N):
        sum=0
        for j in range(n_p,N):
            if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                sum+=1
        sums.append(sum)
    
    if max(sums)==1:
        decoupled = True
    else:
        decoupled = False
    
    while not decoupled:
        for i in range(2,N):
            if i in sums:
                minimum = i
                break
        row = n_p+sums.index(minimum)
        for i in range(n_p,N):
            if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                emit = i
                if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                    protocol.append(['H',emit])
                    target_state.clifford('H',emit)
                elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                    protocol.append(['S',emit])
                    target_state.clifford('S',emit)
                    protocol.append(['H',emit])
                    target_state.clifford('H',emit)
                break
        for i in range(emit+1,N):
            if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                    protocol.append(['H',i])
                    target_state.clifford('H',i)
                elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                    protocol.append(['S',i])
                    target_state.clifford('S',i)
                    protocol.append(['H',i])
                    target_state.clifford('H',i)
                target_state.clifford('CNOT',i,emit)
                protocol.append(['CNOT',i,emit])
        for i in range(n_p,N):
            if target_state.tab[i,emit]!=0 or target_state.tab[i,emit+N]!=0:
                if i!= row:
                    target_state.row_add(row,i)
        sums = []
        for i in range(n_p,N):
            sum=0
            for j in range(n_p,N):
                if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                    sum+=1
            sums.append(sum)
        if max(sums)==1:
            decoupled = True
        else:
            decoupled = False
    rref(target_state)
    for i in range(n_p,N):
        if target_state.tab[i,i]!=0:
            if target_state.tab[i,i]==1 and target_state.tab[i,i+N]==0:
                protocol.append(['H',i])
                target_state.clifford('H',i)
            elif target_state.tab[i,i]==1 and target_state.tab[i,i+N]==1:
                protocol.append(['S',i])
                target_state.clifford('S',i)
                protocol.append(['H',i])
                target_state.clifford('H',i)
    
    for i in range(N):
        if target_state.signvector[i]==1:
            target_state.clifford('X',i)
            protocol.append(['X',i])
    checker_state = Stabilizer(N)
    if np.array_equal(checker_state.tab,target_state.tab) and np.array_equal(checker_state.signvector,target_state.signvector):
        return protocol
    else:
        print('Something went wrong')
        return None

def emitter_cnot(state):
    emitter = num_emitters(state)
    procedure = photonic_circuit_solver(state)
    cnots = 0
    for i in range(len(procedure)):
        if procedure[i][0] == 'CNOT':
            cnots+=1
    data = [emitter,cnots]
    return data

def remove_sign(stabs):
    for i in range(len(stabs)):
        stabs[i] = stabs[i].lstrip('-')
    return stabs