"""
This script generates the opme (operator matrix element) file for operators: l_tau, spin, GT_1, GT_0, Fermi_0, Fermi_1
Fermi: roughly equals to isospin operator with assigned delta isospin (deltaT)
GT: Gamow-Teller operator with assigned delta isospin (deltaT), see https://en.wikipedia.org/wiki/Beta_decay_transition for details
l: orbital angular momentum, tau: isospin
See opem/ folder for the generated files, and also examples/na23 for example
"""


from sympy.physics.wigner import wigner_6j
import numpy as np


# input 2 lists (len=3), output float
def opme(b, a, op, tau=0):
    if len(a) != 3 or len(b) != 3:
        raise ValueError('Input lists must have the length 3')
    na, la, ja = a
    nb, lb, jb = b
    if na != nb or la != lb:
        return 0

    sixj = wigner_6j(1/2,jb,lb,ja,1/2,1).n(6)
    spinME = (-1)**(3/2+la+jb)*np.sqrt(2*ja+1)*np.sqrt(2*jb+1)*np.sqrt(12)*sixj

    if op == 'l_tau':
        return la * np.sqrt(2*tau+1)/np.sqrt(2)

    if op == 'spin':
        return spinME

    if op == 'gt1' or op == 'gt0':
        return spinME * np.sqrt(2*tau+1)/np.sqrt(2)

    if op == 'fermi0' or op == 'fermi1':
        if ja != jb:
            return 0
        return np.sqrt(2*tau+1)/np.sqrt(2) * 2


def generate_str(orbits, op):
    orbits_num = len(orbits)
    out = 'iso\n'
    out += ' ' + str(orbits_num) + '\n'
    for n in range(orbits_num):
        out += ' '*2 + str(n+1) + ' '*4
        out += str(orbits[n][0]) + ' '*4
        out += str(orbits[n][1]) + ' '*4
        out += str(orbits[n][2]) + ' '*4
        out += '\n'

    deltaJ, deltaT = 1, 0
    if op == 'gt1' or op == 'fermi1':
        deltaT = 1
    out += str(deltaJ) + ' '*5 + str(deltaT) + '\n'
    for i in range(orbits_num):
        for j in range(orbits_num):
            me = opme(orbits[i], orbits[j], op=op, tau=deltaT)
            if me != 0:
                out += ' '*4 + str(i+1) + ' '*4 + str(j+1) + ' '*2 + str(me)
                out += '\n'

    return out


# output opme in pn format for GT transition (see p.62 in manual)
def generate_pn(orbits, T0):
    orbits_num = len(orbits)
    out = 'pns\n'
    out += ' ' + str(orbits_num) + '\n'
    for n in range(orbits_num):
        out += ' '*2 + str(n+1) + ' '*4
        out += str(orbits[n][0]) + ' '*4
        out += str(orbits[n][1]) + ' '*4
        out += str(orbits[n][2]) + ' '*4
        out += '\n'

    out += '1     2\n'
    for i in range(orbits_num):
        for j in range(orbits_num):
            me_T0 = opme(orbits[i], orbits[j], op='gt0')
            me_T1 = np.sqrt(2) * me_T0
            p, n =me_pn(me_T0, me_T1, T=T0) # T0 is nucleus isospin
            if me_T0 != 0:
                out += ' '*4 + str(i+1) + ' '*4 + str(j+1) + ' '*2
                out += str(p) + ' '*3 + str(n)
                out += '\n'
    return out


# convert isospin (t0,t1) to pn format for GT operator, T is the nucleus isospin
def me_pn(t0,t1,T):
    a0 = np.sqrt(2*T+1) / np.sqrt(2)
    a1 = np.sqrt(2*T+1) * np.sqrt(T*(T+1)) / T / np.sqrt(6)
    p = t0/(2*a0) + t1/(2*a1)
    n = t0/(2*a0) - t1/(2*a1)
    return p, n


# [n,l,j]
p0 = [[0,1,0.5],[0,1,1.5]]
sd = [[0,2,1.5],[0,2,2.5],[1,0,0.5]]
pf = [[0,3,3.5],[1,1,1.5],[0,3,2.5],[1,1,0.5]]
sdpf = [[0,2,2.5],[1,0,0.5],[0,2,1.5],[0,3,3.5],[1,1,1.5],[0,3,2.5],[1,1,0.5]]
jj55 = [[1,4,3.5],[2,2,2.5],[2,2,1.5],[3,0,0.5],[1,5,5.5]]


# op="spin", "fermi0", "fermi1", "gt0", "gt1", "l_tau"
# orbits=p0, sd, pf, sdpf, jj55
if __name__ == "__main__":
    print(generate_str(orbits=jj55, op='l_tau'))