from __future__ import division
from numpy import *
from scipy import cos, sin, sqrt


def interp_opt_functions_with_orbit(b1, NInterPoints=10, Lmin=1.):
    if hasattr(b1, 'BETX'):
        # Capital
        BETX = b1.BETX
        L = b1.L
        S = b1.S
        DX = b1.DX
        DPX = b1.DPX
        DY = b1.DY
        DPY = b1.DPY
        BETX = b1.BETX
        ALFX = b1.ALFX
        BETY = b1.BETY
        ALFY = b1.ALFY
        ANGLE = b1.ANGLE
        KEYWORD = b1.KEYWORD
        NAME = b1.NAME
        K1L = b1.K1L
        X = b1.X
        Y = b1.Y
        PX = b1.PX
        PY = b1.PY
    elif hasattr(b1, 'betx'):
        # Non capital
        L = b1.l
        S = b1.s
        DX = b1.dx
        DPX = b1.dpx
        DY = b1.dy
        DPY = b1.dpy
        BETX = b1.betx
        ALFX = b1.alfx
        BETY = b1.bety
        ALFY = b1.alfy
        ANGLE = b1.angle
        KEYWORD = b1.keyword
        NAME = b1.name
        K1L = b1.k1l
        X = b1.x
        Y = b1.y
        PX = b1.px
        PY = b1.py
    else:
        raise ValueError('What??!!!')

    j = -1
    s = []
    keyword = []
    name = []
    dx = []
    dpx = []
    dy = []
    dpy = []
    betx = []
    alfax = []
    bety = []
    alfay = []
    x = []
    y = []
    px = []
    py = []

    if L[0] != 0:
        raise(ValueError('The length of the first element must be zero!!!!'))
    
    for i in range(len(BETX)):
        NInterp = NInterPoints+1
        j = j+1
        if (L[i] > Lmin):
            s += (list(S[i-1]+linspace(0, NInterp, NInterp+1) /
                       float(NInterp)*(S[i]-S[i-1])))
            # ~ raise(ValueError)
            keyword += ((NInterp+1)*[KEYWORD[i-1]])
            name += ((NInterp+1)*[NAME[i-1]])
            dx += [DX[i-1]]
            dpx += [DPX[i-1]]
            dy += [DY[i-1]]
            dpy += [DPY[i-1]]
            betx += [BETX[i-1]]
            alfax += [ALFX[i-1]]
            bety += [BETY[i-1]]
            alfay += [ALFY[i-1]]
            x += [X[i-1]]
            y += [Y[i-1]]
            px += [PX[i-1]]
            py += [PY[i-1]]
    
            j = j+1
    
            # transfer matrizes
            if KEYWORD[i] == 'SBEND':
                # see transfer matrix for sbend dipole (Lee, equ.2.37 + equ.2.54)
                M11x = cos(ANGLE[i]/NInterp)
                M12x = L[i]/ANGLE[i]*sin(ANGLE[i]/NInterp)
                M13x = L[i]/ANGLE[i]*(1-cos(ANGLE[i]/NInterp))
                M21x = -ANGLE[i]/L[i]*sin(ANGLE[i]/NInterp)
                M22x = cos(ANGLE[i]/NInterp)
                M23x = sin(ANGLE[i]/NInterp)
    
                M11y = 1
                M12y = L[i]/NInterp
                M13y = 0
                M21y = 0
                M22y = 1
                M23y = 0
    
            elif(KEYWORD[i] == 'RBEND' and ANGLE[i] != 0):
                # see transfer matrix for sbend dipole (Lee, equ.2.37 + equ.2.54 + exercise 2.2.2)
                # edge focusing included ...
    
                # edge focusing:
                delta = ANGLE[i]/2
                M11x = 1
                M12x = 0
                M13x = 0
                M21x = tan(delta)/(L[i]/ANGLE[i])
                M22x = 1
                M23x = 0
    
                M11y = 1
                M12y = 0
                M13y = 0
                M21y = -tan(delta)/(L[i]/ANGLE[i])
                M22y = 1
                M23y = 0
    
                # modify alfa at entrance of Rbend in order to account for edge
                # focusing
                alfax[j-1] = -M11x*M21x*betx[j-1] + \
                    (M11x*M22x+M12x*M21x)*alfax[j-1] - \
                    M12x*M22x*(1+alfax[j-1]**2)/betx[j-1]
                alfay[j-1] = -M11y*M21y*bety[j-1] + \
                    (M11y*M22y+M12y*M21y)*alfay[j-1] - \
                    M12y*M22y*(1+alfay[j-1]**2)/bety[j-1]
    
                # normal transport through bending magnet ...
                M11x = cos(ANGLE[i]/NInterp)
                M12x = L[i]/ANGLE[i]*sin(ANGLE[i]/NInterp)
                M13x = L[i]/ANGLE[i]*(1-cos(ANGLE[i]/NInterp))
                M21x = -ANGLE[i]/L[i]*sin(ANGLE[i]/NInterp)
                M22x = cos(ANGLE[i]/NInterp)
                M23x = sin(ANGLE[i]/NInterp)
    
                M11y = 1
                M12y = L[i]/NInterp
                M13y = 0
                M21y = 0
                M22y = 1
                M23y = 0
    
            elif (K1L[i] == 0):
                # see transfer matrix for drift (Lee, equ.2.37 + equ.2.54)
                M11x = 1
                M12x = L[i]/NInterp
                M13x = 0
                M21x = 0
                M22x = 1
                M23x = 0
    
                M11y = 1
                M12y = L[i]/NInterp
                M13y = 0
                M21y = 0
                M22y = 1
                M23y = 0
    
            # ~ elif(K1L[i]!=0):
            else:
                # quadrupole case ....
                # see transfer matrix for quadrupole (Lee, equ.2.37 + equ.2.54)
                M11x = cos(sqrt(K1L[i]/L[i])*L[i]/NInterp)
                M12x = 1/sqrt(K1L[i]/L[i])*sin(sqrt(K1L[i]/L[i])*L[i]/NInterp)
                M13x = 0
                M21x = -sqrt(K1L[i]/L[i])*sin(sqrt(K1L[i]/L[i])*L[i]/NInterp)
                M22x = cos(sqrt(K1L[i]/L[i])*L[i]/NInterp)
                M23x = 0
    
                M11y = cos(sqrt(-K1L[i]/L[i])*L[i]/NInterp)
                M12y = 1/sqrt(-K1L[i]/L[i])*sin(sqrt(-K1L[i]/L[i])*L[i]/NInterp)
                M13y = 0
                M21y = -sqrt(-K1L[i]/L[i])*sin(sqrt(-K1L[i]/L[i])*L[i]/NInterp)
                M22y = cos(sqrt(-K1L[i]/L[i])*L[i]/NInterp)
                M23y = 0
    
            # propagate optical functions
            for k in range(NInterPoints):
                betx += [M11x**2*betx[j-1]-2*M11x*M12x *
                         alfax[j-1]+M12x**2*(1+alfax[j-1]**2)/betx[j-1]]
                alfax += [-M11x*M21x*betx[j-1]+(M11x*M22x+M12x*M21x)
                          * alfax[j-1]-M12x*M22x*(1+alfax[j-1]**2)/betx[j-1]]
    
                bety += [M11y**2*bety[j-1]-2*M11y*M12y *
                         alfay[j-1]+M12y**2*(1+alfay[j-1]**2)/bety[j-1]]
                alfay += [-M11y*M21y*bety[j-1]+(M11y*M22y+M12y*M21y)
                          * alfay[j-1]-M12y*M22y*(1+alfay[j-1]**2)/bety[j-1]]
    
                dx += [dx[j-1]*M11x + dpx[j-1]*M12x + M13x]
                dpx += [dx[j-1]*M21x + dpx[j-1]*M22x + M23x]
    
                dy += [dy[j-1]*M11y + dpy[j-1]*M12y + M13y]
                dpy += [dy[j-1]*M21y + dpy[j-1]*M22y + M23y]
    
                x += [x[j-1]*M11x + px[j-1]*M12x]
                px += [x[j-1]*M21x + px[j-1]*M22x]
    
                y += [y[j-1]*M11y + py[j-1]*M12y]
                py += [y[j-1]*M21y + py[j-1]*M22y]
    
                j = j+1
    
            alfax += [ALFX[i]]
            alfay += [ALFY[i]]
            betx += [BETX[i]]
            bety += [BETY[i]]
            dx += [DX[i]]
            dpx += [DPX[i]]
            dy += [DY[i]]
            dpy += [DPY[i]]
            x += [X[i]]
            y += [Y[i]]
            px += [PX[i]]
            py += [PY[i]]
    
            # ~ print 'ciao Nicolas'
    
        else:
    
            s += [S[i]]
            keyword += [KEYWORD[i]]
            name += [NAME[i]]
            betx += [BETX[i]]
            bety += [BETY[i]]
            dx += [DX[i]]
            dpx += [DPX[i]]
            dy += [DY[i]]
            dpy += [DPY[i]]
            alfay += [ALFY[i]]
            alfax += [ALFX[i]]
            x += [X[i]]
            y += [Y[i]]
            px += [PX[i]]
            py += [PY[i]]
    
    s = real(array(s))
    betx = real(array(betx))
    bety = real(array(bety))
    dx = real(array(dx))
    dpx = real(array(dpx))
    dy = real(array(dy))
    dpy = real(array(dpy))
    alfay = real(array(alfay))
    alfax = real(array(alfax))
    x = real(array(x))
    y = real(array(y))
    px = real(array(px))
    py = real(array(py))
    
    return s, keyword, name, betx, bety, dx, dy, dpy, alfay, alfax, x, y, px, py
