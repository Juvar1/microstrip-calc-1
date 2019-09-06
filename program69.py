#programs from late 60'

#https://apps.dtic.mil/dtic/tr/fulltext/u2/708720.pdf
#Green's function:
#https://apps.dtic.mil/dtic/tr/fulltext/u2/a117364.pdf
#IBM System 360 Scientific Routine Package
#http://www.ebyte.it/library/codesnippets/IBM_System360_SSP.html

import math

DIST=[0,[0]*101,[0]*101]
XINT=0  #global variables
D1=0
D2=0
D3=0
D4=0
N1=0
N2=0
N3=0
N4=0 #global variables

E1=16.0
WD=0.5
SD=0.2
I=1
D=0.3048 #d board thickness
U1= 4.7#Ur
F= 3000000#f frequency


#input data for MSTRIP2 program
WH1  = 0.1 # starting point
DELW = 0.2 # step size
NT   = 20  # lines in result table
R    = 0.0 # 0.0 = no upper ground plane
DIEK = 9.6 # Ur dielectric constant
SH1  = 0.4 # spacing/height ratio
AIR  = 1   # 0 = single strip, 1 = coupled strips

#function subprogram
def gint(U,R,CO,BO,DIEK):
    V=(R-1.0)*U
    W1=CO*U
    W2=math.cosh(U)
    W3=math.sinh(U)
    W4=math.cos(BO*U)
    if (R != 0):
        W4=W4*math.sinh(V)
        DEN=W3*math.cosh(V)+DIEK*W2*math.sinh(V)
    if (R == 0):
        DEN=W3+DIEK*W2
    GINT=math.sin(W1)/W1*W3/U*W4/DEN
    return GINT

#compute integral
def dqg32(XL,XU,R,CO,BO,DIEK):
    A=0.5*(XU+XL)
    B=XU-XL
    C=0.49863193092474078e0*B
    Y=0.35093050047350483e-2*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.49280575577263417e0*B
    Y=Y+0.8137197365452835e-2*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.48238112779375322e0*B
    Y=Y+0.12696032654631030e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.46745303796886984e0*B
    Y=Y+0.17136931456510717e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.44816057788302606e0*B
    Y=Y+0.21417949011113340e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.42468380686628499e0*B
    Y=Y+0.25499029631188088e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.39724189798397120e0*B
    Y=Y+0.29342046739267774e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.36609105937014484e0*B
    Y=Y+0.32911111388180923e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.33152213346510760e0*B
    Y=Y+0.36172897054424253e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.29385787862038116e0*B
    Y=Y+0.39096947893535153e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.25344995446611470e0*B
    Y=Y+0.41655962113473378e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.21067563806531767e0*B
    Y=Y+0.43826046502201906e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.16593430114106382e0*B
    Y=Y+0.45586939347881942e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.11964368112606854e0*B
    Y=Y+0.46922199540402283e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.7223598079139825e-1*B
    Y=Y+0.47819360039637430e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK))
    C=0.24153832843869158e-1*B
    Y=B*(Y+0.48270044257363900e-1*(gint(A+C,R,CO,BO,DIEK)+gint(A-C,R,CO,BO,DIEK)))
    return Y

# computes the sine and cosine integral
def sici(X):
    Z=abs(X)
    if (Z-4.0 <= 0):
        Y=(4.0-Z)*(4.0+Z)
        SI=X*(((((1.753141e-9*Y+1.568988e-7)*Y+1.374168e-5)*Y+6.939889e-4)*\
            Y+1.964882e-2)*Y+4.395509e-1)
        CI=((5.772156e-1+math.log(Z))/Z-Z*(((((1.386985e-10*Y+1.584996e-8)*\
            Y+1.725752e-6)*Y+1.185999e-4)*Y+4.990920e-3)*Y+1.315308e-1))*Z
    else:
        SI=math.sin(Z)
        Y=math.cos(Z)
        Z=4.0/Z
        U=((((((((4.048069e-3*Z-2.279143e-2)*Z+5.515070e-2)*Z-7.261642e-2)*\
            Z+4.987716e-2)*Z-3.332519e-3)*Z-2.314617e-2)*Z-1.134958e-5)*\
            Z+6.250011e-2)*Z+2.583989e-10
        V=(((((((((-5.108699e-3*Z+2.819179e-2)*Z-6.537283e-2)*Z+\
            7.902034e-2)*Z-4.400416e-2)*Z-7.945556e-3)*Z+2.601293e-2)*Z-\
            3.764000e-4)*Z-3.122418e-2)*Z-6.646441e-7)*Z+2.500000e-1
        CI=Z*(SI*V-Y*U)
        SI=-Z*(SI*U+Y*V)+1.570796
    ret=dict()
    ret["si"]=SI
    ret["ci"]=CI
    return ret

# To solve a system of simultaneous linear equations with symmetric
# coefficient matrix, upper triangular part of which is assumed to be
# stored columnwise.
def dgels(A,M,N,EPS):#,AUX):
    R=[1.0]*(M+1) # return values
    AUX=[0]*M     # auxiliary storage array
    if(M <= 0):
        IER=-1
        return (R,IER)
    #search for greatest main diagonal element
    IER=0
    PIV=0.0
    L=0
    for K in range(1,M+1):
        L=L+K
        TB=abs(A[L])
        if (TB-PIV > 0):
            PIV=TB
            I=L
            J=K
    TOL=EPS*PIV
    #Main diagonal element A[I]=A[0.0] is first pivot element.
    #PIV kontains the absolute value of A[I].
    #Start elimination loop
    LST=0
    NM=N*M
    LEND=M-1
    for K in range(1,M+1):
        #test on usefulness of symmetric algorithm
        if (PIV <= 0):
            IER=-1
            return (R,IER)
        if (IER == 0):
            if (PIV-TOL <= 0):
                IER=K-1
        LT=J-K
        LST=LST+K
        #Pivot row reduction and row interchange in right hand side R
        PIVI=1.0/A[I]
        for  L in range(K,NM+1,M):
            LL=L+LT
            TB=PIVI*R[LL]
            R[LL]=R[L]
            R[L]=TB
        #is elimination terminated
        if (K-M < 0):
            #row and column interchange and pivot row reduction in matrix A.
            #Elements of pivot column are saved in auxiliary vector AUX.
            LR=int(LST+(LT*(K+J-1))/2)
            LL=LR
            L=LST
            for II in range(K,LEND+1):
                L=L+II
                LL=LL+1
                if (L-LR == 0):
                    A[LL]=A[LST]
                    TB=A[L]
                elif (L-LR > 0):
                    LL=L+LT
                    TB=A[LL]
                    A[LL]=A[L]
                elif (L-LR < 0):
                    TB=A[LL]
                    A[LL]=A[L]
                AUX[II]=TB
                A[L]=PIVI*TB
            #save column interchange information
            A[LST]=LT
            #element reduction and search for next pivot
            PIV=0.0
            LLST=LST
            LT=0
            for II in range(K,LEND+1):
                PIVI=-AUX[II]
                LL=LLST
                LT=LT+1
                for LLD in range(II,LEND+1):
                    LL=LL+LLD
                    L=LL+LT
                    A[L]=A[L]+PIVI*A[LL]
                LLST=LLST+II
                LR=LLST+LT
                TB=abs(A[LR])
                if (TB-PIV > 0):
                    PIV=TB
                    I=LR
                    J=II+1
                for LR in range(K,NM+1,M):
                    LL=LR+LT
                    R[LL]=R[LL]+PIVI*R[LR]
    #back substitution and back interchange
    if (LEND < 0):
        IER=-1
        return (R,IER)
    elif (LEND == 0):
        return (R,IER)
    else:
        II=M
        for I in range(2,M+1):
            LST=LST-II
            II=II-1
            L=A[LST]+0.5
            for J in range(II,NM+1,M):
                TB=R[J]
                LL=J
                K=LST
                for LT in range(II,LEND+1):
                    LL=LL+1
                    K=K+LT
                    TB=TB-A[K]*R[LL]
                K=int(J+L)
                R[J]=R[K]
                R[K]=TB
    return (R,IER)

# generates X
def xgen(M,WH,SH1):
    X=[0]*60
    WHM=WH/M
    for I in range(M):
        X[I+1]=I*WHM
    if (SH1 != 0.0):
        for II in range(1,2*M):
            X[II+M]=II*WHM+SH1
    return X

# creates PHI with free air
def mphi(WH,M,S,X):
    PHI=[0]*60
    AM=M
    INDEX1=3*M-1
    if (S == 0):
        INDEX1=M
    EXWON=WH/(2.0*AM)
    WYWON=1.0
    for K in range(1,INDEX1+1):
        XO=X[K]
        EXP=XO+EXWON
        EXN=XO-EXWON
        WYP=2.0*WYWON
        PHI[K]=(EXN/(2.0*EXWON))*math.log((EXN**2)/(EXN**2+WYP**2))-\
            (EXP/(2.0*EXWON))*math.log((EXP**2)/(EXP**2+WYP**2))+\
            (WYP/EXWON)*(math.atan(EXP/WYP)-math.atan(EXN/WYP))
    return PHI

# M = substrips count
# WH = 
# SH1 = 0.4 spacing/height ratio
# DIEK = 9.6 Ur
# R = 0
# creates PHI with dielectric substrate other than free air
def mgreen(M,WH,SH1,DIEK,R,X):
    PHI=[0]*60
    CO=WH/M*0.5
    X1=5.0
    INT=2
    H=X1/float(INT)
    if (SH1 == 0.0):
        MA=M
    else:
        MA=3*M-1
    for MM in range(1,MA+1):
        BO=X[MM]
        YTOT=0.0
        XU=0.0
        XL=0.0
        # compute first integral
        for I in range(INT):
            XU+=H
            YTOT+=dqg32(XL,XU,R,CO,BO,DIEK)
            XL+=H
        AI1=YTOT
        #compute second integral
        S1=(CO+BO)*X1
        S2=(BO-CO)*X1
        AI2A=math.sin(S1)/X1-(CO+BO)*sici(S1)["ci"]
        AI2B=math.sin(S2)/X1-(BO-CO)*sici(S2)["ci"]
        AI2=AI2A-AI2B
        PHI[MM]=4.0*(AI1+1.0/((1.0+DIEK)*CO*2.0)*AI2)
    return PHI

#creates A and B
def amat(S,M,PHI):
    A=[0]*211
    B=[0]*211
    for I in range(1,M+1):
        INDEX1=M+1-I
        for J in range(1,INDEX1+1):
            NF1=0
            for K in range(1,J+1):
                NF1=NF1+K
            NF2=0
            if (I > 2):
                INDEX2=I-2
                for L in range(1,INDEX2+1):
                    NF2=NF2+L
            N=NF1+NF2+(I-1)*J
            INDEX3=3*M+2*(1-J)-I
            A[N]=PHI[I]+S*PHI[INDEX3]
            B[N]=PHI[I]-S*PHI[INDEX3]
    return [A,B]

#calculate final output
def output(CAP1E,CAPKE,CAP1O,CAPKO,WH,AIR):
    C=2.99792458
    EFFKE=CAPKE/CAP1E
    #print("%f %f"%(CAPKE,CAP1E))
    RKE=math.sqrt(EFFKE)
    ZOE=1.0e+4/(C*CAP1E*RKE)
    VELE=(1.0/RKE)*C
    if (AIR == 1):
        EFFKO=CAPKO/CAP1O
        #print("%f %f"%(CAPKO,CAP1O))
        RKO=math.sqrt(EFFKO)
        ZOO=1.0e+4/(C*CAP1O*RKO)
        VELO=(1.0/RKO)*C
        print("W/H1=%.3f ZOE=%.3f ZOO=%.3f VE=%.3f VO=%.3f EFFE=%.3f EFFO=%.3f"%(WH,ZOE,ZOO,VELE,VELO,EFFKE,EFFKO))
    else:
        print("W/H1=%.3f ZOE=%.3f VE=%.3f EFFE=%.3f"%(WH,ZOE,VELE,EFFKE))
    return

def mstrip2(WH1,DELW,NT,R,DIEK,SH1,AIR):
    M=20       # substrips count
    N=1        # number of right hand side vectors for DGELS function
    EPS=1.0e-7 # loss of significance limit level for DGELS function
    for K in range(NT):    # NT = lines on the output file
        WH=WH1+K*DELW      # starting point + K * step size
        X = xgen(M,WH,SH1)
        CAP1E=0
        CAPKE=0
        CAP10=0
        CAPKO=0
        for IJ in range(2):
            if (IJ == 0):
                ADIEK=1.0  # permittivity of air
            if (IJ == 1):
                ADIEK=DIEK # substrate permittivity
            if (R == 0.0 and ADIEK == 1.0): #R = 0 no upper ground plane
                PHI = mphi(WH,M,AIR,X)
            else:
                PHI = mgreen(M,WH,SH1,ADIEK,R,X)
            for JJ in range(AIR+1):
                AB = amat(AIR,M,PHI)
                #
                (V,IER) = dgels(AB[JJ],M,N,EPS)
                if (IER != 0):
                    print("IER= %.0f IN SUBROUTINE DGELS, SO THE CHARGE DENSITY COULD NOT BE CALCULATED TO THE PRECISION OF %.13f DIGITS" % (IER,EPS))
                CAPSUM=0.0
                #print(V)
                for I in range(M):
                    CAPSUM+=V[I+1]
                #print(CAPSUM)
                CC=CAPSUM*111.256
                if (JJ == 0 and IJ == 0):
                    CAP1E=CC
                if (JJ == 0 and IJ == 1):
                    CAPKE=CC
                if (JJ == 1 and IJ == 0):
                    CAP1O=CC
                if (JJ == 1 and IJ == 1):
                    CAPKO=CC
        output(CAP1E,CAPKE,CAP1O,CAPKO,WH,AIR)
    return

#integration subroutine
def dqtfe(H,Y,Z,NDIM):
    SUM2=0.0
    if (NDIM-1 < 0):
        return Z
    elif (NDIM-1 > 0):
        HH=0.500*H
        for I in range(2,NDIM+1):
            SUM1=SUM2
            SUM2=SUM2+HH*(Y[I]+Y[I-1])
            Z[I-1]=SUM1
    Z[NDIM]=SUM2
    return Z

#fourier transform program for coupled strip
def tsum(XN,ARG1):
    global DIST,F,U1,E1,WD,SD,D,XINT,I
    MS=(SD/2.0)/0.025+1
    MW=(SD/2.0+WD)/0.025
    S=0.0
    M=int(MS)-1
    while True:
        M=M+1
        SIGMA=DIST[I][M]
        ARG3=ARG1*M*0.025
        if (I == 1):
            C=math.sin(ARG3)
        if (I == 2):
            C=math.cos(ARG3)
        S=S+SIGMA*C
        if (M >= MW):
            break
    TSUM=S
    return TSUM

#integrand evaluation program for coupled strip
def fn(XN,SQ):
    global F,U1,E1,WD,SD,D,XINT,I
    if (XN == 0):
        FN=0.0
        return FN
    DL=D/30.0
    T=U1*E1-SQ
    Q=(SQ-1.0)/T
    R=(U1*E1-1.0)/T
    ARG1=0.2*math.pi*XN/WD
    ARG2=0.1*math.pi*XN*(SD/WD+1.0)
    DTDEF=math.tanh(ARG1)
    A=tsum(XN,ARG1)
    if (I == 1):
        E=math.sin(ARG2)
    if (I == 2):
        E=math.cos(ARG2)
    if (F == 0.0):
        PART=Q*U1*DTDEF-1.0
        XNUM=1.0/XN*A*PART*E
        XDEN=R**2-E1/SQ*PART*(Q*1.0/DTDEF-1.0/E1)
    else:
        A1=(ARG1/(F*DL))**2
        A2=(2.0*math.pi)**2*T
        SIGN=1.0
        B12=(A2-A1)
        if (B12 < 0):
            SIGN=-1.0
        B1=math.sqrt(SIGN*B12)
        Z=math.tan(F*DL*B1)
        if (B12 < 0):
            Z=math.tanh(F*DL*B1)
        B2=math.sqrt((2.0*math.pi)**2*(SQ-1.0)+A1)
        PART=SIGN*(Q*U1*Z+SIGN*B2/B1)
        XNUM=A*PART*B1*E
        XDEN=A1*R**2+SIGN*E1/SQ*B12*PART*(Q*1.0/Z-1.0/E1*B2/B1)
    if (XDEN == 0.0):
        print("NUM=%.5f DEN=%.5f" % (XNUM,XDEN))
    FN=XNUM/XDEN
    return FN

#integrand evaluation program
def fn1(XN,SQ):
    global F,U1,E1,WD,SD,D,XINT,I
    if (XN == 0):
        FN=0.0
        return FN
    DL=D/30.0
    T=U1*E1-SQ
    Q=(SQ-1.0)/T
    R=(U1*E1-1.0)/T
    ARG1=0.2*math.pi*XN/WD
    ARG3=0.1*math.pi*XN #*(SD/WD+1.0)
    DTDEF=math.tanh(ARG1)
    P3=193.5092066/XN**3
    P1=9.54929668/XN
    P2=60.79271019/XN**2
    A=1.0/XN*(P3*(P1-P3)*math.cos(ARG3)*(2.0-P2)*math.sin(ARG3))
    if (F == 0.0):
        PART=Q*U1*DTDEF-1.0
        XNUM=1.0/XN*A*PART
        XDEN=R**2-E1/SQ*PART*(Q*1.0/DTDEF-1.0/E1)
    else:
        A1=(ARG1/(F*DL))**2
        A2=(2.0*math.pi)**2*T
        SIGN=1.0
        B12=(A2-A1)
        if (B12 < 0):
            SIGN=-1.0
        B1=math.sqrt(SIGN*B12)
        Z=math.tan(F*DL*B1)
        if (B12 < 0):
            Z=math.tanh(F*DL*B1)
        B2=math.sqrt((2.0*math.pi)**2*(SQ-1.0)+A1)
        PART=SIGN*(Q*U1*Z+SIGN*B2/B1)
        XNUM=A*PART*B1
        XDEN=A1*R**2+SIGN*E1/SQ*B12*PART*(Q*1.0/Z-1.0/E1*B2/B1)
    if (XDEN == 0.0):
        print("NUM=%.5f DEN=%.5f" % (XNUM,XDEN))
    FN=XNUM/XDEN
    return FN

#summation program
def summ(SQ):
    global F,U1,E1,WD,SD,D,XINT,I
    global D1,D2,D3,D4,N1,N2,N3,N4
    FF=Z=[0]*801
    XB=0.0
    for J in range(1,N1+1):
        XN=XB+(J-1)*D1
        FF[J]=fn(XN,SQ)
    Z=dqtfe(D1,FF,Z,N1)
    S1=Z[N1]
    XB=25.0
    for J in range(1,N2+1):
        XN=XB+(J-1)*D2
        FF[J]=fn(XN,SQ)
    Z=dqtfe(D2,FF,Z,N2)
    S2=Z[N2]
    XB=50.0
    for J in range(1,N3+1):
        XN=XB+(J-1)*D3
        FF[J]=fn(XN,SQ)
    Z=dqtfe(D3,FF,Z,N3)
    S3=Z[N3]
    XB=75.0
    for J in range(1,N4+1):
        XN=XB+(J-1)*D4
        FF[J]=fn(XN,SQ)
    Z=dqtfe(D4,FF,Z,N4)
    S4=Z[N4]
    XINT=S1+S2+S3+S4
    SUM=XINT
    return SUM

#this subroutine solves transcendental equation of one unknown
#by half interval search method
def trans(A,B,EPS):
    FA=summ(A)
    IER=0
    FB=summ(B)
    if (FA*FB > 0):
        IER=1
        print("IER=%f" % IER)
        return (1,IER)
    if (abs(FA) <= EPS):
        ROOT=A
        return (ROOT,IER)
    if (abs(FB) <= EPS):
        ROOT=B
        return (ROOT,IER)
    for IK in range(1,41):
        X=(A+B)/2.0
        FX=summ(X)
        if (abs(FX) <= EPS):
            ROOT=X
            return (ROOT,IER)
        if (FA*FX < 0):
            B=X
            FB=FX
            continue
        if (FA*FX > 0):
            A=X
            FA=FX
            continue
    
#general program
def generalProgram():
    global F,U1,E1,WD,SD,D,XINT,I
    global D1,D2,D3,D4,N1,N2,N3,N4
    N1=601
    N2=301
    N3=201
    N4=101
    for JP in range(1,3):
        TRY=[0,1.0,((E1+1.0)/2.0),E1]
        D1=25.0/600.0
        D2=25.0/300.0
        D3=25.0/200.0
        D4=25.0/100.0
        for LL in range(1,5):
            if (LL+1 > 3):
                print("ROOT TROUBLE")
                return
            AB=TRY[LL]+.01
            AE=TRY[LL+1]-.01
            EPS=.00001
            (ROOT,IER)=trans(AB,AE,EPS)
            SQ=ROOT
            if (IER == 0):
                break
        print("PSI=%.5f VALUE OF INTEGRAL=%.5f" % (SQ, XINT))



#run programs
mstrip2(WH1,DELW,NT,R,DIEK,SH1,AIR)
generalProgram()
