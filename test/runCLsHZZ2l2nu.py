import os, re, sys, commands, glob

import ROOT
ROOT.gROOT.LoadMacro("fitRvsCLs_C.so")

jobname = "hzz2l2nu"
datacards = "SM4_hzz2l2nu*.txt"
seed=1235
sEstimateTimePerJob = "0:50:00"
nToySB = "5000";
nToyB = "5000";

runjobs = 0

if len(sys.argv)>1 : 
    runjobs = int(sys.argv[1])
    seed = runjobs


def drange(start, stop, step):
    r = start
    while r <= stop:
        yield r
        r += step

def drangex(start, stop, step):
    if start <=0 or step<=0: 
        print "need positive values"
        sys.exit(0)

    if step <= 1: 
        print "step must be > 1"
        sys.exit(0)

    r = start
    while r <= stop:
        yield r
        r *= step

def str_rangex(start, stop, step):
    tmp = ""
    for x in drangex(start, stop, step):
        tmp+=str(x)
        tmp+=" "

    return tmp


def substr_rangex(start, stop, step, xmin, xmax):
    if start <=0 or step<=0: 
        print "need positive values"
        sys.exit(0)

    if step <= 1: 
        print "step must be > 1"
        sys.exit(0)

    z = ""
    r = start 
    while r<=stop:
        if r>=xmin and r<=xmax:
            z+=str(r)
            z+=" "

        r *= step

    return z

z = substr_rangex(1, 100, 1.05, 30, 40)
print z

#sys.exit(0)


MassRangeHZZ2l2nu = [x for x in drange(250, 290, 2)]
MassRangeHZZ2l2nu.extend(x for x in drange(290, 350, 5))
MassRangeHZZ2l2nu.extend(x for x in drange(350, 400, 10))
MassRangeHZZ2l2nu.extend(x for x in drange(400, 600, 20))
MassRangeHZZ2l2nu=sorted(set(MassRangeHZZ2l2nu))

for m in MassRangeHZZ2l2nu:
    #if(m<290 or m>420): continue
    mH = str(m)
    sCFGDirectory = os.getcwd()+"/"+mH+"/"
    lands = "/scratch/hpc/mschen/UserCode/mschen/LandS/test/lands.exe -L /scratch/ufhpc/mschen/CMSSW_4_1_6/lib/slc5_amd64_gcc434/libHiggsAnalysisCombinedLimit.so -m "+mH+" --minuitSTRATEGY 0 --maximumFunctionCallsInAFit 50000 "

    #vRrange = ""

    if m== 250 : vRrange=substr_rangex(0.1, 10, 1.1, 0.147604020492 , 0.268029919436 )
    if m== 252 : vRrange=substr_rangex(0.1, 10, 1.1, 0.137229637446 , 0.267980039649 )
    if m== 254 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0995146066933 , 0.240690621882 )
    if m== 256 : vRrange=substr_rangex(0.1, 10, 1.1, 0.11536421857 , 0.260659199649 )
    if m== 258 : vRrange=substr_rangex(0.1, 10, 1.1, 0.115938380152 , 0.242217613148 )
    if m== 260 : vRrange=substr_rangex(0.1, 10, 1.1, 0.115784922032 , 0.253790685583 )
    if m== 262 : vRrange=substr_rangex(0.1, 10, 1.1, 0.1226382697 , 0.232562700581 )
    if m== 264 : vRrange=substr_rangex(0.1, 10, 1.1, 0.144400955263 , 0.25118599184 )
    if m== 266 : vRrange=substr_rangex(0.1, 10, 1.1, 0.133101951524 , 0.242036892062 )
    if m== 268 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0 , 0.21927784607 )
    if m== 270 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0234766558111 , 0.214902398741 )
    if m== 272 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0787030130244 , 0.210018297327 )
    if m== 274 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0744241347809 , 0.218249172524 )
    if m== 276 : vRrange=substr_rangex(0.1, 10, 1.1, 0.120984804123 , 0.205629214467 )
    if m== 278 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0433454324753 , 0.198133392553 )
    if m== 280 : vRrange=substr_rangex(0.1, 10, 1.1, 0.111306154346 , 0.209403077013 )
    if m== 282 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0957605195999 , 0.205463316338 )
    if m== 284 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0882997692942 , 0.201225743351 )
    if m== 286 : vRrange=substr_rangex(0.1, 10, 1.1, 0.104660220043 , 0.199384977722 )
    if m== 288 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0930008116283 , 0.194545160364 )
    if m== 290 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0999365295519 , 0.197349192274 )
    if m== 295 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0912652351853 , 0.2371227939 )
    if m== 300 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0137981487776 , 0.3504865151 )
    if m== 305 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0538633197024 , 0.486194854 )
    if m== 310 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0941259932623 , 0.3810729143 )
    if m== 315 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0568271391654 , 0.41644864995 )
    if m== 320 : vRrange=substr_rangex(0.1, 10, 1.1, 0.073554125708 , 0.26638696758 )
    if m== 325 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0 , 0.4336629 )
    if m== 330 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0 , 0.4317117959789 )
    if m== 335 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0 , 0.4306865139998 )
    if m== 340 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0 , 0.43127668170849 )
    if m== 345 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0591604246947 , 0.3117608967643 )
    if m== 350 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0 , 0.37637791711 )
    if m== 360 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0 , 0.52435813316 )
    if m== 370 : vRrange=substr_rangex(0.1, 10, 1.1, .00419683866 , 0.4637123893 )
    if m== 380 : vRrange=substr_rangex(0.1, 10, 1.1, 0.089289842944 , 0.42386939 )
    if m== 390 : vRrange=substr_rangex(0.1, 10, 1.1, 0.0807690674087 , 0.30314273174 )
    if m== 400 : vRrange=substr_rangex(0.1, 10, 1.1, 0.9105869833156 , 0.306616705596 )
    if m== 420 : vRrange=substr_rangex(0.2, 20, 1.1, 0.01166414934994 , 0.346145308625 )
    if m== 440 : vRrange=substr_rangex(0.2, 20, 1.1, 0.0 , 0.285653865385 )
    if m== 460 : vRrange=substr_rangex(0.2, 20, 1.1, 0.0117452092396 , 0.382746358838 )
    if m== 480 : vRrange=substr_rangex(0.2, 20, 1.1, 0.182811974341 , 0.397330267401 )
    if m== 500 : vRrange=substr_rangex(0.2, 20, 1.1, 0.167460629325 , 0.440604023495 )
    if m== 520 : vRrange=substr_rangex(0.2, 20, 1.1, 0.32099489527 , 0.592049343064 )
    if m== 540 : vRrange=substr_rangex(0.2, 20, 1.1, 0.347845948819 , 0.621224511089 )
    if m== 560 : vRrange=substr_rangex(0.2, 20, 1.1, 0.424896013885 , 0.752336161652 )
    if m== 580 : vRrange=substr_rangex(0.2, 20, 1.1, 0.593282917301 , 0.858647403881 )
    if m== 600 : vRrange=substr_rangex(0.2, 20, 1.1, 0.677728637322 , 0.951642241311 )

    vRrange = "[0.02,0.3,x1.1]"
    if m>400: vRrange="[0.03,0.3,x1.1]"

    #vRrange = "[0.1,2,x1.1]"
    #if m>400: vRrange="[0.2,4,x1.1]"

    if vRrange == "": 
        print m, " no need to run more toys"
        continue

    print vRrange

    outFile=jobname+mH+"_Seed"+str(seed)+".root"
    combination = lands+" -d "+datacards+" -M Hybrid --freq --bNotCalcCLssbb --bSaveM2lnQ --nToysForCLsb "+nToySB+" --nToysForCLb "+nToyB+" --scanRs 1 -vR "+vRrange+" -n "+outFile+"  --seed "+str(seed)
    #/scratch/hpc/mschen/UserCode/mschen/LandS/test/lands.exe -M Hybrid --freq  --nToysForCLsb 3000 --nToysForCLb 1000  --scanRs 1 -d hzz2l2nu*txt -vR [0.4,20,x1.05]

    print combination

    sPBSfile = sCFGDirectory+"/"+outFile+".jdf"
    fPBSfile = open(sPBSfile, "w") 
    fPBSfile.write('#! /bin/bash \n')
    lines="# \n\
            #PBS -r n \n\
            ####PBS -N combination \n\
            #PBS -N dummy \n\
            #PBS -o "+sPBSfile+".out \n\
            #PBS -e "+sPBSfile+".err \n\
            #PBS -m abe \n\
            ######nomail hi i don't want mail notice again, you -M chenmingshui@gmail.com \n\
            #PBS -l walltime="+sEstimateTimePerJob+" \n\
            #PBS -l nodes=1:ppn=1 \n\
            #PBS -l pmem=800mb\n\
            echo \"Job running on `hostname` at  `date` \" \n\
            source /scratch/ufhpc/mschen/CMSSW_416.sh\n\
            cd "+sCFGDirectory+"\n\
            "+combination+" >& "+outFile+".log \n\
            touch "+sPBSfile+".out \n\
            mkdir "+sCFGDirectory+"/"+jobname+" \n\
            mv "+outFile+"* "+sCFGDirectory+"/"+jobname+" \n\
            echo \"Job done on `hostname` at  `date` \" \n\
            \n" 
    fPBSfile.write(lines)               
    fPBSfile.close()
    if runjobs: os.system("qsub "+sPBSfile)
    else: 
        targetFile = ""

        # get a list of *_m2lnQ.root files
        listfiles = glob.glob(sCFGDirectory+"/"+jobname+"/*.root")

        # if size = 0, then do nothing
        if len(listfiles) ==0 : 
            print listfiles, "empty list, do nothing " 
            continue

        # if size = 1, then do not hadd
        if len(listfiles) ==1 : targetFile = listfiles[0]

        # if size > 1, then do hadd,  rename original files with postfix "_bak"
        if len(listfiles) > 1 : 
            targetFile = listfiles[0]+"_merged.root"
            command = "hadd "+targetFile+" "
            for l in listfiles:
                command+=(l+" ")

            print command
            os.system(command)
            for l in listfiles:
                command="mv "+l+" "+l+"_bak"
                os.system(command)
        
        print targetFile
        #ROOT.run(sCFGDirectory+"/"+jobname+"/"+outFile+"_m2lnQ.root", "bands_"+mH, "bands_"+mH, m)
        ROOT.run(targetFile, "bands_"+mH, "bands_"+mH, m)
        if m<=400: print "LOWERSIDE: if m==", mH, ": vRrange=substr_rangex(0.1, 10, 1.1,", ROOT._m2s,",", ROOT._m1s,")"
        else: print "LOWERSIDE: if m==", mH, ": vRrange=substr_rangex(0.2, 20, 1.1,", ROOT._m2s,",", ROOT._m1s, ")"
