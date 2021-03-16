#! /usr/bin/env python

# TO DO
# 1. 
#
import sys,os,math,argparse
sys.path.append( '..' )
from x4i import exfor_manager, exfor_entry
from nuclear import *
from PoPs import database as databaseModule

# BAD DATA:
rescale_e3 = ['20920']  # Scobel data should be label b, not mb

excludeSE = ['F0005002','10547002','11182002','11164003','10755003']

# excludeSE += ['F0397006']  # EN-MIN but no EN-MAX




lightnuclei = {'n':'N', '2n':'2N', 'H1':'P', 'H2':'D', 'H3':'T', 'He3':'HE3', 'He4':'A', 'photon':'G'}
GNDSnuclei = {'N':'n', '2N':'2n', 'P':'H1', 'D':'H2', 'T':'H3', 'HE3':'He3', 'A':'He4', 'G':'photon'}
lightA      = {'photon':0, 'n':1, '2n':2, 'H1':1, 'H2':2, 'H3':3, 'He3':3, 'He4':4}
lightZ      = {'photon':0, 'n':0, '2n':0, 'H1':1, 'H2':1, 'H3':1, 'He3':2, 'He4':2}


defaultPops = '../ripl2pops_to_Z8.xml'
debug = False
print('Command:',' '.join(sys.argv[:]) ,'\n')

parser = argparse.ArgumentParser(description='Prepare data for Rflow')

parser.add_argument('CN', type=str, help='compound-nucleus   e.g. N-15 or N15.')
parser.add_argument("-E", "--EnergyMax", type=float,  default="19.9", help="Maximum lab incident energy to allow. Default=19.9")
parser.add_argument("-e", "--EnergyMin", type=float,  default="0.1", help="Minimum lab incident energy to allow. Default=0.001")
parser.add_argument("-p", "--projectiles", type=str, default='any', help="Specific projectile only or default 'any'")
parser.add_argument("-n", "--nat", type=str,  help="CN which can be replaced by element as being isotopically dominant (e.g. N-14).")
parser.add_argument("-d", "--dir", type=str, default="Data_X4s", help="Output directory.")

parser.add_argument("-i", "--include", metavar="INCL",  nargs="*", help="Subentries to include, if not all")
parser.add_argument("-x", "--exclude", metavar="EXCL", nargs="*", help="Subentries to exclude if any string within subentry name")
parser.add_argument(      "--pops", type=str, default=defaultPops, help="pops files with all the level information from RIPL3. Default = %s" % defaultPops)
parser.add_argument(      "--allowNeg", action="store_true", help="Allow points if lower error bar goes negative")
parser.add_argument("-t", "--tolerance", type=float,  default="0.1", help="Print warnings for level discrepancies larger than this")


args = parser.parse_args()


pops = databaseModule.database.readFile( args.pops )
CN = args.CN
ElabMin = args.EnergyMin
ElabMax = args.EnergyMax
proj = args.projectiles
nat = args.nat

# print(list(elementsSymbolZ.keys()))

allObs = ['CS','DA','DA/DE','DE']
plots1d = ['CS','NU']
plots2d = ['DA','NU/DE','DA-G']
classes = []
classes = ['INL','NON']
No_neg_errorbars = not args.allowNeg

exclude = args.exclude
include = args.include
if include:
    print("    Include only ",include)
else:
    include = [None]

pi = math.pi
rad = 180/pi
amu    = 931.494013

def getmz(nucl):
    if nucl[0]=='2':
        multiplicity = 2 
        nucl = nucl[1:]
    else: 
        multiplicity = 1
        
    n = pops[nucl.replace('-','')]
    mass = n.getMass('amu')
    if nucl=='n' or nucl=='photon':
        charge = 0
    else:
        charge = n.nucleus.charge[0].value
    return (mass,charge)
    
def PoPsLevelFind(nucname,level,trace):
    if level == 0.0: return ('')
    proximity = 1e6
    nearest = 0
    Enearest = 1e6
    nucl = nucname.replace('-','')
    for l in range(30):
        tex = nucl + ('_e%s' % int(l) if l > 0 else '')
        try:
            E = pops[tex].energy[0].pqu('MeV').value
            if abs(E-level) < proximity:
                proximity = abs(E-level)
                Enearest = E
                nearest = l
                if proximity < 1e-3: break
        except:
            pass
    tag = '_e%s' % nearest if nearest > 0 else ''
    if trace and proximity>args.tolerance : print('Level nearest to E=',level,'in',nucname,'is at',Enearest,':',tag,'missing by%6.3f' % proximity)
#     print('Level nearest to E=',level,'in',nucname,'is at',Enearest,':',tag,'missing by%6.3f' % proximity)
    return(tag)

if __name__ == "__main__":

    if '-' in CN:
        name,Acn = CN.split('-')
    else:
        name = CN[:2] if CN[1].isalpha() else CN[0]
        Acn = CN[len(name):]
        CN = '%s-%s' % (name,Acn)
    gndsName = name+Acn
    CNmass,CNcharge = getmz(gndsName)

    Acn = int(Acn)
    Zcn = elementsSymbolZ[name]

    channels = []
    for p in lightA.keys():
        Ap = lightA[p]
        Zp = lightZ[p]
        At = Acn - Ap
        Zt = Zcn - Zp
        if Ap<=At:
            t = elementsZSymbolName[Zt][0] + '-%i' % At
            channels.insert(0,(p,t))
    print('Binary channels for',CN,':\n',channels)
    
    db = exfor_manager.X4DBManagerPlainFS( )  
    dir = args.dir
    #if nat is not None: dir = 'nat'+name
    os.system('mkdir '+dir)
    excuses = {}

    searched = ''
    levelsSeen = []
    partitions = []
    subentries = []
    for In in channels:
        p,t = In
        pn = lightnuclei[p]   # x4 name
        if pn == '2N': continue  # not a physical projectile
        if pn == 'G': continue  # not a spectroscopic projectile
        if proj != 'any' and pn != proj: continue
        cl = ['TOT'] + classes if p=='n' else classes
        more = [(c,'.') for c in cl]
        pmass,pZ = getmz(p)
        tmass,tZ = getmz(t)
        if len(partitions) == 0:
            qvalue = 0.0
            Tmass = pmass + tmass

        Partmass = pmass + tmass
        qvalue = (Tmass - Partmass) * amu  # relative to the first partition
                
        targets = [t]
        natT = ''
        if nat is not None:  
            if t==nat:    # replace t by nat-t as isotopically dominant
                natT,natA = t.split('-')
                tn = natT + '-0'
                print('Also try replacing',t,'by',tn)
                targets.append(tn)
#             print('Targets=',targets)
        
        for t in targets: 

            this_partition_Listed = False
            for Out in more+channels:
                po,to = Out
                Qvalue_to_gs = 0 # if Out == 'TOT' else -99
                if Out in channels:
                    emass,poZ = getmz(po)
                    rmass,toZ = getmz(to)
                    Qvalue_to_gs = (Partmass - (emass + rmass)) * amu # relative to this reactions incoming channel

                pon = po if po in cl else lightnuclei[po]   # x4 name
                if po == p: pon = 'EL'
                residual = to.replace('-','')

                reaction = pn + ',' + pon
                projectile,ejectile = (pn,pon)
                if ejectile == 'EL': ejectile=projectile
        
#                 if ejectile == 'N': continue
#                 if ejectile == 'G': continue
                if ejectile == 'NON': continue
                if ejectile == 'INL': residual = t.replace('-','')

                for obs in allObs:
                    searchX4 = "\n ***** Search X4 for %s(%s), %s with Q2gs=%7.3f MeV" % (t,reaction,obs,Qvalue_to_gs)
                    print("\n", searchX4)
                    exit_quantity = 'angle' if obs=='DA' else 'exit energy'
                    
                    print("    Include only ",include)
                    subents = {}
                    for incl in include:
                        subs = db.retrieve( target = t, reaction = reaction, quantity = obs , SUBENT = incl)
                        for e in list(subs.keys()) : subents[e] = subs[e]

                    searchesFound = 0

                    if obs == 'DA/DE': continue
                        
                        
                    for e in subents:
                        print('    Entry:',e)
                        try:
                            if isinstance( subents[ e ], exfor_entry.X4Entry ): 
                                ds = subents[ e ].getDataSets( )
                        except KeyError:
                            print("Got KeyError on entry", e, "ignoring")
                            print(list(subents[ e ].keys()))
                            continue
                        for d in ds:
                            print()
                            npoints = 0
                            relErrMin = 1e6
                            relErrMax = 0
                            result = str( ds[ d ] ) 
                            Authors = ', '.join(ds[d].author[:])
                            author1 =  ds[d].author[0].split('.')[-1]
                            Reaction = str(ds[d].reaction[0])
                            Reaction1= str(ds[d].reaction[1])
                            try:
                                frame = '' if obs != 'DA' else ds[d].referenceFrame
                            except:
                                frame = ''
                                print('Assume frame =',frame,'.')
                            ratioToRutherford =  'Rutherford' in Reaction
                            ratioXS = None
                            if 'esonance' in Reaction: continue
                            if 'pin' in Reaction: continue
                            print('       ',d,'from',author1,'('+ds[d].year+') for',Reaction)
                            subent = d[1]
                            ent = subent[:5]
                            x4file = ent + '.x4'
                            n = ds[d].numrows()
                            S_factor = 'SFC' in str(ds[d].reaction[:]) or 'S-factor' in Reaction
                            Legendre = 'Legendre' in Reaction
                            shape_data = None
                        
                            if subent in excludeSE:
                                excuses[subent] = 'Subentry on global excluded list'
                                print(20*' ',excuses[subent])
                                continue
                            # if exclude is not None and subent in exclude:
                            if exclude is not None and any(sub in subent for sub in exclude):
                                excuses[subent] = '## Data set in local exclusion list'
                                print(20*' ',excuses[subent])
                                continue
                        
                            if n < 3: continue  # need 3 or more data points 
                        
                            if len(ds[d].reaction)>1 and ('::' in Reaction1 or ')/(' in Reaction1 or '//' in Reaction1 or 'ratio' in Reaction1):
                                excuses[subent] = 'only ratio data  1 ' + Reaction1
                                print(20*' ',excuses[subent])
                                continue

                            if 'Legendre' in Reaction:
                                excuses[subent] = 'Legendre data not yet extracted'
                                print(20*' ',excuses[subent])
                                continue
                            if '+' in Reaction:
                                excuses[subent] = 'Summed reactions not yet processed'
                                print(20*' ',excuses[subent])
                                continue
                            if 'Delayed' in Reaction:
                                excuses[subent] = 'Delayed fission yields not yet processed'
                                print(20*' ',excuses[subent])
                                continue
                            if 'Cosine' in Reaction:
                                excuses[subent] = 'Cosine data not yet extracted'
                                print(20*' ',excuses[subent])
                                continue
                            if 'resonance' in Reaction or 'Resonance' in Reaction or 'Momentum' in Reaction or 'Spin' in Reaction:
                                excuses[subent] = 'resonance properties not yet extracted'
                                print(20*' ',excuses[subent])
                                continue
                            if 'isotopic' in Reaction:
                                excuses[subent] = 'isotopic data not yet distinguished'
                                print(20*' ',excuses[subent])
                                continue
#                             if 'S-factor' in Reaction:
#                                 excuses[subent] = 'S-factor data not yet plotted'
#                                 print(20*' ',excuses[subent])
#                                 continue
#                             if 'SFC' in str(ds[d].reaction[:]):
#                                 excuses[subent] = 'S-factor data not yet plotted'
#                                 print(20*' ',excuses[subent])
#                                 continue
                            if 'MXW' in str(ds[d].reaction[:]):
                                excuses[subent] = 'Maxwellian data not processed'
                                print(20*' ',excuses[subent])
                                continue
                            if 'Vector' in Reaction:
                                excuses[subent] = 'polarization data not yet plotted'
                                print(20*' ',excuses[subent])
                
                            #if 'Relative' in Reaction or ')/(' in Reaction or '//' in Reaction or 'ratio' in Reaction:
                            if ')/(' in Reaction:
    #                           measure,reference = Reaction.split(')/(')
    #                           ratioXS = None
    #                           for p_std,t_std in standards.keys():
    #                               if p != p_std: continue                    
    #                               try:
    #                                   refChannel = reference.split(t_std+'('+p+',')[1].split(')')[0]
    #                                   print(15*' ','Ratio to %s+%s -> %s' %(p,t_std,refChannel))
    #                                   ratioXS = standards[(p_std,t_std)][0].getReaction(refChannel).crossSection
    #                                   break
    #                               except:
    #                                   continue
                                if ratioXS is None:
                                    excuses[subent] = 'only relative data 0 ' + Reaction
                                    print(20*' ',excuses[subent])
                                    continue
                            if shape_data is None:
                                if 'Relative data' in Reaction:
                                    shape_data = True
                                else:
                                    shape_data = False

                            labels = ds[d].labels
                            units  = ds[d].units 
                            E_index = None
                            included = False
                            emin = 1e3
                            emax = 0
                            energy = True
                            Ein = 'lab' 
                            E_index = E_index_min =  E_index_max = None
                            if 'Energy' in labels:
                                E_index = labels.index('Energy')
                            elif 'EN-CM' in labels:
                                E_index = labels.index('EN-CM')
                                E_scale = 1. # (pMass+tMass)/tMass  # convert to projectile lab energy for comparisons
                                Ein = 'cm' 
                            elif 'EN' in labels:
                                E_index = labels.index('EN')
                            elif 'EN-APRX' in labels:
                                E_index = labels.index('EN-APRX')
                            elif 'KT' in labels:
                                E_index = labels.index('KT')
                            elif 'EN-MEAN' in labels:
                                E_index = labels.index('EN-MEAN')
                            elif 'EN-DUMMY' in labels:
                                E_index = labels.index('EN-DUMMY')
                            elif 'EN-MIN' in labels and 'EN-MIN' in labels:
                                E_index_min = labels.index('EN-MIN')
                                try:
                                    E_index_max = labels.index('EN-MAX')
                                except:
                                    already = subent in excuses.keys()
                                    excuses[subent] = 'EN-MIN but no EN-MAX'
                                    if not already: print(20*' ',excuses[subent])
                                    continue
                            elif 'KT-DUMMY' in labels:
                                E_index = labels.index('KT-DUMMY')
                            if E_index is not None:
                                E_units = units[E_index]
                                E_scale = 1.
                                if E_units.lower() == 'milli-ev':  E_scale *= 1e-9  # convert to MeV
                                if E_units.lower() == 'ev':  E_scale *= 1e-6
                                if E_units.lower() == 'kev': E_scale *= 1e-3
                                if E_units.lower() == 'gev': E_scale *= 1e+3
                                if E_units.lower() == 'gev/c': E_scale *= 1e+3
        
                                for stuff in ds[d].data:
                                    E = stuff[E_index]
                                    if E is None: continue
                                    E *= E_scale
                                    emin = min(emin,E)
                                    emax = max(emax,E)
                                errs = ' '
                                included = True
        
                            elif 'EN-MIN' in labels in labels:
                                E_index_min = labels.index('EN-MIN')
                                E_index_max = labels.index('EN-MAX')
                                for stuff in ds[d].data:
                                    Emin_in = stuff[E_index_min]
                                    Emax_in = stuff[E_index_max]
                                    if Emin_in is None: continue
                                    E = (Emin_in + Emax_in)*0.5
                                    E *= E_scale
                                    emin = min(emin,E)
                                    emax = max(emax,E)
                                errs = ' '
                            elif 'MOM' in labels:
                                E_index = labels.index('MOM')
                                energy = False
                                print('MOMentum index at',E_index)
                            else:
                                emin = ' '
                                emax = ' '
                                errs = ' '
                                included = True
                        
                            if not included: continue
                            
                            level_index = None; levelnum_index = None; Q_index = None
                            if 'E-LVL' in labels:
                                level_index = labels.index('E-LVL')
                            if 'E-EXC' in labels:
                                level_index = labels.index('E-EXC')
                            if 'E' in labels:
                                level_index = labels.index('E')
                            if 'LVL-NUMB' in labels:
                                levelnum_index = labels.index('LVL-NUMB')
                            if 'Q-VAL' in labels:
                                Q_index = labels.index('Q-VAL')
                            levelnum_g = None
                            if 'GROUND STATE'  in Reaction1: levelnum_g = 0
                            if 'FIRST EXCITED' in Reaction1: levelnum_g = 1
                            if levelnum_g is not None: print("        Dataset for final level #",levelnum_g,'from "'+Reaction1+'"')
                                
# FIND RESIDUAL LEVELS, AND SUBDIVIDE DATA OUTPUTS:
                            levels = set()
                            leveltags = set()       
                            for stuff in ds[d].data:
                            
                                E = 0.
                                if E_index is not None:
                                    E = stuff[E_index] 
                                else:
                                    if E_index_min is None and E_index_max is None: 
                                        continue
                                    Emin_in = stuff[E_index_min] if E_index_min is not None else 0.0
                                    Emax_in = stuff[E_index_max]
                                    E = (Emin_in + Emax_in)*0.5
                                if E is None: continue
                                if not (ElabMin < E*E_scale < ElabMax): continue
                                
                                if level_index is not None and obs in ['DA','CS'] and stuff[level_index] is not None:
                                    level = stuff[level_index]
                                    if level is not None:
                                        levelUnits = units[level_index]
                                        if levelUnits.lower() == 'ev':  level *= 1e-6
                                        if levelUnits.lower() == 'kev': level *= 1e-3
                                        levels.add(level)
#                                         print('level=',level,'from level_index=',level_index)
                                    
                                elif levelnum_index is not None and obs in ['DA','CS'] and  stuff[levelnum_index] is not None:
                                    levelnum = stuff[levelnum_index]
                                    if levelnum is not None:
                                        tex = residual + ('_e%s' % int(levelnum) if levelnum > 0 else '')
                                        try:
                                            level = float(pops[tex].energy[0].pqu('MeV').value)
#                                             print('level=',level,'from levelnum=',levelnum)
                                        except:
                                            excuses[subent] = "Residual %s not found in evaluation pops" % tex
                                            print(20*' ',excuses[subent])
                                            continue
                                        levels.add(level)
                                elif levelnum_g is not None:
                                    tex = residual + ('_e%s' % int(levelnum_g) if levelnum_g > 0 else '')
                                    try:
                                        level = pops[tex].energy[0].pqu('MeV').value
                                        levels.add(level)
#                                         print('level=',level,'from levelnum_g=',levelnum_g)
                                    except:
                                        print('Energy of state',tex,'NOT FOUND IN POPS. ERROR')
                                        level = levelnum_g*1.11111
                                        levels.add(level)
               
                                elif Q_index is not None:
                                    Q = stuff[Q_index]
                                    level = Qvalue_to_gs - Q
#                                     print(' level =',level,'from Qgs-Q =',Qvalue_to_gs ,'-', Q)
                                    levels.add(level)
                                else:
                                    level = 0.
                                
                                leveltag = PoPsLevelFind(residual,level,True)   
                                leveltags.add(leveltag)
#                             print('Found leveltags:',leveltags,'from',levels,'for',subent)
# GET DATA                            
                            if 'EN-RSL' in labels:
                                dE_index = labels.index('EN-RSL')
                            elif 'EN-ERR' in labels:
                                dE_index = labels.index('EN-ERR')
                            else:
                                dE_index = None
                            dE_scale = E_scale  # any change in units for energy errors?
                            dE_ISscale = -1
                            if dE_index is not None: 
                                dE_units = units[dE_index] 
                                if dE_units.lower() == 'milli-ev':  dE_scale *= 1e-9
                                if dE_units.lower() == 'ev':  dE_scale *= 1e-6
                                if dE_units.lower() == 'kev': dE_scale *= 1e-3
                                if dE_units.lower() == 'gev': dE_scale *= 1e+3
                                if dE_units.lower() == 'microsec/m': dE_ISscale = 1e-6 # to convert to sec/m
                                if dE_units.lower() == 'nsec/m': dE_ISscale = 1e-9 # to convert to sec/m
                
                         
                            if 'Data' in labels:
                                Data_index = labels.index('Data')
                            elif 'DATA' in labels:
                                Data_index = labels.index('DATA')
                            elif 'DATA-CM' in labels:
                                Data_index = labels.index('DATA-CM')
                            else:
                                excuses[subent] = 'No data column found'
                                print(20*' ',excuses[subent])
                                continue

                            Data_units = units[Data_index]
                            data_scale = 1.0
                            RT_data_scale = None
                            if debug: print('Data_index,Data_units:',Data_index,Data_units)
                            if Data_units.lower() in ['no-dim','prt/fis']: data_scale = 1.
                            if Data_units.lower() in ['arb-units']: 
                                shape_data = True
                                data_scale = 1.
                            if Data_units.lower() in ['mb','mb/sr']: data_scale = 1e-3 
                            if Data_units.lower() in ['micro-b','micro-b/sr','mu-b/sr']: data_scale = 1e-6
                            if Data_units.lower() in ['nb','nb/sr']: data_scale = 1e-9 
                            if Data_units.lower() in ['1/mev']: data_scale = 1
                            if Data_units.lower() in ['1/mev']: data_scale = 1
                            if Data_units.lower() in ['b*rt-ev']: RT_data_scale = 1e6  # divide data by sqrt(E*RT_data_scale) when E is MeV
                            if Data_units.lower() in ['no-dim'] and ratioXS is None and not ratioToRutherford: 
                                excuses[subent] = 'Only dimensionless data'
                                print(20*' ',excuses[subent])
                                continue
                            if 'times 4 pi' in Reaction: 
                                print(10*' ',' Has "times 4 pi", so divide by 4pi')
                                data_scale /= (4*pi)
                            if e in rescale_e3: data_scale *= 1e3
                        
                            if S_factor:
                                if Data_units.lower() in ['mb*mev']: data_scale = 1e-3
                                if Data_units.lower() in ['b*mev']: data_scale = 1.0
                                if Data_units.lower() in ['b*kev']: data_scale = 1e-3
                                if Data_units.lower() in ['mb*kev']: data_scale = 1e-6
                                print('S_factor: data_scale=',data_scale)
           
                            dData_index = None
                            dData_index1 = None
                            dData_index2 = None
                            if 'd(Data)' in labels:
                                dData_index = labels.index('d(Data)')
                            elif 'DATA-ERR' in labels:
                                dData_index = labels.index('DATA-ERR')
                            elif 'ERR-T' in labels:
                                dData_index = labels.index('ERR-T')
                            elif 'ERR-S' in labels:
                                dData_index = labels.index('ERR-S')
                            else:
                                dData_index = None
                            
                            if 'DATA-ERR1' in labels:
                                dData_index1 = labels.index('DATA-ERR1')
                            if 'ERR-1' in labels:
                                dData_index1 = labels.index('ERR-1')
                            if 'ERR-T' in labels:
                                dData_index2 = labels.index('ERR-T')

                            if dData_index is not None and units[dData_index].lower() in ['per-cent']: 
                                dataErr_rel = 0.01
                                dataErr_abs = 0.
                            else:
                                dataErr_rel = 0; 
                                dataErr_abs = 1
                            if dData_index1 is not None and units[dData_index1].lower() in ['per-cent']: 
                                dataErr_rel1 = 0.01
                                dataErr_abs1 = 0.
                            else:
                                dataErr_rel1 = 0; 
                                dataErr_abs1 = 1
                            if dData_index2 is not None and units[dData_index2].lower() in ['per-cent']: 
                                dataErr_rel2 = 0.01
                                dataErr_abs2 = 0.
                            else:
                                dataErr_rel2 = 0; 
                                dataErr_abs2 = 1
    #                         print('Err @',dData_index,dData_index1,dData_index2,'so % are',dataErr_rel,dataErr_rel1,dataErr_rel)

                            ddata_scale = data_scale
                            if dData_index is not None and units[dData_index] != Data_units:
                                dData_units = units[dData_index]
                                if debug: print('dData_index,dData_units:',dData_index,dData_units)
                                if dData_units.lower() in ['no-dim','prt/fis','arb-units']: ddata_scale = 1.
                                if dData_units.lower() in ['mb','mb/sr']: ddata_scale = 1e-3
                                if dData_units.lower() in ['micro-b','micro-b/sr','mu-b/sr']: ddata_scale = 1e-6
                                if dData_units.lower() in ['nb','nb/sr']: ddata_scale = 1e-9
                                if dData_units.lower() in ['1/mev']: ddata_scale = 1
                                if 'times 4 pi' in Reaction:
                                    # print(10*' ',' Has "times 4 pi", so divide by 4pi')
                                    ddata_scale /= (4*pi)
                                if e in rescale_e3: ddata_scale *= 1e3

                            if obs in plots2d:
                                cosines = False
                                if 'Angle' in labels:
                                    Ang_index = labels.index('Angle')
                                elif 'ANG' in labels:
                                    Ang_index = labels.index('ANG')
                                elif 'ANG-CM' in labels:
                                    Ang_index = labels.index('ANG-CM')
                                elif 'COS' in labels:
                                    Ang_index = labels.index('COS')
                                    cosines = True   # data already in cosines
                                elif 'COS-CM' in labels:
                                    Ang_index = labels.index('COS-CM')
                                    cosines = True   # data already in cosines
#                                 elif 'E' in labels and observable!='DA':
#                                     Ang_index = labels.index('E')
                                else:
                                    excuses[subent] = 'No exit '+exit_quantity+' column found'
                                    print(20*' ',excuses[subent])
                                    continue
                                if 'd(Angle)' in labels:
                                    dAng_index = labels.index('d(Angle)')
                                elif 'ANG-ERR-D' in labels:
                                    dAng_index = labels.index('ANG-ERR-D')
                                elif 'E-ERR-DIG' in labels and observable!='DA':
                                    dAng_index = labels.index('E-ERR-DIG')
                                else:
                                    dAng_index = None
                            else:
                                Ang_index = dAng_index = None
                            
                            
                            for leveltag in leveltags:

                        # FIND IF ANY DATA IS IN RANGE, before creating a file! 
                                Datafound = 0
                                for stuff in ds[d].data:
                                
                                    E = None
                                    if E_index is not None:
                                        E = stuff[E_index] 
                                    else:
                                        if E_index_min is None and E_index_max is None: 
                                            continue
                                        Emin_in = stuff[E_index_min] if E_index_min is not None else 0.0
                                        Emax_in = stuff[E_index_max]
                                        E = (Emin_in + Emax_in)*0.5
                                    if E is None: continue
                                    if not (ElabMin < E*E_scale < ElabMax): continue

                                    if E is None: continue
                                    
                                    if level_index is not None and obs in ['DA','CS']  and stuff[level_index] is not None:
                                        level = stuff[level_index]
                                        levelUnits = units[level_index]
                                        if levelUnits.lower() == 'ev':  level *= 1e-6
                                        if levelUnits.lower() == 'kev': level *= 1e-3
                                    elif levelnum_index is not None and obs in ['DA','CS']  and stuff[levelnum_index] is not None:
                                        levelnum = stuff[levelnum_index]
                                        if levelnum is not None:
                                            tex = residual + ('_e%s' % int(levelnum) if levelnum > 0 else '')
                                            try:
                                                level = pops[tex].energy[0].pqu('MeV').value
                                            except:
                                                excuses[subent] = "Residual %s not found in evaluation pops" % tex
                                                print(20*' ',excuses[subent])
                                                continue
                                    elif levelnum_g is not None:
                                        tex = residual + ('_e%s' % int(levelnum_g) if levelnum_g > 0 else '')
                                        try:
                                            level = pops[tex].energy[0].pqu('MeV').value
                                        except:
                                            level = 0
                                    elif Q_index is not None:
                                        Q = stuff[Q_index]
                                        level = Qvalue_to_gs - Q
                                    else:
                                        level = 0.
                                        
                                    if leveltag != PoPsLevelFind(residual,level,False): continue
                                                                    
                                    Datafound += 1                    
                        
                                if Datafound == 0: continue # nothing for this leveltag
                                qual = '-' + d[2] if d[2] != ' ' else ''
                                file = author1 + '-' + subent + qual + leveltag +'.dat'
                                X4_tag = subent + qual
                                file = file.replace(' ','_')
                                data_output = open(dir + '/' + file,'w')
        
                                described = False
                                searchesFound += 1
                                for stuff in ds[d].data:
                                    dData = dData1 = dData2 = 0.0
                                    E = None

                                    if obs in plots2d:
                                        # E,Ang,Data,dE,dAng,dData = stuff
                                        Ediff = 0.0
                                        if E_index is not None:
                                            E = stuff[E_index] 
                                        else:
                                            if E_index_min is None and E_index_max is None: 
                                                print('Neither E_index_min or E_index_max given. Skip *****')
                                                continue
                                            Emin_in = stuff[E_index_min] if E_index_min is not None else 0.0
                                            Emax_in = stuff[E_index_max]
                                            E = (Emin_in + Emax_in)*0.5
                                            Ediff = (Emax_in - Emin_in)*0.5
                                        Ang = stuff[Ang_index]
                                        if cosines: 
                                            if abs(Ang) <= 1.0:
                                                Ang = math.acos(Ang)*rad
                                            else:
                                                excuses[subent] = 'Data abs(cosine=%s) is larger than 1' % Ang
                                                print(20*' ',excuses[subent])
                                                continue
                                        Data = stuff[Data_index]
                                        dE = stuff[dE_index]     if dE_index   is not None else Ediff
                                        dAng = stuff[dAng_index] if dAng_index is not None else 0.0
                                        dData = stuff[dData_index] if dData_index is not None else 0.0
                                        dData1= stuff[dData_index1] if dData_index1 is not None else 0.
                                        dData2= stuff[dData_index2] if dData_index2 is not None else 0.
                                        if dData2 is None: dData2 = 0. 
                                        if debug: print("E,dE,Ang,dAng,Data,dData:",E,dE,Ang,dAng,Data,dData)
                                    else:  # 'CS'
                                        Ang=None; dAng=0; dE=None
                                        # E,Data,dE,dData = stuff
                                        if E_index is not None:
                                            E = stuff[E_index] 
                                        else:
                                            if E_index_min is None and E_index_max is None: 
                                                print('Neither E_index_min or E_index_max given. Skip *****')
                                                continue
                                            Emin_in = stuff[E_index_min] if E_index_min is not None else 0.0
                                            Emax_in = stuff[E_index_max]
                                            E = (Emax_in + Emin_in)*0.5
                                            dE =(Emax_in - Emin_in)*0.5
                                        Data = stuff[Data_index]
                                        # if debug: print("dE_index,dData_index:",dE_index,dData_index)
                                        dE = stuff[dE_index]       if dE_index    is not None else dE
                                        dData = stuff[dData_index] if dData_index is not None else 0.
                                        dData1= stuff[dData_index1] if dData_index1 is not None else 0.
                                        dData2= stuff[dData_index2] if dData_index2 is not None else 0.
                                        if dE_index is not None and dE is not None and units[dE_index].lower() == 'per-cent': 
                                            dE  = E * dE / 100. 
                                        if debug: print("E,dE,Data,dData:",E,dE,Data,dData,"from dE_index,dData_index:",dE_index,dData_index)
                                    if  E    is None: continue
                                    if  Data is None: continue
                                    if dData is None: continue
                                    if dData1 is None: continue
                                    if dE    is None: dE    = 0
                                    if energy:
                                       E  *= E_scale
                                    else:
                                       if debug: print('P =',E,' and m,amU=',masses_in[0] , amu,'so E',(E*E_scale)**2/ (2 * masses_in[0] * amu))
                                       E = (E*E_scale)**2/ (2 * masses_in[0] * amu)
                                    dE *= dE_scale
        #                            print('dData:',dData,dataErr_abs,Data,dataErr_rel)
                                    dData *=  ( dataErr_abs**2 + (Data*dataErr_rel)**2 ) ** 0.5
#                                     print('dData:',dData1,dataErr_abs1,Data,dataErr_rel1)
                                    dData1 *=  ( dataErr_abs1**2 + (Data*dataErr_rel1)**2 ) ** 0.5
#                                     print('dData2:',dData2,dataErr_abs2,Data,dataErr_rel2)
                                    if dData2 is None: dData2 = 0. 
                                    dData2 *=  ( dataErr_abs2**2 + (Data*dataErr_rel2)**2 ) ** 0.5
                                    if debug: print('err:',dData/Data,dData1/Data,dData2/Data)
                                    dData  = (dData**2 + dData1**2 + dData2**2) ** 0.5
                                    if RT_data_scale is not None: Data /= (E*RT_data_scale) ** 0.5

                                    if E is None: continue
                                    if No_neg_errorbars and (Data < 0. or Data - dData < 0): continue
                                    if not (ElabMin < E < ElabMax): continue

                                    level = 0.
                                    if level_index is not None and obs in ['DA','CS'] and stuff[level_index] is not None:
                                        level = stuff[level_index]
                                        levelUnits = units[level_index]
                                        if levelUnits.lower() == 'ev':  level *= 1e-6
                                        if levelUnits.lower() == 'kev': level *= 1e-3
                                    elif levelnum_index is not None and obs in ['DA','CS']  and stuff[levelnum_index] is not None:
                                        levelnum = stuff[levelnum_index]
                                        if levelnum is not None:
                                            tex = residual + ('_e%s' % int(levelnum) if levelnum > 0 else '')
                                            try:
                                                level = pops[tex].energy[0].pqu('MeV').value
                                            except:
                                                pass
                                    elif levelnum_g is not None:
                                        tex = residual + ('_e%s' % int(levelnum_g) if levelnum_g > 0 else '')
                                        try:
                                            level = pops[tex].energy[0].pqu('MeV').value
                                        except:
                                            level = 0
                                    elif Q_index is not None:
                                        Q = stuff[Q_index]
                                        level = Qvalue_to_gs - Q
                                        
                                    if leveltag != PoPsLevelFind(residual,level,False): continue
                                    
                                    if obs in plots2d:
                                        if dAng  is None: dAng  = 0.
                                        if  Ang  is None: continue
                                    if Ang is None: Ang=-1
                                    if 'Legendre' in Reaction: Ang = -10
                            
                                    if not ratioToRutherford:
                                        Data  = Data  * data_scale * 1000.   # mb
                                        dData = dData * ddata_scale * 1000.  # mb
                                        if not described: print('       data_scale=',data_scale,' so  X4*',data_scale * 1000.,' now')
                                        described = True

                                    print(E, Ang, Data , dData, file=data_output)   # MeV, mb units
                                
                                    relErr = dData/max(Data,1e-6)
                                    relErrMin = min(relErrMin, relErr)
                                    relErrMax = max(relErrMax, relErr)
                                    npoints += 1
                            
                                data_output.close()
                                npts = str(npoints)
                         

                                level = '0'
                                if ejectile == 'INL':
                                    ejectile=projectile
                                    level = '*'
                                if leveltag != '':
                                    level = leveltag[2:]
                                sys_error = '5'  if ejectile != 'TOT' else '1' 
                                if shape_data is not None and shape_data: sys_error = '-1'
                                stat_error= '5' if ejectile != 'TOT' else '2' #  for when pointwise errors are 0.0
                                angle_integrated = 'TRUE' if obs=='CS' else 'FALSE'
                                norm = str(1)
                                group = ''
                    
                                splitnorms = 'FALSE'
                                lab = 'TRUE' if frame == 'Lab' else 'FALSE'
                                abserr = 'TRUE' 
                                scale = 'mb'
                                filedir = dir+'/'
                                if Out in channels:
                                    Aflip = Out in channels and emass > rmass*1.1  # because of set of partitions for R-matrix parameters.
        #                                 print('     Aflip =',Out,'in',channels,'=',Out in channels,'and',po,emass,'>',to,rmass*1.1,'::',Aflip)
        #                             print('     Aflip =',Out in channels,'and',po,emass,'>',to,rmass*1.1,'::',Aflip)
                                else:
                                    Aflip = False
                                eshift = str(0)
                                ecalib = str(0)
                                splitshifts = 'FALSE'
                                ratioRuth = 'TRUE' if ratioToRutherford else 'FALSE'
                                Sfactor = 'TRUE' if S_factor else 'FALSE'
                    
                                projectile = GNDSnuclei.get(projectile,projectile)
                                ejectile = GNDSnuclei.get(ejectile,ejectile)
                                if Aflip:   
                                    projectile,ejectile = ejectile,projectile
                                    pmass,tmass,pZ,tZ = tmass,pmass,tZ,pZ
                                else:
                                    if not this_partition_Listed:
                                        partitions.append([p,t.replace('-',''),str(qvalue),str(pmass),str(tmass),str(pZ),str(tZ)])
                                        print('Added partition',partitions[-1])
                                        this_partition_Listed = True

                           
                                target  = t.replace('-','')  # GNDS name
                                if target == natT+'0': target = natT + natA
                                nucname = target  # TEMPORARY ??
                                Aflip = 'TRUE' if Aflip else 'FALSE'
                                s_relErrMin = '%.2e' % relErrMin
                                s_relErrMax = '%.2e' % relErrMax
                                print('shape_data:',shape_data)
                                print(projectile,ejectile,target,residual,level,file,sys_error,stat_error,angle_integrated,norm,group,splitnorms,lab,abserr,scale,filedir,Aflip,Ein,eshift,ecalib,splitshifts,ratioRuth,Sfactor,npoints,X4_tag,s_relErrMin,s_relErrMax)
                                if included:
        #                             subentries.append(', '.join(['"'+part+'"' for part in [subent,x4file,t,reaction,obs,Reaction,str(emin),str(emax),
        #                                 str(n),errs,Authors,ds[d].year,str(ds[d].reference)]]))
                                    subentries.append(','.join([projectile,ejectile,target,residual,level,file,sys_error,stat_error,angle_integrated,norm,group,splitnorms,lab,abserr,scale,filedir,Aflip,Ein,eshift,ecalib,splitshifts,ratioRuth,Sfactor,npts,X4_tag,s_relErrMin,s_relErrMax]))
                                #                             if 0.0 in levels: levels = levels.remove(0.0)
#                             print('Looked for levels',levels,'in',residual) 
                            if levels is not None:
                                if len(levels)>0:
                                    levelE = [float(e) for e in levels]
                                    print('     Subentry',subent,' produces in',reaction,to,'the levels',repr(levelE))
                                    levelsSeen.append(['Subentry %s produces in %s%s levels %s' % (subent,reaction,to,repr(levelE)) ])
                                                                

                    searched +=  searchX4.replace('***** Search X4 for','') + " : %s" % searchesFound

    print('partitions:',partitions)
    print('Searches found:',searched)                

    if len(subentries)>0:
        fn = dir + '/datafile.props.csv'
        print('\nWrite csv file',fn,'for',len(subentries),'subentries')
        f = open ( fn , 'w')
        print('projectile,ejectile,target,residual,level,file,sys-error,stat-error,angle-integrated,norm,group,splitnorms,lab,abserr,scale,filedir,Aflip,Ein,eshift,ecalib,splitshifts,ratioRuth,S_factor,Npoints,EXFOR,minRelErr,maxRelErr',file=f)
        ents = set()
        for subent in sorted(list(subentries)):
            psubent = subent.replace('"','')
            npts,x4_tag = psubent.split(',')[-2:]
            file = psubent.split(',')[3]
            print(psubent)
            if npts == '0':
#                 print('   Not in csv as Npoints =',npts)
                excuses[x4_tag] = "No points. Possibly unprocessed photonuclear"
#                 os.remove(file)
            else:
                print(subent.replace(", ",","), file=f)
        print('\nPartitions:')
    
    if len(list(excuses.keys()))>0: print("\nReasons for exclusions:")
    for sub in list(excuses.keys()):
        if 'Highest energy' in excuses[sub]: print('   ',sub,':',excuses[sub])
    rels = 0
    for sub in list(excuses.keys()):
        if 'relative' in excuses[sub]: 
            print('   ',sub,':',excuses[sub])
            rels += 1
    for sub in list(excuses.keys()):
        if 'Highest energy' not in excuses[sub] and 'relative' not in excuses[sub]:  
            print('   ',sub,':',excuses[sub])
