#! /usr/bin/env python
from __future__ import print_function

import sys
tr = {'barn64': 'Barnard',
      'Elwyn': 'Elwyn',
      'Elwdxs': 'Elwyn_dXS',
      'fas64': 'Fasoli',
      'fiedler': 'Fiedler',
      'har63': 'Harrison',
      'ivan67': 'Ivanovic',
      'mcc63': 'McCray',
      'spi66e': 'Spiger_3He3He',
      'spi66r': 'Spiger_3Hep',
      'Tumino': 'Tumino'
     }

projectiles = {
      'barn64': 3,
      'Elwyn':  1,
      'Elwdxs': 1,
      'fas64': 1,
      'fiedler': 1,
      'har63': 1,
      'ivan67': 4,
      'mcc63': 1,
      'spi66e': 3,
      'spi66r': 3,
      'Tumino': 1,
      'rxn12': 1,
     }

scale = 1e-3 # for mb in, barns out
frame = 'cm'
frame = 'lab'

m3He = 3.01550116	
m4He = 4.001506
e4to3 = m3He/m4He
lab2cm34 = m4He/(m4He + m3He)
lab2cm34 = m4He/(m4He + m3He)
mp = 1.008665
m6Li = 6.01347726
lab2cm16 = m6Li/(mp + m6Li)
if frame=='lab':
    conv = [1. for i in range(5)]
else: # cm
    conv = [0, lab2cm16, 0., lab2cm34, lab2cm43 ]

files = {}
error = 'abs'  # default
for line in open('data','r'):
    pp = line.split()
    if 'label=' in line: 
        label = line.split('label=')[1:]
    if len(pp)==7:
        ang,sig,err = [float(x) for x in pp[:3]]
        set = pp[3]
        en = float(pp[5])
        if set=='ivan67': en *= e4to3
        kind = pp[6]
        en *= conv[projectiles[set])

        error = 'rel'
        esig = sig*err/100.
        find = files.get(set,None)
        if find is None:
            out = tr[set] + '.dat'
            find = open(out,'w')
            files[set] = find
            print('File',out,'with',kind,' ',error,'errors')
        print(' %s %s %s %s' % (en,ang,sig*scale,esig*scale),file=find)

    if len(pp)==2:
        en = None
        if pp[1]=='kev': en = float(pp[0])*1e-3
        if pp[1]=='keV': en = float(pp[0])*1e-3
        if pp[1]=='Mev': en = float(pp[0])
        if pp[1]=='MeV': en = float(pp[0])
        # if en is not None: print("set E=",en,"MeV")
        error = 'abs'  # default
    if 'per cent error' in line: error = 'rel'   

    if len(pp)==4:
        try:
            ang,sig,err = [float(x) for x in pp[:3]]
            set = pp[3]
            esig = err if error=='abs' else err*sig/100.
            find = files.get(set,None)
            if '6Li(p' in label:
                conv = lab2cm16
                enc = en*conv
            else:
                enc = en
                if frame=='cm': print("Unknown frame set for ",set)
            if find is None:
                out = tr[set] + '.dat'
                find = open(out,'w')
                files[set] = find
                print('File',out,'with',kind,' ',error,'errors')
            print(' %s %s %s %s' % (enc,ang,sig*scale,esig*scale),file=find)
        except:
            pass

    if pp[0][:3]=='int':
        sig,err = [float(x) for x in pp[1:3]]
        set = pp[3]
        ang = 0
        unit = opp[1]
        escale = 1e-3 if unit=='keV' else 1.0
        en = float(opp[0]) * escale
        kind = pp[0]
        if frame=='cm': print("Unknown frame set",kind)
        error = 'abs'
        esig = err  # abs error
        find = files.get(set,None)
        if find is None:
            out = tr[set] + '.dat'
            find = open(out,'w')
            files[set] = find
            print('File',out,'with',kind,' E-unit:',unit,' ',error,'errors')
        print(' %s %s %s %s' % (en,ang,sig*scale,esig*scale),file=find)
    opp = pp
