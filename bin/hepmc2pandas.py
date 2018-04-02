#!/usr/bin/env python
# -----------------------------------------------------------------------
# File: hepmc2pandas.py
# Description: write events in HepMC2 format to a pandas file 
#
#
#     status = +- (10 * i + j)
#     + : still remaining particles
#     - : decayed/branched/fragmented/... and not remaining
#     i = 1 - 9 : stage of event generation inside PYTHIA
#     i = 10 -19 : reserved for future expansion
#     i >= 20 : free for add-on programs
#     j = 1 - 9 : further specification
#
# In detail, the list of used or foreseen status codes is:
#
#     11 - 19 : beam particles
#         11 : the event as a whole
#         12 : incoming beam
#         13 : incoming beam-inside-beam (e.g. gamma inside e)
#         14 : outgoing elastically scattered
#         15 : outgoing diffractively scattered
#     21 - 29 : particles of the hardest subprocess
#         21 : incoming
#         22 : intermediate (intended to have preserved mass)
#         23 : outgoing
#         24 : outgoing, nonperturbatively kicked out in diffraction
#
# Created: fall   2017 Harrison B. Prosper
# Updated: 04-Dec-2017 HBP add creation vertex (x,y,z) of particles.
# -----------------------------------------------------------------------
import os, sys
import pandas as pd
from string import split, strip, atoi, atof, upper
from math import sqrt
from time import ctime
from pnames import particleName
# -----------------------------------------------------------------------
def nameonly(s):
    import posixpath
    return posixpath.splitext(posixpath.split(s)[1])[0]

TREENAME= "Events"
MAXPART = 5000
debug = 0

class hepmcstream:
    
    def __init__(self, filename):

        # check that file exists
        
        if not os.path.exists(filename):
            sys.exit("** hepmcstream: can't open file %s" % filename)
        self.inp = open(filename)
        inp = self.inp

        # get version number of HepMC
        
        self.header = [] # cache HepMC header
        version = None
        for line in inp:
            self.header.append(line)
            version = strip(line)
            if version == '': continue
            token = split(version)
            if token[0] == 'HepMC::Version':
                version = token[1]
            break
        else:
            sys.exit("** hepmcstream: format problem in file %s" % filename)
            
        print "HepMC version: %s" % version

        # skip start of listing
        
        for line in inp:
            self.header.append(line)
            break

        # define event structure
        
        erec = '''{
        'Event_number': [],
        'Event_numberMP': [],
        'Event_scale': [],
        'Event_alphaQCD': [],
        'Event_alphaQED': [],
        'Event_processID': [],
        'Event_barcodeSPV': [],
        'Event_numberV': [],
        'Event_barcodeBP1': [],
        'Event_barcodeBP2': [],
        'Event_numberP': [],
        'Xsection_value': [],
        'Xsection_error': [],
        'PDF_parton1': [],
        'PDF_parton2': [],
        'PDF_x1': [],
        'PDF_x2': [],
        'PDF_Q2': [],
        'PDF_x1f': [],
        'PDF_x2f': [],
        'PDF_id1': [],
        'PDF_id2': []
    }'''       
        prec = '''{
        'Particle_name': [],
        'Particle_x': [],
        'Particle_y': [],
        'Particle_z': [],
        'Particle_ctau': [],
        'Particle_barcode': [],
        'Particle_pid': [],
        'Particle_px': [],
        'Particle_py': [],
        'Particle_pz': [],
        'Particle_energy': [],
        'Particle_mass': [],
        'Particle_status': [],
        'Particle_d1': [],
        'Particle_d2': []
    }'''
        self.event = eval(erec)
        self.particle = prec
        self.plist = []

        # indices to vertices        
        self.pvertex = [0]*MAXPART
            
    def __del__(self):
        pass
        
    def __str__(self, index):
        p = self.pbag
        d   = " <%4d, %4d>" % (p['Particle_d1'][index],
                                   p['Particle_d2'][index])
        px  = p['Particle_px'][index]
        py  = p['Particle_py'][index]
        pt  = sqrt(px**2+py**2)
        rec = '%-14s %7d %4d %3d %7.1f (%7.1f, %7.1f, %7.1f, %7.1f)%s' \
          % (particleName(p['Particle_pid'][index]),
                p['Particle_pid'][index],
                p['Particle_barcode'][index],
                p['Particle_status'][index],
                pt, 
                p['Particle_energy'][index],
                p['Particle_px'][index],
                p['Particle_py'][index],
                p['Particle_pz'][index],
                d)
        return rec
    
    def __call__(self):
        inp = self.inp

        # new event
        self.pbag = eval(self.particle)
        
        e = self.event
        p = self.pbag
        
        # find start of event        
        token = None
        for line in inp:
            token = split(line)
            key   = token[0]
            if key != 'E': continue
            if debug > 0:
                print 'BEGIN event'
            break
        else:
            return False

        if token == None:
            sys.exit("** hepmcstream: can't find start of event")

        e['Event_number'].append(atoi(token[1]))
        e['Event_numberMP'].append(atoi(token[2]))
        e['Event_scale'].append(atof(token[3]))
        e['Event_alphaQCD'].append(atof(token[4]))
        e['Event_alphaQED'].append(atof(token[5]))
        e['Event_processID'].append(atoi(token[6]))
        e['Event_barcodeSPV'].append(atoi(token[7]))
        e['Event_numberV'].append(atoi(token[8]))     # number of vertices in event
        e['Event_barcodeBP1'].append(atoi(token[9]))  # barcode beam particle 1 
        e['Event_barcodeBP2'].append(atoi(token[10])) # barcode beam particle 2
        e['Event_numberP'].append(0)                  # number of particles
        
        if debug > 2:
            print "\tbarcode 1: %d" % e['Event_barcodeBP1'][-1]
            print "\tbarcode 2: %d" % e['Event_barcodeBP2'][-1]

        self.vertex = {}

        for line in inp:
            token = split(line)
            key = token[0]

            if key == 'C':
                # CROSS SECTION
                e['Xsection_value'].append(atof(token[1]))
                e['Xsection_error'].append(atof(token[2]))
                if debug > 0:
                    print "\tcross section: %10.3e +\- %10.3e pb" % \
                      (e['Xsection_value'][-1], e['Xsection_error'][-1])
                      
            elif key == 'F':
                # PDF INFO
                e['PDF_parton1'].append(atoi(token[1]))
                e['PDF_parton2'].append(atoi(token[2]))
                e['PDF_x1'].append(atof(token[3]))
                e['PDF_x2'].append(atof(token[4]))
                e['PDF_Q2'].append(atof(token[5]))
                e['PDF_x1f'].append(atof(token[6]))
                e['PDF_x2f'].append(atof(token[7]))
                e['PDF_id1'].append(atoi(token[8]))
                e['PDF_id2'].append(atoi(token[9]))

                if debug > 0:
                    print '\tfound PDF info'

            elif key == 'V':
                # VERTEX
                vbarcode = atoi(token[1])
                self.vertex[vbarcode] = [-1, -1]
                x    = atof(token[3])
                y    = atof(token[4])
                z    = atof(token[5])
                ctau = atof(token[6])
                nout = atoi(token[8])
                if debug > 0:
                    if debug > 2:
                        print "\t%s" % token
                    print '\tvertex(barcode): %10d' % vbarcode
                    print '\tvertex(count):   %10d' % nout

                # particles pertaining to this vertex follow immediately
                # after the vertex
                for ii in xrange(nout):
                    for line in inp:
                        token  = split(line)
                        if debug > 0:
                            print "\t%s" % token
                        key    = token[0]
                        if key != 'P':
                            sys.exit("** hepmcstream: faulty event record\n" + line)

                        if e['Event_numberP'][-1] < MAXPART:
                            index = e['Event_numberP'][-1]
                            e['Event_numberP'][-1] += 1

                            p['Particle_x'].append(x)
                            p['Particle_y'].append(y)
                            p['Particle_z'].append(z)
                            p['Particle_ctau'].append(ctau)
                            
                            p['Particle_barcode'].append(atoi(token[1]))
                            p['Particle_pid'].append(atoi(token[2]))
                            
                            pid = p['Particle_pid'][-1]
                            p['Particle_name'].append(particleName(pid))
                            p['Particle_d1'].append(-1)
                            p['Particle_d2'].append(-1)
                            p['Particle_px'].append(atof(token[3]))
                            p['Particle_py'].append(atof(token[4]))
                            p['Particle_pz'].append(atof(token[5]))
                            p['Particle_energy'].append(atof(token[6]))
                            p['Particle_mass'].append(atof(token[7]))
                            p['Particle_status'].append(atoi(token[8]))
                            self.pvertex[index] = atoi(token[11])

                            if ii == 0:
                                self.vertex[vbarcode][0] = index
                            else:
                                self.vertex[vbarcode][1] = index
                        break
                    else:
                        return False
                    
            if len(self.vertex) >= e['Event_numberV'][-1]:
                for index in xrange(e['Event_numberP'][-1]):
                    code = self.pvertex[index]
                    if self.vertex.has_key(code):
                        d = self.vertex[code]

                        p['Particle_d1'][index] = d[0]
                        p['Particle_d2'][index] = d[1]
                    else:
                        p['Particle_d1'][index] = -1
                        p['Particle_d2'][index] = -1

                self.plist.append(pd.DataFrame(p))
                
                return True
        else:
            return False

    def save(self, filename):
        print "=> saving to file: %s" % filename
        self.event['Particle'] = self.plist
        df = pd.DataFrame(self.event)
        df.to_pickle(filename)
        print "=> done!"
        
    def printTable(self):
        for ii in xrange(self.event['Event_numberP'][-1]):
            print "%4d\t%s" % (ii, self.__str__(ii))
# -----------------------------------------------------------------------    
def main():
    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 1:
        sys.exit('''
    Usage:
        ./hepmc2pandas.py <HepMC-file> [output pickle file = <name>.pkl]
        ''')

    filename = argv[0]
    if argc > 1:
        outfilename = argv[1]
    else:
        outfilename = '%s.pkl' % nameonly(filename)

    stream = hepmcstream(filename)

    ii = 0
    while stream():
        if ii % 100 == 0:
            print ii
        ii += 1
    stream.save(outfilename)
# -----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print '\nciao!'
    
