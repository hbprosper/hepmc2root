#!/usr/bin/env python
# -----------------------------------------------------------------------
# File: hepmc2root.py
# Description: write events in HepMC2 format to a flat ROOT ntuple using
#              variable length arrays
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
#          15-Apr-2019 HBP test that ROOT can be imported
#          31-Jan-2020 HBP make compatible with Python 3
# -----------------------------------------------------------------------
import os, sys
try:
    import ROOT
except:
    sys.exit('\n\t*** This program uses ROOT! ***\n')
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

class hepmc2root:
    
    def __init__(self, filename, outfilename=None, treename=TREENAME, complevel=2):

        # check that file exists
        
        if not os.path.exists(filename):
            sys.exit("** hepmc2root.py: can't open file %s" % \
                         filename)
        self.inp = open(filename)
        inp = self.inp

        # get version number of HepMC
        
        self.header = [] # cache HepMC header
        version = None
        for line in inp:
            self.header.append(line)
            version = str.strip(line)
            if version == '': continue
            token = str.split(version)
            if token[0] == 'HepMC::Version':
                version = token[1]
            break
        else:
            sys.exit("** hepmc2root.py: format problem in file %s" % \
                         filename)
            
        print("HepMC version: %s" % version)

        # skip start of listing
        
        for line in inp:
            self.header.append(line)
            break

        # open output root file

        if outfilename == None:
            outfilename = '%s.root' % nameonly(filename)
            
        self.file = ROOT.TFile(outfilename, "recreate")
        self.tree = ROOT.TTree(treename, 'created: %s HepMC %s' % \
                                   (ctime(), version))

        # define event struct
        
        self.struct = '''struct Bag {
        int    Event_number;
        int    Event_numberMP;
        double Event_scale;
        double Event_alphaQCD;
        double Event_alphaQED;
        int    Event_barcodeSPV;
        int    Event_numberV;
        int    Event_barcodeBP1;
        int    Event_barcodeBP2;
        int    Event_numberP;

        double Xsection_value;
        double Xsection_error;

        int    PDF_parton1;
        int    PDF_parton2;
        double PDF_x1;
        double PDF_x2;
        double PDF_Q2;
        double PDF_x1f;
        double PDF_x2f;
        int    PDF_id1;
        int    PDF_id2;

        double Particle_x[%(size)d];
        double Particle_y[%(size)d];
        double Particle_z[%(size)d];
        double Particle_ctau[%(size)d];

        double Particle_barcode[%(size)d];
        int    Particle_pid[%(size)d];
        double Particle_px[%(size)d];
        double Particle_py[%(size)d];
        double Particle_pz[%(size)d];
        double Particle_energy[%(size)d];
        double Particle_mass[%(size)d];
        int    Particle_status[%(size)d];
        int    Particle_d1[%(size)d];
        int    Particle_d2[%(size)d];
};''' % {'size': MAXPART}

        # indices to vertices
        
        self.pvertex = [0]*MAXPART
        
        # create struct
        
        ROOT.gROOT.ProcessLine(self.struct)
        from ROOT import Bag
        self.bag = Bag()

        # create branches
        
        self.branch = []
        recs = str.split(self.struct, '\n')[1:-1]
        for rec in recs:
            t = str.split(rec)
            if len(t) == 0: continue
                
            fmt, name = t
            T = str.upper(fmt[0])
            name = name[:-1] # skip ";"
            # check for variable length array
            if name[-1] == ']':
                field = str.split(name, '[')[0]
                fmt   = '%s[Event_numberP]/%s' % (field, T)
            else:
                field = name
                fmt   = '%s/%s' % (field, T)
            self.branch.append(self.tree.Branch(field,
                                                ROOT.addressof(self.bag, field),
                                                    fmt))
        # list branches
        
        for ii, b in enumerate(self.branch):
            bname = b.GetName()
            leaves= b.GetListOfLeaves()
            if leaves == None:
                sys.exit("** hepmc2root: no list of leaves found for branch %s" % bname)
            leaf = leaves[0]
            if leaf == None:
                sys.exit("** hepmc2root: no leaf found for branch %s" % bname)
            leafname = leaf.GetName()
            leaftype = leaf.GetTypeName()
            print("%4d\t%-20s\t%s" % (ii+1, bname, leaftype))
            
    def __del__(self):
        self.tree.Write("", ROOT.TObject.kOverwrite)
        
    def __str__(self, index):
        bag = self.bag
        d   = " <%4d, %4d>" % (bag.Particle_d1[index], bag.Particle_d2[index])
        px  = bag.Particle_px[index]
        py  = bag.Particle_py[index]
        pt  = sqrt(px**2+py**2)
        rec = '%-14s %7d %4d %3d %7.1f (%7.1f, %7.1f, %7.1f, %7.1f)%s' \
          % (particleName(bag.Particle_pid[index]),
             bag.Particle_pid[index],
             bag.Particle_barcode[index],
             bag.Particle_status[index],
             pt, 
             bag.Particle_energy[index],
             bag.Particle_px[index],
             bag.Particle_py[index],
             bag.Particle_pz[index],
             d)
        return rec

    def __call__(self):
        inp = self.inp
        bag = self.bag

        self.event = [] # cache HepMC event in original format
        
        # find start of event
        
        token = None
        for line in inp:
            self.event.append(line)
            token = str.split(line)
            key   = token[0]
            if key != 'E': continue
            if debug > 0:
                print('BEGIN event')
            break
        else:
            return False

        if token == None:
            sys.exit("** hepmc2root.py: can't find start of event")

        bag.Event_number     = int(token[1])
        bag.Event_numberMP   = int(token[2])  # number of multi-particle interactions
        bag.Event_scale      = float(token[3])
        bag.Event_alphaQCD   = float(token[4])
        bag.Event_alphaQED   = float(token[5])
        bag.Event_processID  = int(token[6])
        bag.Event_barcodeSPV = int(token[7])
        bag.Event_numberV    = int(token[8])  # number of vertices in event
        bag.Event_barcodeBP1 = int(token[9])  # barcode beam particle 1 
        bag.Event_barcodeBP2 = int(token[10]) # barcode beam particle 2
        bag.Event_numberP    = 0               # number of particles
        
        if debug > 0:
            print("\tbarcode 1: %d" % self.barcode1)
            print("\tbarcode 2: %d" % self.barcode2)

        self.vertex = {}

        for line in inp:
            self.event.append(line)
            token = str.split(line)
            key = token[0]

            if key == 'C':
                # CROSS SECTION
                bag.Xsection_value = float(token[1])
                bag.Xsection_error = float(token[2])
                if debug > 0:
                    print("\tcross section: %10.3e +\- %10.3e pb" % \
                      (bag.Xsection_value, bag.Xsection_error))
                      
            elif key == 'F':
                # PDF INFO
                bag.PDF_parton1  = int(token[1])
                bag.PDF_parton2  = int(token[2])
                bag.PDF_x1       = float(token[3])
                bag.PDF_x2       = float(token[4])
                bag.PDF_Q2       = float(token[5])
                bag.PDF_x1f      = float(token[6])
                bag.PDF_x2f      = float(token[7])
                bag.PDF_id1      = int(token[8])
                bag.PDF_id2      = int(token[9])

                if debug > 0:
                    print('\tfound PDF info')

            elif key == 'V':
                # VERTEX
                vbarcode = int(token[1])
                self.vertex[vbarcode] = [-1, -1]
                x    = float(token[3])
                y    = float(token[4])
                z    = float(token[5])
                ctau = float(token[6])
                nout = int(token[8])
                if debug > 0:
                    if debug > 1:
                        print("\t%s" % token)
                    print('\tvertex(barcode): %10d' % vbarcode)
                    print('\tvertex(count):   %10d' % nout)

                # particles pertaining to this vertex follow immediately
                # after the vertex
                for ii in range(nout):
                    for line in inp:
                        self.event.append(line)
                        token  = str.split(line)
                        if debug > 1:
                            print("\t%s" % token)
                        key    = token[0]
                        if key != 'P':
                            sys.exit("** hepmc2root: faulty event record\n" + line)

                        if bag.Event_numberP < MAXPART:
                            index = bag.Event_numberP
                            bag.Event_numberP += 1

                            bag.Particle_x[index]       = x
                            bag.Particle_y[index]       = y
                            bag.Particle_z[index]       = z
                            bag.Particle_ctau[index]    = ctau
                            
                            bag.Particle_barcode[index] = int(token[1])
                            bag.Particle_pid[index]     = int(token[2])
                            bag.Particle_px[index]      = float(token[3])
                            bag.Particle_py[index]      = float(token[4])
                            bag.Particle_pz[index]      = float(token[5])
                            bag.Particle_energy[index]  = float(token[6])
                            bag.Particle_mass[index]    = float(token[7])
                            bag.Particle_status[index]  = int(token[8])
                            self.pvertex[index]         = int(token[11])

                            if ii == 0:
                                self.vertex[vbarcode][0] = index
                            else:
                                self.vertex[vbarcode][1] = index
                        
                        break
                    else:
                        return False
                    
            if len(self.vertex) >= bag.Event_numberV:
                for index in range(bag.Event_numberP):
                    code = self.pvertex[index]
                    if code in self.vertex:
                        d = self.vertex[code]
                        bag.Particle_d1[index] = d[0]
                        bag.Particle_d2[index] = d[1]
                    else:
                        bag.Particle_d1[index] = -1
                        bag.Particle_d2[index] = -1

                # fill ntuple
                
                self.file.cd()
                self.tree.Fill()
                
                return True
        else:
            return False
        
    def printTable(self):
        for ii in xrange(self.bag.Event_numberP):
            print("%4d\t%s" % (ii, self.__str__(ii)))
# -----------------------------------------------------------------------    
def main():
    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 1:
        sys.exit('''
    Usage:
        ./hepmc2root.py <HepMC-file> [output root file = <name>.root]
        ''')

    filename = argv[0]
    if argc > 1:
        outfilename = argv[1]
    else:
        outfilename = '%s.root' % nameonly(filename)

    stream = hepmc2root(filename, outfilename)

    ii = 0
    while stream():
        if ii % 1000 == 0:
            print(ii)
        ii += 1
# -----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print('\nciao!')
    
