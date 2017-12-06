#!/usr/bin/env python
# -----------------------------------------------------------------------
# File: hepmcfilter.py
# Description: read events from one HepMC file to another with a simple
#              filter between them.
#
#              Example 1: You have a HepMC file called susy200.hepmc and
#              you wish to select events in which the heavy neutral Higgs
#              boson decays to tau-/tau+. You can select these events
#              using the command
#
#              hepmcfilter.py susy200.hepmc 35 15 -15
#
#              Example 2: You now want to select events in which the heavy
#              neutral Higgs boson decays either to tau-/tau+ or top/antitop.
#              The command is
#
#              hepmcfilter.py susy200.hepmc 35 15 -15, 35 6 -6
# -----------------------------------------------------------------------
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
# Updated: 05-Dec-2017 HBP allow for selection of multiple decays
# -----------------------------------------------------------------------
import os, sys, ROOT
from string import split, strip, atoi, atof, upper, joinfields
# -----------------------------------------------------------------------
def nameonly(s):
    import posixpath
    return posixpath.splitext(posixpath.split(s)[1])[0]

TREENAME= "Events"
MAXPART = 5000
debug   = 0

class hepmcfilter:
    
    def __init__(self, filename, decays):

        # check that file exists
        
        if not os.path.exists(filename):
            sys.exit("** hepmcfilter: can't open file %s" % filename)
        self.inp = open(filename)
        inp = self.inp

        # get version number of HepMC
        
        self.header = [] # cache HepMC header
        for line in inp:
            self.header.append(line)
            line = strip(line)
            if line == '': continue
            token = split(line)
            if token[0] == 'HepMC::Version':
                break
        else:
            sys.exit("** hepmcfilter: format problem in file %s" % filename)

        # skip start of listing
        
        for line in inp:
            self.header.append(line)
            break

        self.eventsIn  = 0
        self.eventsOut = 0
        self.parents = set(decays.keys())
        self.decays  = decays
        self.outfilename = "filtered_%s.hepmc" %  nameonly(filename)
        self.out = open(self.outfilename, 'w')
        self.out.writelines(self.header)
        
    def __del__(self):
        self.out.close()
        
    def __call__(self):
        inp     = self.inp
        parents = self.parents
        decays  = self.decays
        
        self.event = [] # cache HepMC event in original format
        
        # find start of event
        
        found = False
        for line in inp:
            self.event.append(line)
            if line[0] != 'E': continue
            found = True
            self.eventsIn += 1
            break
        else:
            return False

        if not found:
            sys.exit("** hepmcfilter: can't find start of event")

        token = split(line)
        self.eventNumber = atoi(token[1])
        self.numberV     = atoi(token[8])  # number of vertices in event

        self.vertex  = {} # cache for event vertices
        self.parent  = {} # cache for final instance of parents in event

        # loop over vertices in event
        for line in inp:
            self.event.append(line)
            if line[0] == 'V':
                
                # found a vertex
                
                token = split(line)
                
                # get vertex barcode
                vbarcode = atoi(token[1]) 
                self.vertex[vbarcode] = None

                # get number of particles out of this vertex
                nout = atoi(token[8])
                if debug > 2:
                    print '\tvertex(barcode): %10d' % vbarcode
                    print '\tvertex(count):   %10d' % nout

                # particles pertaining to this vertex follow immediately
                # after the vertex
                d = [] # cache for PDG id of particles
                # loop over particles pertaining to current vertex
                for ii in xrange(nout):
                    for line in inp:
                        self.event.append(line)
                        if line[0] != 'P':
                            sys.exit("** hepmcfilter: faulty event record\n" + line)

                        token = split(line)
                        pid   = atoi(token[2])
                        d.append(pid)

                        # if this is not one of the desired parents
                        # continue to the next particle by breaking out
                        # of the inner loop
                        if not (pid in parents): break
                        
                        # found parent particle, cache its pid and vertex
                        self.parent[pid] = atoi(token[11])
                        if debug > 1:
                            print '\t\tpid(%d) vertex(%d)' % (pid, self.parent[pid])
                        # break out of inner loop and move to next particle
                        break
                    else:
                        return False # End of File
                    
                # cache PDG ids of daughters associated with current vertex
                self.vertex[vbarcode] = d
                
                if debug > 1:
                    print "%8d\t%s" % (vbarcode, self.vertex[vbarcode])

            # if we have not exhausted the number of vertices for this
            # event, continue to the next vertex
            if len(self.vertex) < self.numberV: continue

            # ----------------------------------------------------------                
            # we've found all vertices for current event
            # ----------------------------------------------------------
            # For each specified parent particle check if it matches a
            # parent ID in self.parent. If there is no match, we skip the
            # event because the desired particle is not present within the
            # event. If there is a match, we check if the decay of the
            # parent particle matches one of its required decays. An event
            # is selected if the decay of each required parent particle
            # matches at least one of its desired decays.

            keepEvent = True 

            # loop over required parents
            for pid in self.parents:

                # if required parent particle not present in the event
                # skip event.
                if not self.parent.has_key(pid):
                    keepEvent = False
                    if debug > 0:
                        print "**  parent with pid %d not present in event **" % pid
                    break

                # required parent is within event, so get "pointer"
                # to the vertex containing its daughters
                v  = self.parent[pid]  # get "pointer" to vertex
                if debug > 1:
                    print '\tPARENT(%d): vertex(%d)' % (pid, v)

                # if this parent has no daughters skip event
                if not self.vertex.has_key(v):
                    keepEvent = False
                    break

                # found vertex for current parent; get daughters
                daughters = set(self.vertex[v])
                if debug > 1:
                    print "\tVERTEX(%d): %s" % (v, daughters)

                # loop over required decays for current parent
                # and check if there is a match with the actual
                # decay.
                keep = False
                for d in decays[pid]:
                    c = 0
                    for p in d:
                        if p in daughters: c += 1
                    keep = keep or (c == len(d))
                    if keep: break

                keepEvent = keepEvent and keep

            if keepEvent:
                self.eventsOut += 1
                if debug > 0:
                    print "\tKEEP EVENT(%d)" % self.eventNumber
                self.out.writelines(self.event)

            # continue to next event
            return True
        else:
            self.out.write('HepMC::IO_GenEvent-END_EVENT_LISTING\n')
            return False # End of File
# -----------------------------------------------------------------------    
def main():
    argv = sys.argv[1:]
    argc = len(argv)
    if argc < 2:
        sys.exit('''
    Usage:
        hepmc2filter.py <HepMC-file> parent d1 .. [, ...]

    Example:
        hepmc2filter.py pMSSM_H0.hepmc 35 15 -15, 35 6 -6

        This will select events in which the heavy neutral Higgs boson decays
        either to tau-/tau+ or to top/~top.
        ''')

    filename = argv[0]

    # decode decays
    rec = joinfields(argv[1:], ' ')
    dd  = split(rec, ',')
    decays = {}
    for x in dd:
        t = map(atoi, split(x))
        parent = t[0]
        daughters = t[1:]
        if not decays.has_key(parent): decays[parent] = []
        decays[parent].append(daughters)
    
    # apply filter
    stream = hepmcfilter(filename, decays)
    ii = 0
    while stream():
        if ii % 100 == 0:
            print "%10d %10d" % (stream.eventsIn, stream.eventsOut)
        ii += 1

    print '''
Summary
    events(in):    %10d
    events(out):   %10d
    fraction:      %10.3e
''' % (stream.eventsIn, stream.eventsOut,
           float(stream.eventsOut)/stream.eventsIn)
# -----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print '\nciao!'
    
