#!/bin/env python

# Example script to validate LFV topo triggers
#
# Example input file:
# /afs/cern.ch/user/g/gerbaudo/public/tmp/for_marco/user.olya.11640704._000003.tmptrig.root
#
# davide.gerbaudo@gmail.com
# Jul 2017

import itertools
import optparse
import re
import os
import inspect
import string
import sys
from collections import defaultdict
from pprint import pprint
from math import pi, sqrt

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../L1TopoValidation/L1TopoCheck/python/")))
if cmd_subfolder not in sys.path:
     sys.path.insert(0, cmd_subfolder)


import ROOT as R
R.PyConfig.IgnoreCommandLineOptions = True # don't let root steal your cmd-line options
R.gROOT.SetBatch(1)                        # go batch!
R.gErrorIgnoreLevel = 9999                 # suppress messages about missing dict
                                           # (can't get rid of the 'duplicate' ones?)
R.gROOT.ProcessLine('#include "L1TopoCheck/TriggerBits.h"')
R.gROOT.ProcessLine('#include "L1TopoCheck/AlgorithmBits.h"')
R.gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')
EmTob = R.L1Topo.offline.EmTOB
JetTob = R.L1Topo.offline.JetTOB
# MuonTob = R.MuonTOB
# MuonTob = R.L1Topo.EnhancedMuonTOB
MuonTob = R.L1Topo.MuonTOB
TauTob = R.L1Topo.offline.TauTOB

import utils

def main():
    difference_dphi = 0   #CHANGE
    totat_dphi = 0
    usage = ("Usage : %prog [options] filename"
             "\n Examples :"
             "\n %prog  -v tmptrig.root"
             )

    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-n', '--num-events', default=None, type=int, help='number of events to process (default all)')
    parser.add_option('-s', '--skip-events', default=None, type=int, help='number of events to skip (default none)')
    parser.add_option('-v', '--verbose', default=False, action='store_true')
    parser.add_option('-d', '--debug', default=False, action='store_true')
    parser.add_option('-t', '--treename', default='trig')

    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    verbose = options.verbose
    debug = options.debug
    if verbose:
        utils.print_running_conditions(parser, options)

    input_filenames = utils.read_filename_arguments(args[0], options)
    if verbose:
        print 'Input files:'
        print '\n'.join(input_filenames)
    chain = R.TChain(options.treename)
    for input_filename in input_filenames:
        chain.Add(input_filename)  #chain beomes an array with the various files .root introduced
    num_available = chain.GetEntries()
    num_skip = options.skip_events
    num_toprocess = number_of_entries_to_process(num_available, options)
    if verbose:
        print "About to process %s (out of %d) entries: " % (num_toprocess, num_available)
    # # test: print all branches available in the tree
    # print 'chain ',chain
    # print 'branches:'
    # print 'list of branches ',chain.GetListOfBranches()
    # print '\n'.join([k.GetName() for k in chain.GetListOfBranches()])
    # print
    # return

    iEntry = 0
    possible_outcomes = ['pass_em_pass_hw', 'pass_em_fail_hw', 'fail_em_pass_hw', 'fail_em_fail_hw']
    algo_counters = defaultdict(int)
    item_counters = defaultdict(int)
    valid_counters = {k:0 for k in ['overflow'] + possible_outcomes}

    histos = {}
    histos2= {}
    lvl1item_name  = 'L1_LFV-MU6'
    algorithm_name = '0DR15-2MU6ab'

    l = lvl1item_name
    dr_binning    = (31 , -0.05 , 3.05)
    histos = {
            'dr_min_all'      : R.TH1F('dr_min_all'      , l+'; min #DeltaR 2mu6 pair'              , *dr_binning),
            'dr_min_hdw_pass' : R.TH1F('dr_min_hdw_pass' , l+'; min #DeltaR 2mu6 pair and hdw pass' , *dr_binning),
            'ratio'           : R.TH1F('ratio'           , l+'; efficienciy 0DR15-2MU6ab algorithm' , *dr_binning),
            }

    histo_names = [name for name, histo in histos.items()]

    for iEvent, event in enumerate(chain):
        if num_skip and iEvent<num_skip: continue
        if iEntry > num_toprocess: break
        # see how branches are created in TopoNtuple.py
        item_bits = item_tbits2bits(getattr(event,
                                            lvl1item_name.replace('-','_')+'_iBits'))
        increment_counters(item_counters, item_bits)
        # # These are not filled, so don't bother for now # DG-2016-06-23
        algo_bits = algo_tbits2bits(getattr(event, algorithm_name.replace('-','_')+'_0_aBits'))
        increment_counters(algo_counters, algo_bits)
        pass_hw = item_bits['TDT_TBP']
        pass_sim = item_bits['L1SIMULATED']

        overflown = algo_bits['OVERFLOWN']
        if overflown:
            valid_counters['overflow'] += 1
            continue


        muons = [Muon(tob.pt, tob.eta, tob.phi) for tob in event.emuMuonTOB
                 if tob.bcn==0] # only pick the ones from bunch crossing number 0
#        muons = [Muon(tob.Pt(), tob.Eta(), tob.Phi()) for tob in event.recomuon]
#                 if tob.Bcn()==0] # only pick the ones from bunch crossing number 0

        list_mu6ab = algo_MU6ab(muons) #mu6 list
        if len(list_mu6ab)<2: continue

        dr_min = algo_0DR15(list_mu6ab) #get minimum dr

        if dr_min >3: continue

        fill_histos(histos, dr_min, pass_hw) #fill histograms

    histos['ratio'].Divide(histos['dr_min_all'])

    c = R.TCanvas('c')    
    for name in histo_names:
        c.Clear()
        h = histos[name]
        h.Draw('h text')
        c.Update()
        c.SaveAs(name+'.png')
        c.SaveAs(name+'.root')



def algo_0DR15(muons): #retuns dr_min of all possible pairs
    dr_list = []
    n_mu = len(muons)
    for i in range(n_mu-1): #check all muon couples to see if Delta r is lower than 1.5
        for j in range(i+1,n_mu):
            dr = muons[i].p4.DeltaR(muons[j].p4)
            dr_list.append(dr)

    dr_list.sort() #sort list

    return dr_list[0]


def algo_MU6ab(muons): #returns list with all muons satifying MU6 sorted by energy
    mu6ab_list = []

    for muon in muons:
        pt = muon.p4.Pt()
        #some error when doing >=6000 instead of >5999 I think it is due to rounding errors changed to >5000
        #no difference noticed between >5999 and >5000. condition >5000 is the same as in official simulation.
        if pt>5000:
            muon          
            mu6ab_list.append(muon)
    
    mu6ab_list.sort(key = lambda mu6ab: mu6ab.p4.Pt()) #sort list
    return mu6ab_list 

    
def fill_histos(histos, dr, pass_hdw): #fills histograms
    
    histos['dr_min_all'].Fill(dr)
    if pass_hdw:
        histos['ratio'].Fill(dr)
        histos['dr_min_hdw_pass'].Fill(dr)

class Muon(object):
    def __init__(self, pt, eta, phi):
        muon_mass = 105.65
        self.p4 = R.TLorentzVector() # four-momentum
        self.p4.SetPtEtaPhiM(pt, eta, phi, muon_mass)

def algo_bit_names_and_numbers():
    "Bits stored in the TBits for each L1 algorithm"
    return (('FIRED'        , R.L1Topo.FIRED       ),
            ('OVERFLOWN'    , R.L1Topo.OVERFLOWN   ),
            ('TOB_EMULATED' , R.L1Topo.TOB_EMULATED),
            ('ROI_EMULATED' , R.L1Topo.ROI_EMULATED),
            ('CTP_TIP'      , R.L1Topo.CTP_TIP     ),
            ('TOPOSIM'      , R.L1Topo.TOPOSIM     ))

def item_bit_names_and_numbers():
    "Bits stored in the TBits for each L1 item"
    return (('TDT_TBP'    , R.L1Topo.TDT_TBP    ),
            ('TDT_TAP'    , R.L1Topo.TDT_TAP    ),
            ('TDT_TAV'    , R.L1Topo.TDT_TAV    ),
            ('L1EMULATED' , R.L1Topo.L1EMULATED ),
            ('L1SIMULATED', R.L1Topo.L1SIMULATED),
            ('CTP_TBP'    , R.L1Topo.CTP_TBP    ))

def algo_tbits2bits(bits=None):
    "convert TBits to dict"
    return dict((k, bits.TestBitNumber(b)) for k, b in algo_bit_names_and_numbers())

def item_tbits2bits(bits=None):
    "convert TBits to dict"
    return dict((k, bits.TestBitNumber(b)) for k, b in item_bit_names_and_numbers())

def increment_counters(counters={}, bits={}):
    for k, b in bits.items():
        counters[k] += 1 if b else 0


def number_of_entries_to_process(available_entries, options=None):
    N = available_entries
    n = options.num_events
    s = options.skip_events
    to_process = (min([N, n, N-s]) if n and s else
                  min([N, n]) if n else
                  N-s if s else
                  N)
    to_process = to_process if to_process > 0 else 0
    return to_process

if __name__=='__main__':
    main()
