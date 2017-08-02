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
from math import pi, sqrt, cos, cosh

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
    lvl1item_name  = 'L1_BPH-2M9-MU6MU4_BPH-0DR15-MU6MU4'
    algorithm1_name = '2INVM9-MU6ab-MU4ab'
    algorithm2_name = '0DR15-MU6ab-MU4ab'

    l = lvl1item_name
    num_binning   = (9 , -0.5, 8.5)
    dr_binning    = (30 , 0.0 , 6.0)
    m_binning     = (15 , -500.0 , 14500.0)
    pt_binning    = (8, 3500.0 , 11500.0) 
    angle_binning = (28, -3.5, 3.5)
    for k in possible_outcomes: #initialize the histograms, they will still be empty after 
        histos[k] = {
            'n_mu'                     : R.TH1F('n_mu'+'_'+k                     , l+'; N input l1mus'                , *num_binning),
            'n_mu4'                    : R.TH1F('n_mu4'+'_'+k                    , l+'; N mu4 muons'                  , *num_binning),
            'n_pairs_mu6mu4_2m9_0dr15' : R.TH1F('n_pairs_mu6mu4_2m9_0dr15'+'_'+k , l+'; N mu6mu4_2m9_0dr15 pairs'     , *num_binning),
            'n_pairs_mu6mu4_2m9'       : R.TH1F('n_pairs_mu6mu4_2m9'+'_'+k       , l+'; N mu6mu4_2m9 pairs'           , *num_binning),
            'n_pairs_mu6mu4_0dr15'     : R.TH1F('n_pairs_mu6mu4_0dr15'+'_'+k     , l+'; N mu6mu4_0dr15 pairs'         , *num_binning),
            'n_pairs_mu6mu4'           : R.TH1F('n_pairs_mu6mu4'+'_'+k           , l+'; N mu6mu4 pairs'               , *num_binning),
            'n_pairs_any'              : R.TH1F('n_pairs_any'+'_'+k              , l+'; N any pair'                   , *num_binning),
            'dr_any'                   : R.TH1F('dr_any'+'_'+k                   , l+'; #DeltaR any pair'             , *dr_binning),
            'dr_mu6mu4'                : R.TH1F('dr_mu6mu4'+'_'+k                , l+'; #DeltaR mu6mu4 pairs'         , *dr_binning),
            'dr_min_mu6mu4'            : R.TH1F('dr_min_mu6mu4'+'_'+k            , l+'; min #DeltaR mu6mu4 pairs'     , *dr_binning),
            'dr_mu6mu4_0dr15'          : R.TH1F('dr_mu6mu4_0dr15'+'_'+k          , l+'; #DeltaR mu6mu4_0dr15'         , *dr_binning),
            'dr_min_mu6mu4_0dr15'      : R.TH1F('dr_min_mu6mu4_0dr15'+'_'+k      , l+'; min #DeltaR mu6mu4_0dr15'     , *dr_binning),
            'dr_mu6mu4_2m9'            : R.TH1F('dr_mu6mu4_2m9'+'_'+k            , l+'; #DeltaR mu6mu4_2m9'           , *dr_binning),
            'dr_min_mu6mu4_2m9'        : R.TH1F('dr_min_mu6mu4_2m9'+'_'+k        , l+'; min #DeltaR mu6mu4_2m9'       , *dr_binning),
            'dr_mu6mu4_2m9_0dr15'      : R.TH1F('dr_mu6mu4_2m9_0dr15'+'_'+k      , l+'; #DeltaR mu6mu4_2m9_0dr15'     , *dr_binning),            'dr_mu6mu4'                : R.TH1F('dr_mu6mu4'+'_'+k                , l+'; #DeltaR mu6mu4 pairs'         , *dr_binning),
            'dr_min_mu6mu4_2m9_0dr15'  : R.TH1F('dr_min_mu6mu4_2m9_0dr15'+'_'+k  , l+'; min #DeltaR mu6mu4_2m9_0dr15' , *dr_binning),            'dr_mu6mu4'                : R.TH1F('dr_mu6mu4'+'_'+k                , l+'; #DeltaR mu6mu4 pairs'         , *dr_binning),
            'm_any'                    : R.TH1F('m_any'+'_'+k                    , l+'; #InvMass any pair'            , *m_binning),
            'm_mu6mu4'                 : R.TH1F('m_mu6mu4'+'_'+k                 , l+'; #InvMass mu6mu4 pairs'        , *m_binning),
            'm_mu6mu4_0dr15'           : R.TH1F('m_mu6mu4_0dr15'+'_'+k           , l+'; #InvMass mu6mu4_0dr15'        , *m_binning),
            'm_mu6mu4_2m9'             : R.TH1F('m_mu6mu4_2m9'+'_'+k             , l+'; #InvMass mu6mu4_2m9'          , *m_binning),
            'm_mu6mu4_2m9_0dr15'       : R.TH1F('m_mu6mu4_2m9_0dr15'+'_'+k       , l+'; #InvMass mu6mu4_2m9_0dr15'    , *m_binning),
            'Phi_mu4'                  : R.TH1F('Phi_mu6mu4'+'_'+k               , l+'; Phi angle any mu4 muon'       , *angle_binning),
            'Eta_mu4'                  : R.TH1F('Eta_mu6mu4'+'_'+k               , l+'; Eta angle any mu4 muon'       , *angle_binning),
            'pt_any'                   : R.TH1F('pt_any'+'_'+k                   , l+'; #Pt any muon'                 , *pt_binning),
            }
        histos2[k] = R.TH2F('PhiEta_mu6mu4'+'_'+k   , l+'; Phi angle any mu6mu4; Eta angle any mu6mu4' , *2*angle_binning)

    histo_names = [name for name, histo in histos[possible_outcomes[0]].items()]

    for iEvent, event in enumerate(chain):
        if num_skip and iEvent<num_skip: continue
        if iEntry > num_toprocess: break
        # see how branches are created in TopoNtuple.py
        item_bits = item_tbits2bits(getattr(event,
                                            lvl1item_name.replace('-','_')+'_iBits'))
        increment_counters(item_counters, item_bits)
        # # These are not filled, so don't bother for now # DG-2016-06-23
        algo1_bits = algo_tbits2bits(getattr(event, algorithm1_name.replace('-','_')+'_0_aBits'))
        increment_counters(algo_counters, algo1_bits)
        algo2_bits = algo_tbits2bits(getattr(event, algorithm2_name.replace('-','_')+'_0_aBits'))
        increment_counters(algo_counters, algo2_bits)
        pass_hw = item_bits['TDT_TBP']
        pass_sim = item_bits['L1SIMULATED']
        pass_hw = item_bits['TDT_TBP']
        
        overflown1 = algo1_bits['OVERFLOWN']
        overflown2 = algo2_bits['OVERFLOWN']
        if overflown1 and overflown2:
            valid_counters['overflow'] += 1
            continue
        # emTobs = [EmTob(w) for w in event.emTobs]
        # jetTobs = [JetTob(w) for w in event.jetTobs]
        # tauTobs = [TauTob(w) for w in event.tauTobs]
        # if debug:
        #     print 'emTobs[%d]' % len(emTobs)
        #     for i, et in enumerate(emTobs):
        #         print "[%d] (%f, %f)"%(i, et.eta(), et.phi())
        #     print 'jetTobs[%d]' % len(jetTobs)
        #     for i, jt in enumerate(jetTobs):
        #         print "[%d] (%f, %f)"%(i, jt.eta(), jt.phi())
        #     print 'tauTobs[%d]' % len(tauTobs)
        #     for i, tt in enumerate(tauTobs):
        #         print "[%d] (%f, %f)"%(i, tt.eta(), tt.phi())

        # these are EnhancedMuonTOB objects
        muons = [Muon(tob.pt, tob.eta, tob.phi) for tob in event.hdwMuonTOB
                 if tob.bcn==0] # only pick the ones from bunch crossing number 0
#        muons = remove_equal_muons(muons)
#        muons_emu = [Muon(tob.pt, tob.eta, tob.phi) for tob in event.emuMuonTOB
#                 if tob.bcn==0] # only pick the ones from bunch crossing number 0

#        muons = muons_emu

        list_mu4 = sorted(muons, key = lambda muon: muon.p4.Pt()) #all muons satisfy mu4
        list_mu6mu4_2m9_0dr15_pairs, list_mu6mu4_2m9_pairs, list_mu6mu4_0dr15_pairs, list_mu6mu4_pairs = algo_2M9_0DR15(list_mu4, pass_hw = pass_hw) #2im9_0dr15 couplelist
        list_pairs = make_all_pairs(muons)

        pass_emul = len(list_mu6mu4_2m9_0dr15_pairs)  #returns true if mu6mu4, 2m9 or 0dr15
#        pass_emul = len(list_mu6mu4_2m9_pairs) and len(list_mu6mu4_0dr15_pairs)  #returns true if mu6mu4, 2m9 or 0dr15

        outcome = ('pass_em_pass_hw' if pass_hw and pass_emul else
                   'pass_em_fail_hw' if pass_emul else
                   'fail_em_pass_hw' if pass_hw else
                   'fail_em_fail_hw')

#        if pass_hw and not len(list_mu6mu4_2m9_pairs):
#        if outcome == 'fail_em_pass_hw':
#        if not pass_hw:
        if False:
#            for muon in muons:
#                mu = muon.p4
#                print("Pt = {0:.0f}  \tPhi = {1:.0f}  \tEta = {2:.0f}".format(mu.Pt()/1000., mu.Phi()*10., mu.Eta()*10.))

#            print("")
#            for pair in list_mu6mu4_pairs:
#                mu1 = pair.muon1.p4
#                mu2 = pair.muon2.p4
#                print("muon1:  Pt = {:.0f}  \tPhi = {:.0f}  \tEta = {:.0f}".format(mu1.Pt()/1000., mu1.Phi()*10., mu1.Eta()*10.))
#                print("muon2:  Pt = {:.0f}  \tPhi = {:.0f}  \tEta = {:.0f}".format(mu2.Pt()/1000., mu2.Phi()*10., mu2.Eta()*10.))
#                print('invm2 = {:.0f}'.format(pair.invm2/1000000.))
            print('runNumber = {}  eventNumber = {}  lumiBlock = {}'.format(event.runNumber, event.eventNumber, event.lumiBlock))
#            print("-------------------------------------------")
#        else: continue

        if False:
            for  muon in muons:
                mu = muon.p4
                print("Pt = {0:.0f}\t\tPhi = {1:.0f}  \tEta = {2:.0f}".format(mu.Pt()/1000., mu.Phi()*10, mu.Eta()*10))
            print("")
            for  muon in muons_emu:
                mu = muon.p4
                print("Pt = {0:.0f}\t\tPhi = {1:.0f}  \tEta = {2:.0f}".format(mu.Pt()/1000., mu.Phi()*10, mu.Eta()*10))
            print("")
            print("")
        
#        if outcome == 'pass_em_fail_hw':
#        if outcome == 'fail_em_pass_hw':
        if False:
            print("all muons in event")
            for  muon in muons:
                mu = muon.p4
                print("Pt = {0:.0f}\t\tPhi = {1:.0f}  \tEta = {2:.0f}".format(mu.Pt()/1000., mu.Phi()*10, mu.Eta()*10))
#            print("pairs with a pass in emulation")
            print("pairs with a pass in 0DR15-MU6ab-MU4ab emulation")
            for pair in list_mu6mu4_0dr15_pairs:
                mu1 = pair.muon1.p4
                mu2 = pair.muon2.p4
                print("Pt = {0:.0f}\t\tPhi = {1:.0f}  \tEta = {2:.0f}".format(mu1.Pt()/1000., mu1.Phi()*10, mu1.Eta()*10))
                print("Pt = {0:.0f}\t\tPhi = {1:.0f}  \tEta = {2:.0f}".format(mu2.Pt()/1000., mu2.Phi()*10, mu2.Eta()*10))
                print("dr = {0:.2f}  \t\tinvm = {1:.2f}".format(pair.dr*10, pair.invm/1000.))
            print("pairs with a pass in 2INVM9-MU6ab-MU4ab emulation")
            for pair in list_mu6mu4_2m9_pairs:
                mu1 = pair.muon1.p4
                mu2 = pair.muon2.p4
                print("Pt = {0:.0f}\t\tPhi = {1:.0f}  \tEta = {2:.0f}".format(mu1.Pt()/1000., mu1.Phi()*10, mu1.Eta()*10))
                print("Pt = {0:.0f}\t\tPhi = {1:.0f}  \tEta = {2:.0f}".format(mu2.Pt()/1000., mu2.Phi()*10, mu2.Eta()*10))
                print("dr = {0:.2f}  \t\tinvm = {1:.2f}".format(pair.dr*10, pair.invm/1000.))
            print("-----------")


        valid_counters[outcome] += 1
        fill_histos(histos[outcome], histos2[outcome], muons, list_mu4, list_mu6mu4_2m9_pairs,
                     list_mu6mu4_0dr15_pairs, list_mu6mu4_2m9_0dr15_pairs,
                     list_mu6mu4_pairs, list_pairs) #fill histograms
        
        if debug and pass_hw:
            print "passed, %d muons" % len(muons)
        iEntry += 1


    print 'algo_counters:'
    pprint(dict(algo_counters))
    print 'item_counters:'
    pprint(dict(item_counters))
    print 'valid counters:'
    pprint(dict(valid_counters))

    if True:
        #print errors
        p_p=valid_counters['pass_em_pass_hw']
        p_f=valid_counters['pass_em_fail_hw']
        f_p=valid_counters['fail_em_pass_hw']
        f_f=valid_counters['fail_em_fail_hw']
    
        total_imputs = p_p+p_f+f_p+f_f
        total_pass_em = p_p+p_f
        total_pass_hw = f_p+p_p
        total_fail_em = f_p+f_f
        total_fail_hw = p_f+f_f
        total_discordance = 100.*(f_p+p_f)/total_imputs
        pass_em_discordance = 100.*p_f/total_pass_em
        fail_em_discordance = 100.*f_p/total_fail_em
        pass_hw_discordance = 100.*f_p/total_pass_hw
        fail_hw_discordance = 100.*p_f/total_fail_hw
        print('  total   error {:.2f}%'.format(total_discordance))
        print('  em pass error {:.2f}%'.format(pass_em_discordance))
        print('  em fail error {:.2f}%'.format(fail_em_discordance))
        print('  hw pass error {:.2f}%'.format(pass_hw_discordance))
        print('  hw fail error {:.2f}%'.format(fail_hw_discordance))

    c = R.TCanvas('c')
    order = [2,4,3,1]
    
    for name in histo_names:
        i = 0
        c.Clear()
        c.Divide(2,2)
        for outcome, hs in histos.items():
            h = histos[outcome][name]
            c.cd(order[i])
            h.Draw('h text')
            c.Update()
            i+=1
#        h = histos['fail_em_pass_hw'][name]
#        h.Draw('h text')
#        c.Update()
        
        c.SaveAs(name+'.png')
        c.SaveAs(name+'.root')
    
    i=0
    c.Clear()
    c.Divide(2,2)
    for outcome, h in histos2.items(): 
        c.cd(order[i])
        h.Draw('Colz')
        c.Update()
        i+=1
#    h = histos2['fail_em_pass_hw']
#    h.Draw('Colz')
#    c.Update()

#        if verbose:
#            h.Print("all")
#    if verbose:
#        print('\n')

    c.SaveAs('PhiEta_mu6mu4.png')
    c.SaveAs('PhiEta_mu6mu4.root')



def algo_2M9_0DR15(list_mu4, pass_hw): #retuns ordered list with couples of mu6mu4 muons satisfying 2M9_0DR15
    couples_any   = []
    list_mu6 = [muon for muon in list_mu4 if muon.p4.Pt()>5500]
    n_mu4 = len(list_mu4)
    n_mu6 = len(list_mu6)
    n_dif = n_mu4-n_mu6
    for i in range(n_dif): #build all mu4 mu6 pairs
        for muon in list_mu6:
            pair = MuonPair(list_mu4[i], muon)
            couples_any.append(pair)

    for i in range(n_mu6-1): #build all mu6 mu6 pairs
        for j in range(i+1,n_mu6):
            pair = MuonPair(list_mu6[i], list_mu6[j])
            couples_any.append(pair)

    couples_any.sort(key = lambda couple: couple.dr) #sort list
    couples_0dr15 = [couple for couple in couples_any if couple.dr<=1.5] #take only 0dr15
#    couples_0dr15 = [couple for couple in couples_any if (not couple.isPhi3 and couple.dr<1.505) or (couple.isPhi3 and (couple.dr<1.505 or (couple.dr>2.75 and couple.dr<3.2)))] #take only 0dr15, different algorithm for Phi>3
    couples_2m9   = [couple for couple in couples_any if couple.invm2>=4000000 and couple.invm<=81000000] #take only 2m9
#    couples_2m9   = [couple for couple in couples_any if couple.invm>=2000 and couple.invm<=9000] #take only 2m9
#    couples_2m9_0dr15 = [couple for couple in couples_2m9 if couple.dr<=1.5] #take only 2mu9_0dr15
    couples_2m9_0dr15 = [couple for couple in couples_0dr15 if couple.invm2>=4000000 and couple.invm2<=81000000] #take only 2mu9_0dr15
#    couples_2m9_0dr15 = [couple for couple in couples_0dr15 if couple.invm>=2000 and couple.invm<=9000] #take only 2mu9_0dr15
#    couples_2m9_0dr15 = couples_2m9
#    if pass_hw and not len(couples_2m9_0dr15):
    if False:
        for couple in couples_any:
           #print some data to check
            mu1 = couple.muon1.p4
            mu2 = couple.muon2.p4
            dr = couple.dr
            m = couple.invm
            DPhi = mu1.DeltaPhi(mu2)
            print("muon1: Pt = {0:.2f} Phi = {1:.2f} Eta = {2:.2f}".format(mu1.Pt(), mu1.Phi(), mu1.Eta()))
            print("muon2: Pt = {0:.2f} Phi = {1:.2f} Eta = {2:.2f}".format(mu2.Pt(), mu2.Phi(), mu2.Eta()))
            print("Dr  = {:.2f}   m = {:.2f}   DPhi = {:.2f}".format(dr, m, DPhi))
#            DPhiBad = (couple[0].p4.Phi()-couple[1].p4.Phi())%(2*pi)
#            print("DPhibad = {:.2f}\n\n".format(DPhiBad))
        print("#events {:d}".format(len(couples_any)))
        print("--------------")

    return (couples_2m9_0dr15, couples_2m9, couples_0dr15, couples_any)


def algo_MU4(muons): #returns sorted list of muons satistfying MU4
    mu4_list = []#

    for muon in muons:
        pt = muon.p4.Pt()
        #some error when doing >=4000 solved with >3500
        if pt>3500:                   
            mu4_list.append(muon)

    mu4_list.sort(key = lambda muon: muon.p4.Pt())
    
    return mu4_list

    
def fill_histos(histos, histos2, muons, list_mu4, list_mu6mu4_2m9_pairs,
                list_mu6mu4_0dr15_pairs, list_mu6mu4_2m9_0dr15_pairs,
                list_mu6mu4_pairs, list_pairs): #fills histograms
    n_mu = len(muons)
    n_mu4= len(list_mu4)
    n_pairs = len(list_pairs)
    n_mu6mu4_pairs = len(list_mu6mu4_pairs)
    n_mu6mu4_2m9_pairs = len(list_mu6mu4_2m9_pairs)
    n_mu6mu4_0dr15_pairs = len(list_mu6mu4_0dr15_pairs)
    n_mu6mu4_2m9_0dr15_pairs = len(list_mu6mu4_2m9_0dr15_pairs)
    histos['n_mu'                    ].Fill(n_mu)
    histos['n_mu4'                   ].Fill(n_mu4)
    histos['n_pairs_mu6mu4_2m9_0dr15'].Fill(n_mu6mu4_2m9_0dr15_pairs)  #number of mu6mu4_2m9_0dr15 pairs
    histos['n_pairs_mu6mu4_2m9'      ].Fill(n_mu6mu4_2m9_pairs)        #number of mu6mu4_2m9 pairs
    histos['n_pairs_mu6mu4_0dr15'    ].Fill(n_mu6mu4_0dr15_pairs)      #number of mu6mu4_0dr15 pairs
    histos['n_pairs_mu6mu4'          ].Fill(n_mu6mu4_pairs)            #number of mu6mu4 pairs
    histos['n_pairs_any'             ].Fill(n_pairs)                   #number of pairs
    
    
    if n_mu6mu4_2m9_0dr15_pairs:
        histos['dr_min_mu6mu4_2m9_0dr15'].Fill(list_mu6mu4_2m9_0dr15_pairs[0].dr)
        for pair in list_mu6mu4_2m9_0dr15_pairs:
            histos['dr_mu6mu4_2m9_0dr15'].Fill(pair.dr)
            histos['m_mu6mu4_2m9_0dr15'].Fill(pair.invm)
    
    if n_mu6mu4_2m9_pairs:
        histos['dr_min_mu6mu4_2m9'].Fill(list_mu6mu4_2m9_pairs[0].dr)
        for pair in list_mu6mu4_2m9_pairs:
            histos['dr_mu6mu4_2m9'].Fill(pair.dr)
            histos['m_mu6mu4_2m9'].Fill(pair.invm)
    
    if n_mu6mu4_0dr15_pairs:
        histos['dr_min_mu6mu4_0dr15'].Fill(list_mu6mu4_0dr15_pairs[0].dr)
        for pair in list_mu6mu4_0dr15_pairs:
            histos['dr_mu6mu4_0dr15'].Fill(pair.dr)
            histos['m_mu6mu4_0dr15'].Fill(pair.invm)

    if n_mu6mu4_pairs: 
        histos['dr_min_mu6mu4'].Fill(list_mu6mu4_pairs[0].dr) #lowest dr
        for pair in list_mu6mu4_pairs:
            histos['dr_mu6mu4'].Fill(pair.dr)
            histos['m_mu6mu4'].Fill(pair.invm)

    if n_pairs: 
        for pair in list_pairs:
            histos['dr_any'].Fill(pair.dr)
            histos['m_any'].Fill(pair.invm)
   
    for muon in muons: #fill histogram of momentums
        histos['pt_any'].Fill(muon.p4.Pt())
  
    for muon in list_mu4: #fill histograms of angles
        Phi = muon.p4.Phi()
        Eta = muon.p4.Eta()
        histos['Phi_mu4'].Fill(Phi)
        histos['Eta_mu4'].Fill(Eta)
        histos2.Fill(Phi, Eta)


class Muon(object):
    def __init__(self, pt, eta, phi):
        muon_mass = 105.65
        #tlv = R.TLorentzVector() # four-momentum
        #self.p4 = tlv.SetPtEtaPhiE(pt, eta, phi, energy)
        self.p4 = R.TLorentzVector() # four-momentum
        self.p4.SetPtEtaPhiM(pt, eta, phi, muon_mass)

class MuonPair(object):
    def __init__(self, muon1, muon2):
        pi = 3.141592653
        tau = 2*pi
        self.muon1  = muon1
        self.muon2  = muon2
        self.dr     = muon1.p4.DeltaR(muon2.p4)
#        DPhi = muon1.p4.DeltaPhi(muon2.p4)
#        Phi1 = muon1.p4.Phi()
#        if Phi1 < 0:
#            Phi1 -=.1
#        Phi2 = muon2.p4.Phi()
#        if Phi2 < 0:
#            Phi2 -=.1

#        if Phi1>=3.09:
#            Phi1=0
#        Phi2 = muon2.p4.Phi()
#        if Phi2>=3.09:
#            Phi2=0

#        if Phi1>3:
#            Phi1+=0.1
#        if Phi2>3:
#            Phi2+=0.1
#        DPhi = abs(Phi1-Phi2)
#        DPhi = min((Phi1-Phi2)%tau, (Phi1-Phi2+tau/2)%tau) 
#        if DPhi>pi:
#            DPhi=2*pi-DPhi
        DEta = abs(muon1.p4.Eta()-muon2.p4.Eta())
        DPhi = abs(muon1.p4.Phi()-muon2.p4.Phi())
        if DPhi > pi:
            DPhi = tau - DPhi
#        self.invm2 = round(2*round(muon1.p4.Pt()/1000.)*round(muon2.p4.Pt()/1000.)*(cosh(DEta)-cos(DPhi)))*1000000
        DPhi = muon1.p4.DeltaPhi(muon2.p4)
        Eta1 = muon1.p4.Eta()
        Eta2 = muon2.p4.Eta()
        DEta = Eta1-Eta2
        self.invm2 = round(2*muon1.p4.Pt()*muon2.p4.Pt()*(cosh(DEta)-cos(DPhi)))
        self.invm = sqrt(2*muon1.p4.Pt()*muon2.p4.Pt()*(cosh(DEta)-cos(DPhi)))
#        if muon1.p4.Pt()>=9999:
#            self.invm = sqrt(2*(muon1.p4.Pt()+1000)*muon2.p4.Pt()*(cosh(DEta)-cos(DPhi)))
#        elif muon2.p4.Pt()>=9999:
#            self.invm = sqrt(2*muon1.p4.Pt()*(muon2.p4.Pt()+1000)*(cosh(DEta)-cos(DPhi)))
#        if muon1.p4.Pt()>=9999 and muon2.p4.Pt()>=9999:
#            self.invm = sqrt(2*(muon1.p4.Pt()+1000)*(muon2.p4.Pt()+1000)*(cosh(DEta)-cos(DPhi)))

#        self.invm      = (muon1.p4+muon2.p4).M()
        self.isPhi3 = muon1.p4.Phi()>3 or muon2.p4.Phi()>3 #true if the Phi value of any of the muons is greater than 3.

def remove_equal_muons(muons): #remove repeated muons of a list
    ZERO = 1.e-2
    i = 0
    while i<(len(muons)-1):
        replaced = 0
        j = i+1
        while j<len(muons):
            if (abs(muons[i].p4.Pt()-muons[j].p4.Pt())+abs(muons[i].p4.Eta()-muons[j].p4.Eta())+abs(muons[i].p4.Phi()-muons[j].p4.Phi()))<ZERO:
                muons.pop(j)
                replaced+=1
                j-=1
            j+=1
        if replaced and abs(muons[i].p4.Pt()/1000.-10)<ZERO:
            muons[i] = Muon(2*muons[i].p4.Pt(), muons[i].p4.Eta(), muons[i].p4.Phi())
        i +=1

    return muons

def make_all_pairs(muons): #given a muon list returns all possible non repeated pairs of muons sorted by dr
    list_pairs = []
    n_mu = len(muons)
    for i in range(n_mu-1):
        for j in range(i+1, n_mu):
            pair = MuonPair(muons[i],muons[j])
            list_pairs.append(pair)

    list_pairs.sort(key = lambda pair: pair.dr)
    return list_pairs

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
