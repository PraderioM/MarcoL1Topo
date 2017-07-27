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
    num_binning   = (9 , -0.5, 8.5)
    dr_binning    = (30 , 0 , 6)
    binning_0dr15 = (15 , 0.0 , 1.5)
    pt_binning    = (8, 3500.0 , 11500.0) 
    angle_binning = (28, -3.5, 3.5)
    for k in possible_outcomes: #initialize the histograms, they will still be empty after 
        histos[k] = {
            'n_mu'              : R.TH1F('n_mu'+'_'+k              , l+'; N input l1mus'                   , *num_binning),
            'n_mu6ab'           : R.TH1F('n_mu6ab'+'_'+k           , l+'; N mu6 muons'                     , *num_binning),
            'n_pairs_mu6_0dr15' : R.TH1F('n_pairs_mu6_0dr15'+'_'+k , l+'; N mu6_0dr15 pairs'               , *num_binning),
            'n_pairs_0dr15'     : R.TH1F('n_pairs_0dr15'+'_'+k     , l+'; N 0dr15 pairs'                   , *num_binning),
            'n_pairs_mu6ab'     : R.TH1F('n_pairs_mu6ab'+'_'+k     , l+'; N mu6 pairs'                     , *num_binning),
            'n_cand_pairs'      : R.TH1F('n_cand_pairs'+'_'+k      , l+'; N candidate pairs'               , *num_binning),
            'dr_min'            : R.TH1F('dr_min'+'_'+k            , l+'; min #DeltaR best candidate pair' , *dr_binning),
            'dr_0dr15'          : R.TH1F('dr_0dr15'+'_'+k          , l+'; #DeltaR 0dr15 pairs'             , *binning_0dr15),
            'dr_mu6_min'        : R.TH1F('dr_mu6_min'+'_'+k        , l+'; min #DeltaR best mu6 pair'       , *dr_binning),
            'dr_mu6'            : R.TH1F('dr_mu6'+'_'+k            , l+'; #DeltaR mu6 pairs'               , *dr_binning),
            'dr_mu6_0dr15'      : R.TH1F('dr_mu6_0dr15'+'_'+k      , l+'; #DeltaR mu6_0dr15'               , *binning_0dr15),
            'dr_any'            : R.TH1F('dr_any'+'_'+k            , l+'; #DeltaR any candidate pair'      , *dr_binning),
            'Phi_any'           : R.TH1F('Phi_any'+'_'+k           , l+'; Phi angle any muon'              , *angle_binning),
            'Phi_mu6'           : R.TH1F('Phi_mu6'+'_'+k           , l+'; Phi angle any mu6 muon'          , *angle_binning),
            'Eta_any'           : R.TH1F('Eta_any'+'_'+k           , l+'; Eta angle any muon'              , *angle_binning),
            'Eta_mu6'           : R.TH1F('Eta_mu6'+'_'+k           , l+'; Eta angle any mu6 muon'          , *angle_binning),
            'pt_0dr15'          : R.TH1F('pt_0dr15'+'_'+k          , l+'; #Pt 0dr15 muons'                 , *pt_binning),            
            'pt_any'            : R.TH1F('pt_any'+'_'+k            , l+'; #Pt any muon'                    , *pt_binning),
            }
        histos2[k] = R.TH2F('PhiEta_mu6'+'_'+k   , l+'; Phi angle any mu6; Eta angle any mu6' , *2*angle_binning)

    histo_names = [name for name, histo in histos[possible_outcomes[0]].items()]

    for iEvent, event in enumerate(chain):
        if num_skip and iEvent<num_skip: continue
        if iEntry > num_toprocess: break
        # see how branches are created in TopoNtuple.py
        item_bits = item_tbits2bits(getattr(event,
                                            lvl1item_name.replace('-','_')+'_iBits'))
        increment_counters(item_counters, item_bits)
        # # These are not filled, so don't bother for now # DG-2016-06-23
        #algo_bits = algo_tbits2bits(getattr(event, algorithm_name+'_0_aBits'))
        #increment_counters(algo_counters, algo_bits)
        pass_hw = item_bits['TDT_TBP']
        pass_sim = item_bits['L1SIMULATED']

        #overflown = algo_bits['OVERFLOWN']
        #if overflown:
        #    valid_counters['overflow'] += 1
        #    continue
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
#        muons = [Muon(tob.pt, tob.eta, tob.phi) for tob in event.emuMuonTOB
#                 if tob.bcn==0] # only pick the ones from bunch crossing number 0
        muons = [Muon(tob.pt, tob.eta, tob.phi) for tob in event.hdwMuonTOB
                 if tob.bcn==0] # only pick the ones from bunch crossing number 0
#        muons = remove_equal_muons(muons)

        list_mu6ab = algo_MU6ab(muons) #mu6 list

#        if pass_sim!=pass_hw:
#            continue
#        Phi3 = algo_PHI3(muons)
#        if not Phi3:
#            continue

        list_0dr15_pairs, list_pairs, list_0dr15= algo_0DR15(muons, muonList=True) #0dr15 couplelist
        #aux = [muon for muon, Pt in list_mu6ab]
#        list_mu6_0dr15_pairs, list_mu6_pairs = algo_0DR15([muon for muon in list_mu6ab], printout = pass_hw) #mu6_0dr15 couplelist
        list_mu6_0dr15_pairs, list_mu6_pairs = algo_0DR15(list_mu6ab) #mu6_0dr15 couplelist
        list_0dr15_mu6 = algo_MU6ab_pairs(list_0dr15_pairs)

        pass_emul = len(list_mu6_0dr15_pairs)   #returns true if 2mu6ab and 0dr15
        #pass_emul = pass_sim
        #pass_emul = len(list_0dr15_mu6) or len(list_mu6_0dr15_pairs)
        #pass_emul = len(list_0dr15_mu6)

        outcome = ('pass_em_pass_hw' if pass_hw and pass_emul else
                   'pass_em_fail_hw' if pass_emul else
                   'fail_em_pass_hw' if pass_hw else
                   'fail_em_fail_hw')
        valid_counters[outcome] += 1
        fill_histos(histos[outcome], histos2[outcome], muons, list_mu6ab, list_0dr15,
                    list_0dr15_pairs, list_pairs, list_mu6_0dr15_pairs, list_mu6_pairs) #fill histograms

#        if outcome == 'pass_em_fail_hw':
#        if outcome == 'pass_em_fail_hw':
#        if len(list_mu6_pairs):
        if True:
#            mu1 = list_mu6_pairs[0].muon1.p4
#            mu2 = list_mu6_pairs[0].muon2.p4
#            print ("runNumber = {0:d}  eventNumber = {1:d}".format(event.runNumber, event.eventNumber))
#            print("muon1:  Pt = {0:.2f}\t\tPhi = {1:.0f} \tEta = {2:.0f}".format(mu1.Pt()/1000, mu1.Phi()*10, mu1.Eta()*10))
#            print("muon2:  Pt = {0:.2f}\t\tPhi = {1:.0f} \tEta = {2:.0f}\n".format(mu2.Pt()/1000, mu2.Phi()*10, mu2.Eta()*10))
            for  muon in muons:
                mu = muon.p4
                print("Pt = {0:.0f}\t\tPhi = {1:.2f}  \tEta = {2:.2f}".format(mu.Pt()/1000., mu.Phi()*10, mu.Eta()*10))
            print("-------------")
        

#        if outcome == 'fail_em_pass_hw':
#        if True:
        if False:
#            print ("runNumber = {0:d}  eventNumber = {1:d}".format(event.runNumber, event.eventNumber))
#            print("emuMuonTOB")
#            for i, mu in enumerate(list(set(event.emuMuonTOB))):
#                print("<{:d}>: Pt = {:.2f} Phi = {:.2f} Eta = {:.2f}".format(i, mu.pt, mu.phi, mu.eta))

            for i, mu in enumerate(event.l1muon):
                print("Pt = {:.0f}\t\tPhi = {:.2f}  \tEta = {:.2f}  \tisVetoed = {:}".format(mu.p4.Pt()/1000., mu.p4.Phi()*10, mu.p4.Eta()*10, mu.isVetoed))

#            print("hdwMuonTOB")
#            for i, mu in enumerate(list(set(event.hdwMuonTOB))):
#                print("<{:d}>: Pt = {:.2f} Phi = {:.2f} Eta = {:.2f}".format(i, mu.pt, mu.phi, mu.eta))
            print("-------------")

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
        total_discordance = 100.*(f_p+p_f)/total_imputs
        pass_em_discordance = 100.*p_f/total_pass_em
        pass_hw_discordance = 100.*f_p/total_pass_hw
        print('  total   error {:.2f}%'.format(total_discordance))
        print('  em pass error {:.2f}%'.format(pass_em_discordance))
        print('  hw pass error {:.2f}%'.format(pass_hw_discordance))

#    print('Difference Dphi = {:.2f}%'.format(100.*difference_dphi/totat_dphi))

    c = R.TCanvas('c')
    order = [2,4,3,1] #is used to reorganize the histograms so that the firs line corresponds to pass emu and the first column to pass hdw
    
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
            if verbose:
                print('\n')
                h.Print("all")
                print("- Mean     = %.3f +- %.3f"%(h.GetMean(),h.GetMeanError()))
                print("- Std Dev  = %.3f +- %.3f"%(h.GetStdDev(), h.GetStdDevError()))
                print("- Skewness = %.3f"%(h.GetSkewness()))
                print("- Kurtosis = %.3f"%(h.GetKurtosis()))
        if verbose:
            print('\n\n')
        
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
        if verbose:
            h.Print("all")
    if verbose:
        print('\n')

    c.SaveAs('PhiEta_mu6.png')
    c.SaveAs('PhiEta_mu6.root')

def algo_MU6ab_pairs(list_0dr15_pairs):
    list_0dr15_mu6 = []
    for pair in list_0dr15_pairs:
        if pair.muon1.p4.Pt()>5000 and pair.muon2.p4.Pt()>5000:
            list_0dr15_mu6.append(pair)

    return list_0dr15_mu6

#def algo_0DR15(muons, muonList=False, printout = False): #retuns ordered list with any couple of muons satisfying 0DR15
def algo_0DR15(muons, muonList=False, pass_hw = False, pass_sim = False): #retuns ordered list with any couple of muons satisfying 0DR15
    couples_any   = []
    n_mu = len(muons)
    for i in range(n_mu-1): #check all muon couples to see if Delta r is lower than 1.5
        for j in range(i+1,n_mu):
            pair = MuonPair(muons[i],muons[j])
            couples_any.append(pair)

    couples_any.sort(key = lambda couple: couple.dr) #sort list
#    couples_0dr15 = [couple for couple in couples_any if (not couple.isPhi3 and couple.dr<1.505) or (couple.isPhi3 and (couple.dr<1.505 or (couple.dr>2.75 and couple.dr<3.2)))] #take only 0dr15, different algorithm for Phi>3
    couples_0dr15 = [couple for couple in couples_any if couple.dr<=1.5] #take only 0dr15
#    couples_0dr15 = [couple for couple in couples_any if couple.dr<=1.5 or couple.dr1<=1.5] #take only 0dr15
#    couples_0dr15 = [couple for couple in couples_any if ((couple.muon1.p4.Pt()<=7 and couple.muon2.p4.Pt()<=7) and couple.dr1<=1.5) or ((couple.muon1.p4.Pt()>7 or couple.muon2.p4.Pt()>7) and couple.dr<=1.5)] #take only 0dr15 and something else
#    couples_0dr15 = [couple for couple in couples_any if couple.dr2<=2.25] #take only 0dr15 and something else
#    if False:
#    if not len(couples_0dr15) and (printout and False):
    if pass_hw == pass_sim and (pass_sim != len(couples_0dr15) and False):#not muonList):
        for couple in couples_any:
           #print some data to check
            mu1 = couple.muon1.p4
            mu2 = couple.muon2.p4
            dr = couple.dr
            dr2= couple.dr2
            DPhi = mu1.DeltaPhi(mu2)
            print("muon1: Pt = {0:.2f} Phi = {1:.2f} Eta = {2:.2f}".format(mu1.Pt(), mu1.Phi(), mu1.Eta()))
            print("muon2: Pt = {0:.2f} Phi = {1:.2f} Eta = {2:.2f}".format(mu2.Pt(), mu2.Phi(), mu2.Eta()))
            print("Dr  = {0:.2f}   Dr2 = {1:.3f}   DPhi = {2:.2f}".format(dr, dr2, DPhi))
#            DPhiBad = (couple[0].p4.Phi()-couple[1].p4.Phi())%(2*pi)
#            print("DPhibad = {:.2f}\n\n".format(DPhiBad))
        print('pass_hw = {0:d}   pass_sim = {1:d}   pass_em = {2:d}'.format(pass_hw,pass_sim,len(couples_0dr15)))
        print("#events {:d}".format(len(couples_any)))
        print("--------------")
    
    if muonList: #return a list of the muons that belong to at least one couple of couples_0dr15
        list_0dr15 = []
        for couple in couples_any:
            if couple.muon1 in list_0dr15:
                pass
            else:
                list_0dr15.append(couple.muon1)
            
            if couple.muon2 in list_0dr15:
                pass
            else:
                list_0dr15.append(couple.muon2)
        return (couples_0dr15, couples_any, list_0dr15)

    return (couples_0dr15, couples_any)


def algo_PHI3(muons): #returns 1 if there is a muon satisfying PHI3 0 otherwise
    for muon in muons:
        if muon.p4.Phi()>3:
            return 1
    return 0



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

    
def fill_histos(histos, histos2, muons, list_mu6ab, list_0dr15, list_0dr15_pairs,
                list_pairs, list_mu6_0dr15_pairs, list_mu6_pairs): #fills histograms
    n_mu = len(muons)
    n_mu6= len(list_mu6ab)
    n_0dr15_pairs = len(list_0dr15_pairs)
    n_mu6_pairs = len(list_mu6_pairs)
    n_mu6_0dr15_pairs = len(list_mu6_0dr15_pairs)
    n_pairs = len(list_pairs)
    n_0dr15 = len (list_0dr15)
    histos['n_mu'             ].Fill(n_mu)
    histos['n_mu6ab'          ].Fill(n_mu6)
    histos['n_pairs_mu6_0dr15'].Fill(n_mu6_0dr15_pairs)  #number of mu6_0dr15 pairs
    # same number obtained here than using the formula n*(n-1)/2
    histos['n_pairs_0dr15'    ].Fill(n_0dr15_pairs)      #number of 0dr15 pairs
    histos['n_pairs_mu6ab'    ].Fill(n_mu6_pairs)        #number of mu6 pairs
    histos['n_cand_pairs'     ].Fill(n_pairs)            #number of candidate pairs
   
    if n_0dr15_pairs: #fill histograms of dr
        for pair in list_0dr15_pairs:
            histos['dr_0dr15'].Fill(pair.dr)
   
    if n_mu6_0dr15_pairs:
        for pair in list_mu6_0dr15_pairs:
            histos['dr_mu6_0dr15'].Fill(pair.dr)
    
    if n_mu6_pairs: 
        histos['dr_mu6_min'].Fill(list_mu6_pairs[0].dr)
        for pair in list_mu6_pairs:
            histos['dr_mu6'].Fill(pair.dr)

    if n_pairs:
        histos['dr_min'].Fill(list_pairs[0].dr) #takes the lower dr possible
        for pair in list_pairs:
            histos['dr_any'].Fill(pair.dr)

    for muon in muons: #fill histogram of momentums
        histos['pt_any'].Fill(muon.p4.Pt())
        histos['Phi_any'].Fill(muon.p4.Phi())
        histos['Eta_any'].Fill(muon.p4.Eta())
   
    for muon in list_0dr15: #fill histogram of momentums
        histos['pt_0dr15'].Fill(muon.p4.Pt())
  
    for muon in list_mu6ab: #fill histograms of angles
        Phi = muon.p4.Phi()
        Eta = muon.p4.Eta()
        histos['Phi_mu6'].Fill(Phi)
        histos['Eta_mu6'].Fill(Eta)
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
#        pi = 3.141592653
#        tau = 2*pi
        self.muon1  = muon1
        self.muon2  = muon2
#        Phi1 = muon1.p4.Phi()
#        Phi2 = muon2.p4.Phi()
#        if Phi1>3.09:
#            Phi1=0
#        Phi2 = muon2.p4.Phi()
#        if Phi2>3.09:
#            Phi2=0
#        DPhi = min((Phi1-Phi2)%tau, (Phi1-Phi2+tau/2)%tau) 
#        DPhi = abs(Phi1-Phi2)
#        if DPhi>pi:
#            DPhi=2*pi-DPhi
#        Eta1 = muon1.p4.Eta()
#        Eta2 = muon2.p4.Eta()
#        DEta = (Eta1-Eta2)
#        self.dr2=DEta**2+DPhi**2
#        self.dr = sqrt(self.dr2)
#        self.dr1 = muon1.p4.DeltaR(muon2.p4)
        self.dr = muon1.p4.DeltaR(muon2.p4)
        self.isPhi3 = muon1.p4.Phi()>3 or muon2.p4.Phi()>3 #true if the Phi value of any of the muons is greater than 3.

def remove_equal_muons(muons): #remove repeated muons of a list
    ZERO = 0.1
    i = 0
#    muons = [muon for muon in muonscopy]
    while i<(len(muons)-1):
        j = i+1
        while j<len(muons):
            if (abs(muons[i].p4.Pt()-muons[j].p4.Pt())/1000.+abs(muons[i].p4.Eta()-muons[j].p4.Eta())*10+abs(muons[i].p4.Phi()-muons[j].p4.Phi())*10)<ZERO:
                muons.pop(j)
                j-=1
            j+=1
        i+=1

    return muons

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
