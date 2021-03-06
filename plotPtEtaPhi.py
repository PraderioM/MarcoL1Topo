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
from array import array

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
             "\n %prog  -v tmptrig.root input.txt"
             )

    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-n', '--num-events', default=None, type=int, help='number of events to process (default all)')
    parser.add_option('-s', '--skip-events', default=None, type=int, help='number of events to skip (default none)')
    parser.add_option('-v', '--verbose', default=False, action='store_true')
    parser.add_option('-d', '--debug', default=False, action='store_true')
    parser.add_option('-t', '--treename', default='trig')

    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")

    Events_hdw = get_events_hdw_emu_offline_ef(args[0], options, 'hdw')
    Events_emu = get_events_hdw_emu_offline_ef(args[0], options, 'emu')
    Events_offline = get_events_hdw_emu_offline_ef(args[0], options, 'of')
    Events_ef = get_events_hdw_emu_offline_ef(args[0], options, 'ef')

    Events_sim = get_events_sim(args[1])

    Pt_sim, Pt_hdw_sim, Eta_sim, Eta_hdw_sim, Phi_sim, Phi_hdw_sim, matched_sim, total_sim = generate_pt_eta_phi(Events_sim, Events_hdw)
    Pt_emu, Pt_hdw_emu, Eta_emu, Eta_hdw_emu, Phi_emu, Phi_hdw_emu, matched_emu, total_emu = generate_pt_eta_phi(Events_emu, Events_hdw)
    Pt_emu_sim, Pt_sim_emu, Eta_emu_sim, Eta_sim_emu, Phi_emu_sim, Phi_sim_emu, matched_emu_sim, total_emu_sim = generate_pt_eta_phi(Events_emu, Events_sim)

    #Pt_offline, Pt_hdw_offline, Eta_offline, Eta_hdw_offline, Phi_offline, Phi_hdw_offline, matched_offline, total_offline = generate_pt_eta_phi(Events_offline, Events_hdw)
    #Pt_ef, Pt_hdw_ef, Eta_ef, Eta_hdw_ef, Phi_ef, Phi_hdw_ef, matched_ef, total_ef = generate_pt_eta_phi(Events_ef, Events_hdw)


    names = ['sim', 'emu']
    magnitudes = ['Pt', 'Eta', 'Phi']
    values = {  'sim'     : {'Pt' : {'sim': Pt_sim,      'hdw': Pt_hdw_sim},
                             'Eta': {'sim': Eta_sim,     'hdw': Eta_hdw_sim},
                             'Phi': {'sim': Phi_sim,     'hdw': Phi_hdw_sim}
                            },
                'emu'     : {'Pt' : {'emu': Pt_emu,      'hdw': Pt_hdw_emu},
                             'Eta': {'emu': Eta_emu,     'hdw': Eta_hdw_emu},
                             'Phi': {'emu': Phi_emu,     'hdw': Phi_hdw_emu}
                            },
                'emu_sim' : {'Pt' : {'emu': Pt_emu_sim,  'sim': Pt_sim_emu},
                             'Eta': {'emu': Eta_emu_sim, 'sim': Eta_sim_emu},
                             'Phi': {'emu': Phi_emu_sim, 'sim': Phi_sim_emu}
                            }
             }
    matches = {'sim': {'matched': matched_sim, 'total': total_sim},
               'emu': {'matched': matched_emu, 'total': total_emu},
               'emu_sim': {'matched': matched_emu_sim, 'total': total_emu_sim}
               }

    binnings = {'Pt': (22, -0.5 , 21.5), 'Eta': (28, -3.5, 3.5), 'Phi': (28, -3.5, 3.5)}
    histos = {}
    for name in names:
        histos[name] = {}
        for mag in magnitudes:
            histos[name][mag] = R.TH2F('hdw_'+name, mag+' hdw_'+name+'; '+mag+' hdw; '+mag+' '+name, *2*binnings[mag])
            for k in range(len(values[name][mag][name])):
                histos[name][mag].Fill(values[name][mag]['hdw'][k], values[name][mag][name][k])

    histos['emu_sim'] = {}
    for mag in magnitudes:
        histos['emu_sim'][mag] = R.TH2F('emu_sim'+name, mag+' emu_sim; '+mag+' emu; '+mag+' sim', *2*binnings[mag])
        for k in range(len(values['emu_sim'][mag]['emu'])):
            histos['emu_sim'][mag].Fill(values['emu_sim'][mag]['emu'][k], values['emu_sim'][mag]['sim'][k])


    c = R.TCanvas('c', '')
    draw_opt = ['colz3', 'box same']

    for name in names:
        c.Clear()
        c.Divide(2, 2, 0.01, 0.01)
        for i in range(len(magnitudes)):
            mag = magnitudes[i]
            c.cd(i+1)
            c.GetPad(i+1).SetGrid()
            for opt in draw_opt:
                histos[name][mag].Draw(opt)
        c.Update()
        c.SaveAs('hdw_'+name+'.png')
        c.SaveAs('hdw_'+name+'.root')
        n_match = matches[name]['matched']
        n_total = matches[name]['total']
        percent = 100.*n_match/n_total
        print('matches for hdw_'+name+' are {0:d} out of {1:d}\nMeaninig {2:.2f}%\n'.format(n_match, n_total, percent))

    c.Clear()
    c.Divide(2, 2, 0.01, 0.01)
    for i in range(len(magnitudes)):
        mag = magnitudes[i]
        c.cd(i+1)
        c.GetPad(i+1).SetGrid()
        for opt in draw_opt:
            histos['emu_sim'][mag].Draw(opt)
    c.Update()
    c.SaveAs('emu_sim.png')
    c.SaveAs('emu_sim.root')
    n_match = matches['emu_sim']['matched']
    n_total = matches['emu_sim']['total']
    percent = 100.*n_match/n_total
    print('matches for emu_sim are {0:d} out of {1:d}\nMeaninig {2:.2f}%\n'.format(n_match, n_total, percent))





class Muon(object):
    def __init__(self, pt, eta, phi):
        muon_mass = 105.65
        #tlv = R.TLorentzVector() # four-momentum
        #self.p4 = tlv.SetPtEtaPhiE(pt, eta, phi, energy)
        self.p4 = R.TLorentzVector() # four-momentum
        self.p4.SetPtEtaPhiM(pt, eta, phi, muon_mass)

class Event(object):
    def __init__(self, run_number, event_number, muons):
        self.run_number = run_number
        self.event_number = event_number
        self.muons = muons

class Candidate(object):
    def __init__(self, muon, position, difference):
        self.muon = muon
        self.position = position
        self. difference = difference

def remove_equal_muons(muons): #remove repeated muons of a list
    ZERO = 1.e-2
    i = 0
    while i<(len(muons)-1):
        j = i+1
        while j<len(muons):
            if (abs(muons[i].p4.Pt()-muons[j].p4.Pt())+abs(muons[i].p4.Eta()-muons[j].p4.Eta())+abs(muons[i].p4.Phi()-muons[j].p4.Phi()))<ZERO:
                muons.pop(j)
                j-=1
            j+=1
        i+=1

    return muons

def remove_repeated_events(event_list): #remove repeated events of a list
    i = 0
    while i<len(event_list)-1:
        j = i+1
        while j<len(event_list):
            if event_list[i].event_number == event_list[j].event_number:
                event_list.pop(j)
                j-=1
            j+=1
        i+=1

    return event_list

def get_events_hdw_emu_offline_ef(filenames, options, hdw_emu_of_ef):
    verbose = options.verbose
    debug = options.debug
    if verbose:
        utils.print_running_conditions(parser, options)

    input_filenames = utils.read_filename_arguments(filenames, options)
    if verbose:
        print 'Input files:'
        print '\n'.join(input_filenames)
    chain = R.TChain(options.treename)
    for input_filename in input_filenames:
        chain.Add(input_filename)  #chain beomes an array with the various files .root introduced
    num_available = chain.GetEntries()
    num_skip = options.skip_events
    num_toprocess = number_of_entries_to_process(num_available, options)

    iEntry = 0
    Events = []
    for iEvent, event in enumerate(chain):
        if num_skip and iEvent<num_skip: continue
        if iEntry > num_toprocess: break

        if hdw_emu_of_ef == 'hdw':
            muons = [Muon(tob.pt/1000., tob.eta, tob.phi) for tob in event.hdwMuonTOB
                    if tob.bcn==0] # only pick the ones from bunch crossing number 0
        elif hdw_emu_of_ef == 'emu':
            muons = [Muon(tob.pt/1000., tob.eta, tob.phi) for tob in event.emuMuonTOB
                    if tob.bcn==0] # only pick the ones from bunch crossing number 0
        elif hdw_emu_of_ef == 'of':
            muons = [Muon(tob.Pt()/1000., tob.Eta(), tob.Phi()) for tob in event.recomuon]
        else:
            muons = [Muon(tob.Pt()/1000., tob.Eta(), tob.Phi()) for tob in event.efmuon]

#        muons = remove_equal_muons(muons)
        muons.sort(key = lambda x: x.p4.Eta())

        entry = Event(event.runNumber, event.eventNumber, muons)
        Events.append(entry)
        iEntry +=1

    return Events



def get_events_sim(file_name):
    Events_sim = []
    sim_file = open(file_name, 'r')
    read_muons = 0
    read_event = 0
    run_number = 0
    event_number = 0

    for line in sim_file:
        if line == '</muon>\n':
            read_muons = 0
            continue

        if read_muons:
            line = line.split('  ')
            Pt = int(line[0])
            if Pt>10: Pt=10
            muon  = Muon(Pt, int(line[1])/10., int(line[2])/10.)
            Muons.append(muon)
            continue

        if line == '<muon>\n':
            Muons = []
            read_muons = 1
            continue
        
        if line == '</info>\n':
            read_event = 0
            Muons.sort(key = lambda x: x.p4.Eta())
            entry = Event(run_number, event_number, Muons)
            Events_sim.append(entry)
            continue

        if read_event:
            line = line.split('  ')
            run_number = int(line[0])
            event_number = int(line[1])
            continue

        if line == '<info>\n':
            read_event = 1
            continue

    Events_sim = remove_repeated_events(Events_sim)
    return Events_sim

def generate_pt_eta_phi(Events1, Events2, min_dPhi = 0): #this function takes 2 events lists and returns values of the found matches of 
    Pt1 = []
    Pt2 = []
    Eta1 = []
    Eta2 = []
    Phi1 = []
    Phi2 = []

    matched = 0
    total = 0
    for event1 in Events1:
        for event2 in Events2:
            if event1.event_number == event2.event_number and event1.run_number == event2.run_number:
                if len(event1.muons)<len(event2.muons):
                    for muon1 in event1.muons:
                        candidates = []
                        for imuon, muon2 in enumerate(event2.muons):
                            #dPt = (abs(muon1.p4.Pt()-muon2.p4.Pt()) if muon1.p4.Pt()<10 and muon2.p4.Pt()<10 else
                            dPt = abs(muon1.p4.Pt()-muon2.p4.Pt())
                            dPt /=10. #make the value be of the same order as the others.
                            dEta = abs(muon1.p4.Eta()-muon2.p4.Eta())
                            dPhi = (abs(muon1.p4.Phi()-muon2.p4.Phi()) if (abs(muon1.p4.Phi()-3.1)>1.e-2 and abs(muon2.p4.Phi()-3.1)>1.e-2) or (abs(muon1.p4.Phi()-3.1)<=1.e-2 and abs(muon2.p4.Phi()-3.1)<=1.e-2)
                                    else min(abs(muon1.p4.Phi()), abs(muon2.p4.Phi())))
                            #dPhi = abs(muon1.p4.Phi()-muon2.p4.Phi())

                            candidate = Candidate(muon1, imuon, dPt+dEta+dPhi)
#                            if max(dPt, dEta,)<=0.1:
#                                candidate = Candidate(muon1, imuon, dPhi)
#                                candidates.append(candidate)
#                            candidate = Candidate(muon1, imuon, dPhi+dEta+dPt)
                            candidates.append(candidate)
                        
                        if len(candidates):
                            candidates.sort(key = lambda candidate: candidate.difference)
#                            if candidates[0].difference>=min_dPhi:
                            if True:
                                muon2 = candidates[0].muon

                                Pt1.append(muon1.p4.Pt())
                                Pt2.append(muon2.p4.Pt())

                                Eta1.append(muon1.p4.Eta())
                                Eta2.append(muon2.p4.Eta())
                                        
                                Phi1.append(muon1.p4.Phi())
                                Phi2.append(muon2.p4.Phi())
                                matched += 1
                        total += 1
                        
                        #event2.muons.pop(candidates[0].position)

                else:
                    for muon2 in event2.muons:
                        candidates = []
                        for imuon, muon1 in enumerate(event1.muons):
                            #dPt = (abs(muon1.p4.Pt()-muon2.p4.Pt()) if muon1.p4.Pt()<10 and muon2.p4.Pt()<10 else
                            #        min(abs(muon1.p4.Pt()-muon2.p4.Pt()), 1))
                            dPt = abs(muon1.p4.Pt()-muon2.p4.Pt())
                            dPt /=10. #make the value be of the same order as the others.
                            dEta = abs(muon1.p4.Eta()-muon2.p4.Eta())
                            dPhi = (abs(muon1.p4.Phi()-muon2.p4.Phi()) if (abs(muon1.p4.Phi()-3.1)>1.e-2 and abs(muon2.p4.Phi()-3.1)>1.e-2) or (abs(muon1.p4.Phi()-3.1)<=1.e-2 and abs(muon2.p4.Phi()-3.1)<=1.e-2)
                                    else min(abs(muon1.p4.Phi()), abs(muon2.p4.Phi())))
#                            dPhi = abs(muon1.p4.Phi()-muon2.p4.Phi())

                            candidate = Candidate(muon1, imuon, dPt+dEta+dPhi)
#                            if max(dPt, dEta)<=0.2:
#                                candidate = Candidate(muon1, imuon, dPhi)
#                                candidates.append(candidate)
#                            candidate = Candidate(muon1, imuon, dPhi+dEta+dPt)
                            candidates.append(candidate)

                        if len(candidates):
                            candidates.sort(key = lambda candidate: candidate.difference)
#                            if candidates[0].difference>=min_dPhi:
                            if True:
                                muon1 = candidates[0].muon

                                Pt1.append(muon1.p4.Pt())
                                Pt2.append(muon2.p4.Pt())

                                Eta1.append(muon1.p4.Eta())
                                Eta2.append(muon2.p4.Eta())
                                    
                                Phi1.append(muon1.p4.Phi())
                                Phi2.append(muon2.p4.Phi())
                                matched += 1
                        total += 1
                        
                        #event1.muons.pop(candidates[0].position)

                #Events2.remove(event2)
                break

    
#    Pt1 = array('f', tuple(Pt1))
#    Eta1 = array('f', tuple(Eta1))
#    Phi1 = array('f', tuple(Phi1))

#    Pt2 = array('f', tuple(Pt2))
#    Eta2 = array('f', tuple(Eta2))
#    Phi2 = array('f', tuple(Phi2))
    return (Pt1, Pt2, Eta1, Eta2, Phi1, Phi2, matched, total)


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
