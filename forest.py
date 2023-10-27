import sys
import argparse
import time
import numpy as np
from nu_reactions.analytic_generator import analytic_generator
from nu_reactions.nu_event_generator import nu_event_generator
from nu_reactions.inverse_beta_decay import inverse_beta_decay
from nu_osc_models.no_osc import no_osc
from nu_osc_models.MSW_normal import MSW_normal
from nu_osc_models.MSW_inverted import MSW_inverted
from detectors.super_kamiokande import super_kamiokande
from analytic_spectra import analytic_spectra
from sn_spectra import SNspectra

forest_version = "1.0.0"
ev_format = "{0:d}\t{1:e}\t{2:e}\t{3:e}\t{4:e}\t{5:e}\t{6:e}\t{7:e}\t{8:e}\t{9:d}\t{10:d}\n"  #ev_no, time[s], en_ene [MeV], nu_eve [MeV], theta, phi, x [cm], y [cm], z [cm], event id, fv 


def write_events(ev_list, output_file):
    text = ""
    text += "#FOREST version: " + forest_version + "\n"
    text += "#Python version: " + sys.version.replace('\n', '') + "\n"
    text += "#"
    for cmd in sys.argv:
        text += cmd + " "
    text += "\n"
    text += "#Seed: " + str(seed) + "\n"
    text += ""
    

    if(output_file != None):
        with open(output_file, "w") as f:
            f.write(text)
            for index, ev_line in enumerate(ev_list):
                f.write(ev_format.format(index, ev_line["time"], ev_line["ev_ene"], ev_line["nu_ene"],
                                          ev_line["theta"], ev_line["phi"], ev_line["x"], ev_line["y"], ev_line["z"], ev_line["id"], ev_line["fv"]))
    else:
         for index, ev_line in enumerate(ev_list):
                print(ev_format.format(index, ev_line["time"], ev_line["ev_ene"], ev_line["nu_ene"],
                                          ev_line["theta"], ev_line["phi"], ev_line["x"], ev_line["y"], ev_line["z"], ev_line["id"], ev_line["fv"]),end='')        



def gen_events(times):
    event_list = []
    pre_output_ev = 0
    for time1, time2 in zip(times[:-1], times[1:]):
        event_list.extend(detector.get_events(time1, time2))
        if(len(event_list)%100 == 0 and len(event_list) != pre_output_ev):
            pre_output_ev = len(event_list)
            print("Events number: ", len(event_list)//100*100)
    return event_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='forest',
        description='FOREST: Forecasting Events from Supernovae Theoretical modeling\
                                     This program predicts supernova neutrino events on earth conducting realistic detector simulations.')
    parser.add_argument('-model', choices=['analytic_rate', 'analytic_spectra', 'numerical_spectra'])
    parser.add_argument('-etot', '-total_energy', help='Total energy for the analytic supernova model in erg',required=False, default=1.0e53, type=float)
    parser.add_argument('-pns_m',  help='Proto-neutron star mass for the analytic supernova model in solar mass',required=False, default=1.4, type=float)
    parser.add_argument('-pns_r',  help='Proto-neutron star radius for the analytic supernova model in km',required=False, default=10.0, type=float)
    parser.add_argument('-gbeta', help='gbeta',required=False, default=3.0, type=float)
    parser.add_argument('-end_time', help='End time of neutrino emission for the analytic supernova model in second',required=False, default=100, type=float)
    parser.add_argument('-spectra_file', help='Neutrino spectra file',required=False, type=str)
    parser.add_argument('-distance', default=10.0, type=float, required=True)
    parser.add_argument('-detector',choices=['superk', 'hyperk'], type=str, required=True)
    parser.add_argument('-o', '-output', help='Output filename', type=str)
    parser.add_argument('-if', '-input_format', choices=['nakazato'], help='Input file format. Default is the Nakazato format', default='nakazato', type=str)
    parser.add_argument('-no', '-neutrino_oscillation', choices=['no', 'msw_normal', 'msw_inverted'], type=str, default='msw_normal')
    parser.add_argument('-seed', type=int, help="if None, the unix time is used.")
    parser.add_argument('-v', '-visualize', action='store_true', help='Do you want to see a result after generating.')
    args = parser.parse_args()
    
    if args.seed == None:
        seed = int(time.time())
    else:
        seed = args.seed
    
    np.random.seed(seed)

    print('##### FOREST #####')
    print('Model: ', args.model)
    print('Seed: ', seed)

    if args.model == 'analytic_rate':
        if args.detector == 'superk':
            volume = super_kamiokande.VOLUME
            gen = analytic_generator(args.pns_m, args.pns_r, args.gbeta, args.etot, volume, args.distance)
            detector = super_kamiokande(None, gen)
        elif argparse.detector == 'hyperk':
            print('Not implemented with hyperk yet.')
            sys.exit(0)
        
        start_time = gen.get_t0()
        end_time = args.end_time
        times = np.arange(start_time, end_time, 0.01)

    elif args.model == 'analytic_spectra':
        spectra = analytic_spectra(args.pns_m, args.pns_r, args.gbeta, args.etot, args.distance, args.end_time, 1, 300, 100)
        times = spectra.get_times()
        flux_earth =  no_osc(spectra)
        if args.detector == 'superk':
            ibd = inverse_beta_decay()
            gen = nu_event_generator(flux_earth, [ibd], [super_kamiokande.PROTONS], 200)
            detector = super_kamiokande(flux_earth, gen)
        elif args.dector == 'hyperk':
            print('Not implemented with hyperk yet.')
            sys.exit(0)
            
    elif args.model == 'numerical_spectra':
        spectra = SNspectra(args.spectra_file, args.distance)
        times = spectra.get_times()
        if args.no == 'no':
            flux_earth =  no_osc(spectra)
        elif args.no == 'msw_normal':
            flux_earth =  MSW_normal(spectra)
        elif args.no == 'msw_inverted':
            flux_earth = MSW_inverted(spectra)

        if args.detector == 'superk':
            ibd = inverse_beta_decay()
            gen = nu_event_generator(flux_earth, [ibd], [super_kamiokande.PROTONS], 200)
            detector = super_kamiokande(flux_earth, gen)
        elif args.dector == 'hyperk':
            print('Not implemented with hyperk yet.')
            sys.exit(0)


    print('Generting events:\n')
    ev_list = gen_events(times)
    write_events(ev_list,args.o)
    
    sys.exit()