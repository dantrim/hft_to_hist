#!/usr/bin/env python

from __future__ import print_function # just incase
import os, sys
import argparse

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)
r.TH1F.__init__._creates = False
r.TH1.__init__._creates = False

class Region :
    def __init__(self) :
        self.name = ""
        self.tcut = ""
    def __str__(self) :
        return " {} : {}".format(self.name, self.tcut)

class Blinder :
    def __init__(self) :
        self.h_counts = None
        self.h_disc = None

def get_regions(args) :

    reg_dict = {}

    trigger = "(( year == 2015 && trig_tight_2015 == 1 ) || ( year == 2016 && trig_tight_2016 == 1 ))"

    reg_dict["crtt"] = "{trig} && nBJets>=2 && mbb>140 && NN_d_tt>1.5".format(trig=trigger)
    reg_dict["crwt"] = "{trig} && nBJets>=2 && mbb>140 && NN_d_tt<2.5 && NN_d_wt>2.2".format(trig=trigger)
    reg_dict["sr"] = "{trig} && mbb>110 && mbb<140 && nBJets>=2 && mt2_bb>65".format(trig=trigger)

    regions = []
    for reg_name, tcut in reg_dict.iteritems() :
        r = Region()
        r.name = reg_name
        r.tcut = tcut
        regions.append(r)

    if args.region != "" :
        requested_regions = [x.strip() for x in args.region.split(",")]
        loaded_regions = [reg.name for reg in regions]

        request_ok = True
        for req_reg in requested_regions :
            if req_reg not in loaded_regions :
                print('ERROR requested region (={}) not found in list of configured regions'.format(req_reg))
                request_ok = False
        if not request_ok :
            sys.exit()
        
        tmp = []
        for reg_req in requested_regions :
            for reg in regions :
                if reg.name == reg_req :
                    tmp.append(reg)
        regions = tmp
            
    return regions

def find_unique_process_and_sys_names(hft_root, verbose) :

    '''
    Assumes the basic structure of a HFT where it contains
    only TTrees whose names are <Process>_<Systematic>.
    Get a list of the Process names.

    Args:
        hft_root : opened TFile of the HFT

    Returns:
        python list of strings
    '''

    process_names = set()
    sys_names = set()
    name_of_first_nominal_tree = "" # it will actually be the last
    for key in hft_root.GetListOfKeys() :
        key_name = key.GetName()
        key_split = key_name.split('_')
        process_name = key_split[0]
        sys_name = '_'.join(key_split[1:])

        if process_name is not "Data" and sys_name == "CENTRAL" :
            name_of_first_nominal_tree = key.GetName()

        process_names.add(process_name)
        sys_names.add(sys_name)

    # we need always the CENTRAL tree
    if 'CENTRAL' not in sys_names :
        print('ERROR Did not find the nominal systematic in the input file')
        sys.exit()

    # the weight-based systematics are in the CENTRAL TTree, grab the first non-Data TTree to get their names
    weight_sys_names = set()
    nom_tree = hft_root.Get(name_of_first_nominal_tree)
    if not nom_tree :
        print('ERROR Loading nominal TTree (={}) from input (attempting to grab weight-based sys'.format(name_of_first_nominal_tree))
        sys.exit()
    for branch_name in nom_tree.GetListOfBranches() :
        branch_name = branch_name.GetName()
        if not branch_name.startswith('syst_') : continue
        weight_sys_names.add(branch_name)

    if not process_names :
        print('ERROR Did not find any processes in the input file')
        sys.exit()

    if not sys_names :
        print('ERROR Did not find any systematics in the input file')
        sys.exit()

    n_processes = len(process_names)
    print('Found {} unique processes'.format(n_processes))
    if verbose :
        print(25 * '- ')
        for iprocess, process in enumerate(process_names) :
            print('[{0:02d}/{1:02d}] {2}'.format(iprocess+1, n_processes, process))

    n_sys = len(sys_names)
    print('Found {} unique shape systematics'.format(n_sys))
    if verbose :
        print(25 * '- ')
        for isys, sys_name in enumerate(sys_names) :
            print('[{0:02d}/{1:02d}] {2}'.format(isys+1, n_sys, sys_name))

    n_w_sys = len(weight_sys_names)
    print('Found {} unique weight systematics'.format(n_w_sys))
    if verbose :
        print(25 * '- ')
        for isys, sys_name in enumerate(weight_sys_names) :
            print('[{0:02d}/{1:02d}] {2}'.format(isys+1, n_w_sys, sys_name))

    for w_sys in weight_sys_names :
        sys_names.add(w_sys)

    return process_names, sys_names

def get_user_sys_names(sys_list, selected, verbose) :

    if selected == '' :
        return sys_list

    sys_selected = [s.strip() for s in selected.split(',')]

    sys_ok = True
    sys_out = set()
    for s in sys_selected :
        if s not in sys_list :
            print('ERROR A requested systematic (={}) was not found in input'.format(s))
            sys_ok = False
        else :
            sys_out.add(s)

    if not sys_ok :
        sys.exit()

    return sys_out

def get_tree_name(process_name, sys_name) :

    tree_name = process_name

    if sys_name.startswith("syst_") :
        tree_name += "_CENTRAL"
    else :
        tree_name += "_{}".format(sys_name)

    return tree_name

def make_histogram_file(hft_root, process_names, sys_names, regions, args) :

    outfilename = "test_ws.root"
    rout = r.TFile.Open(outfilename, 'RECREATE') # should think of a way to append new things to the file if not found already

    n_reg = len(regions)
    n_proc = len(process_names)
    n_sys = len(sys_names)
    for iregion, region in enumerate(regions) :

        # first move to the top level
        rout.cd()

        print('Region [{0:02d}/{1:02d}] {2}'.format(iregion+1, n_reg, region.name))

        region_dir = r.TDirectoryFile(region.name, region.name, "", rout)
        region_dir.cd()

        histos_for_blinding = {}

        for iprocess, process_name in enumerate(process_names) :
            print('  - Process [{0:02d}/{1:02d}] {2}'.format(iprocess+1, n_proc, process_name))

            # blind
            if process_name == "Data" and region.name == "sr" : continue

            process_dir = r.TDirectoryFile(process_name, process_name, "", region_dir)
            process_dir.cd()

            if region.name == "sr" :
                histos_for_blinding[process_name] = {}

            # don't add data histogram to SR (we could perhaps add the sum of the bkg so that there is at least a histogram)
            for isys, sys_name in enumerate(sys_names) :

                if process_name == "Data" and sys_name == "CENTRAL" :
                    print('     > Systematic [1/1] CENTRAL')
                else :
                    print('     > Systematic [{0:02d}/{1:02d}] {2}'.format(isys+1, n_sys, sys_name))

                if process_name == "Data" and sys_name != "CENTRAL" : continue

                sys_dir = r.TDirectoryFile(sys_name, sys_name, "", process_dir)
                sys_dir.cd()


                h = r.TH1F("h_{}_{}_{}_counts".format(region.name, process_name, sys_name), "", 4, 0, 4)
                h.Sumw2()
                sel = r.TCut("1")
                cut = region.tcut


                weight = "eventweightbtag"
                cutstr = "({cut}) * {weight} * 36.1".format(cut = cut, weight = weight)
                if process_name == "Data" :
                    cutstr = "({cut})".format(cut = cut)

                cut = r.TCut(cutstr)

                tree_to_get = get_tree_name(process_name, sys_name)
                tree = hft_root.Get(tree_to_get)
                cmd = "isMC>>{}".format(h.GetName())
                tree.Draw(cmd, cut * sel, "goff")

                h_out = r.TH1F("h_{}_{}_{}".format(region.name, process_name, sys_name), "{}_{}_{};;".format(region.name, process_name, sys_name), 1, 0, 1)
                h_out.Sumw2()
                integral = h.Integral()

                h.Delete()

                h_out.SetBinContent(1, integral)
                h_out.Draw("goff")

                sys_dir.cd()
                h_out.Write()

                if region.name == "sr" :
                    # add finely binned histogram for the NN discriminant we wish to fit
                    h_nn = r.TH1F("h_{}_{}_{}_disc".format(region.name, process_name, sys_name), "{}_{}_{}_disc;;".format(region.name, process_name, sys_name), 300, -30, 30)
                    h_nn.Sumw2()

                    cmd = "NN_d_hh>>{}".format(h_nn.GetName())
                    tree.Draw(cmd, cut * sel, "goff")
                    h_nn.Write()

                sys_dir.Write()

            process_dir.Write()

        if region.name == "sr" :

            # add blinded data
            process_dir_data = r.TDirectoryFile("Data", "Data", "", region_dir)
            process_dir_data.cd()


            sys_dir_data = r.TDirectoryFile("CENTRAL", "CENTRAL", "", process_dir_data)
            sys_dir_data.cd()
            
            h_counts_data = r.TH1F("h_{}_{}_CENTRAL".format("sr", "Data"), "{}_{}_CENTRAL;;".format("sr", "Data"), 1, 0, 1)
            h_counts_data.Sumw2()
            h_disc_data = r.TH1F("h_{}_{}_CENTRAL_disc".format("sr", "Data"), "{}_{}_CENTRAL_disc;;".format("sr", "Data"), 300, -30, 30)
            h_disc_data.Sumw2()
            
            n_counts_data = 0
            rdir = rout.Get("sr")
            rdir.cd()
            for iprocess, process_name in enumerate(process_names) :

                if "Data" in process_name :
                    continue

                pdir = rdir.Get(process_name)
                pdir.cd()
                sdir = pdir.Get("CENTRAL")
                sdir.cd()

                hname = "h_{}_{}_{}".format("sr", process_name, "CENTRAL")
                h_counts = sdir.Get(hname)
                n_counts_data += h_counts.Integral()

                hname += "_disc"
                h_disc = sdir.Get(hname)
                h_disc_data.Add(h_disc_data, h_disc)

            h_counts_data.SetBinContent(1,n_counts_data)

            h_counts_data.Write()
            h_disc_data.Write()
            sys_dir_data.Write()

        region_dir.Write()

    rout.Write()
                

def main(args) :

    rfile = r.TFile.Open(args.input)
    print('Opening {}...'.format(args.input))

    process_names, sys_names = find_unique_process_and_sys_names(rfile, args.verbose)

    if args.print_sys or args.print_processes :
        if args.print_processes :
            print(50 * '=')
            print('Processes found in provided input:')
            print(25 * '- ')
            n_proc = len(process_names)
            for iproc, proc_name in enumerate(process_names) :
                print('[{0:02d}/{1:02d}] {2}'.format(iproc+1, n_proc, proc_name))
            print(50 * '=')
        if args.print_sys :
            print(50 * '=')
            print('Systematics found in provided input:')
            print(25 * '- ')
            n_sys = len(sys_names)
            for isys, sys_name in enumerate(sys_names) :
                print('[{0:02d}/{1:02d}] {2}'.format(isys+1, n_sys, sys_name))
            print(50 * '=')
        sys.exit()

    sys_names = get_user_sys_names(sys_names, args.sys, args.verbose)

    regions = get_regions(args)

    make_histogram_file(rfile, process_names, sys_names, regions, args)


if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description = "Convert HistFitter Trees to Histogram Files for Workspace Building")
    parser.add_argument("-i", "--input", required = True,
        help = "Input HistFitter tree file to convert"
    )
    parser.add_argument("--sys", default = "",
        help = "Select specific systematics to use (comma-separated-list)"
    )
    parser.add_argument("--region", default = "",
        help = "Select a specific region to consider (comma-separated-list)"
    )
    parser.add_argument('-v', '--verbose', default = False, action = 'store_true',
        help = 'Print out more information as it runs'
    )
    parser.add_argument('--print-sys', default = False, action = 'store_true',
        help = 'Only print out the systematics in the input'
    )
    parser.add_argument('--print-processes', default = False, action = 'store_true',
        help = 'Only print out the processes in the input'
    )
    parser.add_argument('--suffix', default = '',
        help = 'Provide a suffix to append to the output(s)'
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input) :
        print('ERROR Provided input (={}) not found'.format(args.input))
        sys.exit()

    main(args)


