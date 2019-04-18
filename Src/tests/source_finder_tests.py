import os
import shutil
import random
import math
import xml.etree.ElementTree as Et
import time

import gammalib
import ctools
import cscripts

import matplotlib.pyplot as plt
import numpy as np

import sys

sys.path.append("..")  # ADDED TO BE EXECUTED AS A SCRIPT
from source_finder import SourceFinder
from util import Util


def random_coords(start_coords, radius):
    module = random.uniform(0, radius)
    angle = random.uniform(0, 2 * math.pi)
    ra = start_coords[0] + math.cos(angle) * module
    dec = start_coords[1] + math.sin(angle) * module
    return ra, dec


def padded_number(n, max_n):
    digits = math.ceil(math.log10(max_n))
    n_digits = math.ceil(math.log10(n + 1))
    if n_digits == 0:
        n_digits = 1
    result = ""
    for i in range(0, digits - n_digits):
        result += "0"
    result += str(n)
    return result


def update_xml_file(file_name, flow, src_coords):
    tree = Et.parse(file_name)
    root = tree.getroot()
    root[0][0][0].attrib['value'] = str(flow)
    root[0][1][0].attrib['value'] = str(src_coords[0])
    root[0][1][1].attrib['value'] = str(src_coords[1])
    tree.write(file_name)


def generate_src_xml_files(flow, n, start_coords, radius, models_dir, start_model):
    """"""
    true_coords = []

    for i in range(0, n):
        num = padded_number(i, n)
        src_coords = random_coords(start_coords, radius)
        model_name = models_dir + 'sigma4_' + num + '.xml'
        shutil.copy(start_model, model_name)
        update_xml_file(model_name, flow, src_coords)
        true_coords.append(src_coords)
        print("Generated " + model_name + " with flow =", flow, "and source coordinates =", src_coords)

    Util.write_list_to_json_file(true_coords, models_dir + "coordinates.json")

    return true_coords


def generate_events(inmodel, outevents, ra=221.0, dec=46.0, rad=5.0, tmin='2020-01-01T00:00:00',
                    tmax='2020-01-01T00:15:00', emin=0.1, emax=100.0, caldb='prod2', irf='South_0.5h', seed=1):
    sim = ctools.ctobssim()
    sim['ra'] = ra
    sim['dec'] = dec
    sim['rad'] = rad
    sim['tmin'] = tmin
    sim['tmax'] = tmax
    sim['emin'] = emin
    sim['emax'] = emax
    sim['caldb'] = caldb
    sim['irf'] = irf
    sim['inmodel'] = inmodel
    sim['outevents'] = outevents
    sim['seed'] = seed
    sim.execute()
    print("Generated " + outevents + " using seed: " + str(seed))


def generate_skymap(inobs, outmap, xref=221.0, yref=46.0, proj='CAR', coordsys='CEL', binsz=0.02,
                    nxpix=200, nypix=200, emin=0.1, emax=100.0, bkgsubtract='NONE'):
    smap = ctools.ctskymap()
    smap['inobs'] = inobs
    smap['xref'] = xref
    smap['yref'] = yref
    smap['proj'] = proj
    smap['coordsys'] = coordsys
    smap['binsz'] = binsz
    smap['nxpix'] = nxpix
    smap['nypix'] = nypix
    smap['emin'] = emin
    smap['emax'] = emax
    smap['bkgsubtract'] = bkgsubtract
    smap['outmap'] = outmap
    smap.execute()
    print("Generated " + outmap)


def generate_skymaps(n, models_dir, events_dir, skymaps_dir, start_seed):
    """"""
    seed = start_seed
    for model in sorted(os.listdir(models_dir)):
        if model.endswith('.xml'):
            i = padded_number(int(model[0:-4].split('_')[1]), n)
            model_filename = models_dir + 'sigma4_' + i + '.xml'
            event_filename = events_dir + i + '.fits'
            skymap_filename = skymaps_dir + i + '.fits'
            generate_events(inmodel=model_filename, outevents=event_filename, seed=seed)
            generate_skymap(inobs=event_filename, outmap=skymap_filename)
            seed += 1


def generate_bkg_only_skymaps(n, bkg_model, events_dir, skymaps_dir, start_seed):
    seed = start_seed
    for i in range(0, n):
        num = padded_number(i, n)
        event_filename = events_dir + num + '.fits'
        skymap_filename = skymaps_dir + num + '.fits'
        generate_events(inmodel=bkg_model, outevents=event_filename, seed=seed)
        generate_skymap(inobs=event_filename, outmap=skymap_filename)
        seed += 1


def generate_dirs():
    """"""
    timestr = time.strftime("%Y%m%d-%H%M%S")
    dirs = {'models': timestr + '/Models/', 'events': timestr + '/Events/',
            'skymaps': timestr + '/Skymaps/'}
    for d in dirs.values():
        os.makedirs(d)
    return dirs


def generate_data(flow, n, start_coords, radius, start_model, start_seed):
    """"""
    os.chdir("../Tests")
    dirs = generate_dirs()
    if start_model != 'background_only.xml':
        generate_src_xml_files(flow, n, start_coords, radius, dirs['models'], start_model)
        generate_skymaps(n, dirs['models'], dirs['events'], dirs['skymaps'], start_seed)

    else:
        generate_bkg_only_skymaps(n, start_model, dirs['events'], dirs['skymaps'], start_seed)
    os.chdir("../Src")
    return dirs


def set_conf(skymaps_dir):
    conf_params = Util.read_list_from_json_file('conf.json')
    conf_params['dir'] = '../Tests/' + skymaps_dir
    Util.write_list_to_json_file(conf_params, 'conf.json')


def analyse_data(dir_name, bkg_only=False):
    os.chdir("../")
    computed_coords = Util.read_list_from_json_file('../Tests/' + dir_name + '/Skymaps/computed_coordinates.json')
    if not bkg_only:
        true_coords = Util.read_list_from_json_file('../Tests/' + dir_name + '/Models/coordinates.json')
        print(true_coords)
        print(computed_coords)
        for i in range(0, len(computed_coords)):
            if computed_coords[i] and true_coords[i]:
                distance = Util.distance_eu(true_coords[i], computed_coords[i])
                print(i, distance)
            else:
                print(i, "None")
    else:
        print(computed_coords)
    os.chdir("tests")


def source_finder_tests(flow=2.0, n=10, start_coords=(221, 46), radius=1, start_model='default.xml', start_seed=1,
                        generate=True):
    os.chdir("../")
    if generate:
        dirs = generate_data(flow, n, start_coords, radius, start_model, start_seed)
        set_conf(dirs['skymaps'])
    sf = SourceFinder('conf.json')
    sf.compute_coords()
    os.chdir("tests")

    timestr = Util.read_list_from_json_file('../conf.json')['dir'].split('/')[-3]
    if start_model != 'background_only.xml':
        analyse_data(timestr)
    else:
        analyse_data(timestr, bkg_only=True)


base_seed = random.randint(1, 1000000)
# source_finder_tests(n=5, start_seed=base_seed, start_model='background_only.xml', generate=True)
source_finder_tests(n=5, start_seed=base_seed, start_model='default.xml', generate=False, flow=1.0)
plt.show()
