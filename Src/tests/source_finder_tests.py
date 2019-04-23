import os
import shutil
import random
import math
import xml.etree.ElementTree as Et
import time

import ctools

import matplotlib.pyplot as plt
from astropy.wcs import WCS

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


def ctobssim(inmodel, outevents, ra=221.0, dec=46.0, rad=5.0, tmin='2020-01-01T00:00:00',
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


def ctskymap(inobs, outmap, xref=221.0, yref=46.0, proj='CAR', coordsys='CEL', binsz=0.02,
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
            model_filename = models_dir + model
            events_filename = events_dir + i + '.fits'
            skymap_filename = skymaps_dir + i + '.fits'
            ctobssim(inmodel=model_filename, outevents=events_filename, seed=seed)
            ctskymap(inobs=events_filename, outmap=skymap_filename)
            seed += 1


def generate_bkg_only_skymaps(n, bkg_model, events_dir, skymaps_dir, start_seed):
    seed = start_seed
    for i in range(0, n):
        num = padded_number(i, n)
        events_filename = events_dir + num + '.fits'
        skymap_filename = skymaps_dir + num + '.fits'
        ctobssim(inmodel=bkg_model, outevents=events_filename, seed=seed)
        ctskymap(inobs=events_filename, outmap=skymap_filename)
        seed += 1


def ctlike(inobs, inmodel, outmodel, caldb='prod2', irf='South_0.5h'):
    """"""
    ctl = ctools.ctlike()
    ctl['inobs'] = inobs
    ctl['caldb'] = caldb
    ctl['irf'] = irf
    ctl['inmodel'] = inmodel
    ctl['outmodel'] = outmodel
    ctl.execute()
    print("Generated " + outmodel)


def generate_fit_results(n, models_dir, events_dir, fit_results_dir):
    for model in sorted(os.listdir(models_dir)):
        if model.endswith('.xml'):
            i = padded_number(int(model[0:-4].split('_')[1]), n)
            model_filename = models_dir + model
            events_filename = events_dir + i + '.fits'
            fit_results_filename = fit_results_dir + i + '.xml'
            ctlike(inobs=events_filename, inmodel=model_filename, outmodel=fit_results_filename)


def generate_dirs(timestr, bkg_only=False):
    dirs = {'events': timestr + '/Events/',
            'skymaps': timestr + '/Skymaps/',
            'fit_results': timestr + '/FitResults/'}
    if not bkg_only:
        dirs['models'] = timestr + '/Models/'
    for d in dirs.values():
        os.makedirs(d)
    return dirs


def set_conf(skymaps_dir):
    conf_params = Util.read_list_from_json_file('conf.json')
    conf_params['dir'] = '../Tests/' + skymaps_dir
    Util.write_list_to_json_file(conf_params, 'conf.json')


def generate_src_data(start_model='default2.xml', flow=2.0, n=10, start_coords=(221, 46), radius=1, start_seed=1):
    os.chdir("../Tests")
    timestr = time.strftime("%Y%m%d-%H%M%S_") + start_model.split('.')[0] + '_' + str(flow)
    src_dirs = generate_dirs(timestr)
    generate_src_xml_files(flow, n, start_coords, radius, src_dirs['models'], start_model)
    generate_skymaps(n, src_dirs['models'], src_dirs['events'], src_dirs['skymaps'], start_seed)
    generate_fit_results(n, src_dirs['models'], src_dirs['events'], src_dirs['fit_results'])
    os.chdir("../Src")
    return src_dirs['skymaps']
    

def generate_bkg_only_data(bkg_only_model='background_only.xml', n=10, start_seed=1):
    os.chdir("../Tests")
    timestr = time.strftime("%Y%m%d-%H%M%S_") + bkg_only_model.split('.')[0]
    bkg_only_dirs = generate_dirs(timestr, bkg_only=True)
    generate_bkg_only_skymaps(n, bkg_only_model, bkg_only_dirs['events'], bkg_only_dirs['skymaps'], start_seed)
    os.chdir("../Src")
    return bkg_only_dirs['skymaps']


def analyse_src_data(computed_coords, dir_name):
    true_coords = Util.read_list_from_json_file('../Tests/' + dir_name + '/Models/coordinates.json')
    print("true coordinates:    ", true_coords)
    print("computed coordinates:", computed_coords)
    for i in range(0, len(computed_coords)):
        if computed_coords[i] and true_coords[i]:
            distance = Util.distance_eu(true_coords[i], computed_coords[i])
            print(i, distance)
        else:
            print(i, "None")


def analyse_bkg_only_data(computed_coords):
    print(computed_coords)


def analyse_dir(dir_name, bkg_only=False):
    computed_coords = Util.read_list_from_json_file('../Tests/' + dir_name + '/Skymaps/computed_coordinates.json')
    if not bkg_only:
        analyse_src_data(computed_coords, dir_name)
    else:
        analyse_bkg_only_data(computed_coords)


def run_source_finder():
    """"""
    sf = SourceFinder('conf.json')
    sf.compute_coords()


def plot_data(dir_name, measures):
    """"""
    skymaps_path = '../Tests/' + dir_name + '/Skymaps/'
    measures_dir = 'Measures/'
    path = skymaps_path + measures_dir
    data = {}

    for measure_file in sorted(os.listdir(path)):
        data[measure_file.split('.')[0]] = Util.read_list_from_json_file(path + measure_file)

    for skymap in sorted(os.listdir(skymaps_path)):
        if skymap.endswith('.fits'):
            wcs = WCS(skymaps_path + skymap)
            for measure in measures:
                for point in data[skymap.split('.fits')[0] + '_' + measure]:
                    eq = Util.from_pix_to_wcs(point['original'], wcs)
                    # eq = point['original']
                    parameters = Util.read_list_from_json_file('conf.json')
                    plt.plot(parameters['sigma_array'], point['values'],
                             label=str(round(float(eq[0]), 2)) + ", " + str(round(float(eq[1]), 2)))
                    plt.title(measure)
                    plt.legend()
                plt.figure()
    plt.show()


def run(bkg_only=False, dir_name=None, compute=False):
    """"""
    os.chdir("../")

    measures = ['isolatedness', 'intensity']

    if not dir_name:
        base_seed = random.randint(1, 1000000)
        if not bkg_only:
            # CHANGE SRC PARAMETERS HERE
            skymaps_dir = generate_src_data(flow=2.0, n=5, start_seed=base_seed)
        else:
            # CHANGE BKG_ONLY PARAMETERS HERE
            skymaps_dir = generate_bkg_only_data(n=5, start_seed=base_seed)
        set_conf(skymaps_dir)
        run_source_finder()
        dir_name = Util.read_list_from_json_file('conf.json')['dir'].split('/')[-3]

    elif compute:
        skymaps_dir = dir_name+'/Skymaps'
        set_conf(skymaps_dir)
        run_source_finder()

    analyse_dir(dir_name, bkg_only)
    plot_data(dir_name, measures)

    os.chdir("tests")


# UNCOMMENT WHAT YOU NEED
# run()  # GENERATE NEW SRC SKYMAPS
run(dir_name='20190420-091900_default_2.0', compute=True)  # ANALYSE SRC SKYMAPS IN DIR_NAME/SKYMAPS/
# run(bkg_only=True, compute=False)  # GENERATE NEW BKG_ONLY SKYMAPS
# run(bkg_only=True, dir_name='20190420-111028_background_only', compute=True)  # ANALYSE BKG_ONLY SKYMAPS IN DIR_NAME/SKYMAPS
