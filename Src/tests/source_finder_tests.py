import os
import shutil
import random
import math
import xml.etree.ElementTree as Et
from source_finder import SourceFinder
import gammalib
import ctools
import cscripts
import matplotlib.pyplot as plt
import numpy as np
from util import Util


def random_coords(center):
    module = random.uniform(0, 1)
    angle = random.uniform(0, 2 * math.pi)
    x = center[0] + math.cos(angle) * module
    y = center[1] + math.sin(angle) * module
    return x, y


def update_xml_file(file_name, flow, coords):
    tree = Et.parse(file_name)
    root = tree.getroot()
    root[0][0][0].attrib['value'] = str(flow)
    root[0][1][0].attrib['value'] = str(coords[0])
    root[0][1][1].attrib['value'] = str(coords[1])
    tree.write(file_name)


def generate_events(i):
    sim = ctools.ctobssim()
    sim['ra'] = 221
    sim['dec'] = 46
    sim['rad'] = 5.0
    sim['tmin'] = '2020-01-01T00:00:00'
    sim['tmax'] = '2020-01-01T00:15:00'
    sim['emin'] = 0.1
    sim['emax'] = 100.0
    sim['caldb'] = 'prod2'
    sim['irf'] = 'South_0.5h'
    sim['inmodel'] = 'Models/sigma4_'+str(i)+'.xml'
    sim['outevents'] = 'Events/'+str(i)+'.fits'
    sim.execute()
    print("Generated Events/" + str(i) + ".fits")


def generate_skymap(i):
    smap = ctools.ctskymap()
    smap['inobs'] = 'Events/'+str(i)+'.fits'
    smap['xref'] = 221
    smap['yref'] = 46
    smap['proj'] = 'CAR'
    smap['coordsys'] = 'CEL'
    smap['binsz'] = 0.02
    smap['nxpix'] = 200
    smap['nypix'] = 200
    smap['emin'] = 0.1
    smap['emax'] = 100
    smap['bkgsubtract'] = 'NONE'
    smap['outmap'] = 'Skymaps/' + str(i) + '.fits'
    smap.execute()
    print("Generated Skymap/" + str(i) + ".fits")


def generate_xml_files(n, central_coords, flow):
    os.chdir("../Tests")
    shutil.rmtree("Models")
    os.mkdir("Models")

    true_coords = []

    for i in range(0, n):
        coords = random_coords(central_coords)
        shutil.copy("default.xml", "Models/sigma4_" + str(i) + ".xml")
        update_xml_file("Models/sigma4_" + str(i) + ".xml", flow, coords)
        true_coords.append(coords)
        print("Generated Models/sigma4_" + str(i) + ".xml, with flow =", flow, "and coordinates =", coords)

    Util.write_list_to_json_file(true_coords, "Models/coordinates.json")

    os.chdir("../Src")
    return true_coords


def generate_skymaps():
    os.chdir("../Tests")
    shutil.rmtree("Events")
    shutil.rmtree("Skymaps")
    os.mkdir("Events")
    os.mkdir("Skymaps")

    for model in sorted(os.listdir('./Models')):
        if model == "coordinates.json":
            continue
        i = int(model[0:-4].split('_')[1])
        generate_events(i)
        generate_skymap(i)

    os.chdir("../Src")


def generate_data(n, coords, flow):
    generate_xml_files(n, coords, flow)
    generate_skymaps()


def source_finder_tests(flow=2, n=10):
    os.chdir("../")
    # generate_data(n, (221, 46), flow)
    true_coords = Util.read_list_from_json_file("../Tests/Models/coordinates.json")
    sf = SourceFinder("conf.json")
    sf.compute_coords()
    computed_coords = Util.read_list_from_json_file("../Tests/Skymaps/computed_coordinates.json")
    print(true_coords)
    print(computed_coords)
    for i in range(0, n):
        if computed_coords[i]:
            print(i, Util.distance_eu(true_coords[i], computed_coords[i]))
        else:
            print(i, "None")
    os.chdir("tests")


source_finder_tests(n=100)
# plt.show()
