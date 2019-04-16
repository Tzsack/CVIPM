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


def generate_events(i, n):
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
    sim['inmodel'] = 'Models/sigma4_' + padded_number(i, n) + '.xml'
    sim['outevents'] = 'Events/' + padded_number(i, n) + '.fits'
    sim.execute()
    print("Generated Events/" + padded_number(i, n) + ".fits")


def generate_skymap(i, n):
    smap = ctools.ctskymap()
    smap['inobs'] = 'Events/' + padded_number(i, n) + '.fits'
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
    smap['outmap'] = 'Skymaps/' + padded_number(i, n) + '.fits'
    smap.execute()
    print("Generated Skymap/" + padded_number(i, n) + ".fits")


def generate_xml_files(n, central_coords, flow):
    os.chdir("../Tests")
    shutil.rmtree("Models")
    os.mkdir("Models")

    true_coords = []

    for i in range(0, n-1):
        coords = random_coords(central_coords)
        shutil.copy("default.xml", "Models/sigma4_" + padded_number(i, n) + ".xml")
        update_xml_file("Models/sigma4_" + padded_number(i, n) + ".xml", flow, coords)
        true_coords.append(coords)
        print("Generated Models/sigma4_" + padded_number(i, n) + ".xml, with flow =", flow, "and coordinates =", coords)

    shutil.copy("default.xml", "Models/sigma4_" + str(n-1) + ".xml")
    update_xml_file("Models/sigma4_" + str(n-1) + ".xml", 0.001, central_coords)
    true_coords.append(None)
    print("Generated Models/sigma4_" + str(n-1) + ".xml, with flow =", 0.001, "and coordinates =", central_coords)

    Util.write_list_to_json_file(true_coords, "Models/coordinates.json")

    os.chdir("../Src")
    return true_coords


def generate_skymaps(n):
    os.chdir("../Tests")
    shutil.rmtree("Events")
    shutil.rmtree("Skymaps")
    os.mkdir("Events")
    os.mkdir("Skymaps")

    for model in sorted(os.listdir('./Models')):
        if model == "coordinates.json":
            continue
        i = int(model[0:-4].split('_')[1])
        generate_events(i, n)
        generate_skymap(i, n)

    os.chdir("../Src")


def generate_data(n, coords, flow):
    generate_xml_files(n, coords, flow)
    generate_skymaps(n)


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


def source_finder_tests(flow=2.0, n=10):
    os.chdir("../")
    generate_data(n, (221, 46), flow)
    true_coords = Util.read_list_from_json_file("../Tests/Models/coordinates.json")
    sf = SourceFinder("conf.json")
    sf.compute_coords()
    computed_coords = Util.read_list_from_json_file("../Tests/Skymaps/computed_coordinates.json")
    print(true_coords)
    print(computed_coords)
    for i in range(0, n):
        if computed_coords[i] and true_coords[i]:
            print(i, Util.distance_eu(true_coords[i], computed_coords[i]))
        else:
            print(i, "None")
    os.chdir("tests")


source_finder_tests(n=100, flow=1.0)