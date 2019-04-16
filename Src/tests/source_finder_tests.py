import os
import shutil
import random
import math
import xml.etree.ElementTree as Et

import gammalib
import ctools
import cscripts

import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("..") # ADDED TO BE EXECUTED AS A SCRIPT
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


def generate_xml_files(flow, n, start_coords, radius):
    os.chdir("../Tests")
    shutil.rmtree("Models")
    os.mkdir("Models")

    true_coords = []

    for i in range(0, n-1):
        num = padded_number(i, n)
        src_coords = random_coords(start_coords, radius)
        shutil.copy("default.xml", "Models/sigma4_" + padded_number(i, n) + ".xml")
        update_xml_file("Models/sigma4_" + num + ".xml", flow, src_coords)
        true_coords.append(src_coords)
        print("Generated Models/sigma4_" + num + ".xml, with flow =", flow, "and source coordinates =", src_coords)

    shutil.copy("default.xml", "Models/sigma4_" + str(n-1) + ".xml")
    update_xml_file("Models/sigma4_" + str(n-1) + ".xml", 0.001, start_coords)
    true_coords.append(None)
    print("Generated Models/sigma4_" + str(n-1) + ".xml, with flow =", 0.001, "and source coordinates =", start_coords)

    Util.write_list_to_json_file(true_coords, "Models/coordinates.json")

    os.chdir("../Src")
    return true_coords


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
    sim['inmodel'] = 'Models/sigma4_' + i + '.xml'
    sim['outevents'] = 'Events/' + i + '.fits'
    sim.execute()
    print("Generated Events/" + i + ".fits")


def generate_skymap(i):
    smap = ctools.ctskymap()
    smap['inobs'] = 'Events/' + i + '.fits'
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
    smap['outmap'] = 'Skymaps/' + i + '.fits'
    smap.execute()
    print("Generated Skymaps/" + i + ".fits")


def generate_skymaps(n):
    os.chdir("../Tests")
    shutil.rmtree("Events")
    shutil.rmtree("Skymaps")
    os.mkdir("Events")
    os.mkdir("Skymaps")

    for model in sorted(os.listdir('./Models')):
        if model == "coordinates.json":
            continue
        i = padded_number(int(model[0:-4].split('_')[1]), n)
        generate_events(i)
        generate_skymap(i)

    os.chdir("../Src")


def generate_data(flow, n, start_coords, radius):
    generate_xml_files(flow, n, start_coords, radius)
    generate_skymaps(n)


def source_finder_tests(flow=2.0, n=10, start_coords=(221, 46), radius=1):
    os.chdir("../")
    # generate_data(flow, n, start_coords, radius)
    true_coords = Util.read_list_from_json_file("../Tests/Models/coordinates.json")
    sf = SourceFinder("conf.json")
    # sf.compute_coords()
    computed_coords = Util.read_list_from_json_file("../Tests/Skymaps/computed_coordinates.json")
    print(true_coords)
    print(computed_coords)
    for i in range(0, n):
        if computed_coords[i] and true_coords[i]:
            distance = Util.distance_eu(true_coords[i], computed_coords[i])
            print(i, distance)
        else:
            print(i, "None")
    os.chdir("tests")


source_finder_tests(flow=1.0, n=100)
