'''
    A utility script to help the analysis of imaging calorimeter data

    Author: Chao Peng (ANL)
    Date: 04/30/2021

    Added all mc particle info to obtain decaying particles
    Author: Jihee Kim (ANL)
    Date: 08/06/2021
'''

import os
import ROOT
import numpy as np
import pandas as pd
import matplotlib
import DDG4
from ROOT import gROOT, gInterpreter
from lxml import etree as ET


class ReadoutDecoder:
    def __init__(self, compact, readout):
        self.readouts = self.getReadouts(compact)
        self.changeReadout(readout)

    def changeReadout(self, readout):
        self.fieldsmap = self.decomposeIDs(self.readouts[readout])

    def get(self, idvals, field):
        start, width = self.fieldsmap[field]
        if width >= 0:
            return np.bitwise_and(np.right_shift(idvals, start), (1 << width) - 1)
        # first bit is sign bit
        else:
            width = abs(width) - 1
            vals = np.bitwise_and(np.right_shift(idvals, start), (1 << width) - 1)
            return np.where(np.bitwise_and(np.right_shift(idvals, start + width), 1), vals - (1 << width), vals)

    def mask(self, field):
        start, width = self.fieldsmap[field]
        return np.uint64((2**abs(width) - 1) << start)

    def decode(self, idvals):
        return {field: self.get(idvals, field) for field, _ in self.fieldsmap.items()}

    @staticmethod
    def getReadouts(path):
        res = dict()
        ReadoutDecoder.__getReadoutsRecur(path, res)
        return res

    @staticmethod
    def __getReadoutsRecur(path, res):
        if not os.path.exists(path):
            print('Xml file {} not exist! Ignored it.'.format(path))
            return
        lccdd = ET.parse(path).getroot()
        readouts = lccdd.find('readouts')
        if readouts is not None:
            for readout in readouts.getchildren():
                ids = readout.find('id')
                if ids is not None:
                    res[readout.attrib['name']] = ids.text
        for child in lccdd.getchildren():
            if child.tag == 'include':
                root_dir = os.path.dirname(os.path.realpath(path))
                ReadoutDecoder.__getReadoutsRecur(os.path.join(root_dir, child.attrib['ref']), res)

    @staticmethod
    def decomposeIDs(id_str):
        res = dict()
        curr_bit = 0
        for field_bits in id_str.split(','):
            elements = field_bits.split(':')
            field_name = elements[0]
            bit_width = int(elements[-1])
            if len(elements) == 3:
                curr_bit = int(elements[1])
            res[field_name] = (curr_bit, bit_width)
            curr_bit += abs(bit_width)
        return res


# read from RDataFrame and flatten a given collection, return pandas dataframe
def flatten_collection(rdf, collection, cols=None):
    if not cols:
        cols = [str(c) for c in rdf.GetColumnNames() if str(c).startswith('{}.'.format(collection))]
    else:
        cols = ['{}.{}'.format(collection, c) for c in cols]
    if not cols:
        print('cannot find any branch under collection {}'.format(collection))
        return pd.DataFrame()
    # print(rdf.GetColumnNames())
    data = rdf.AsNumpy(cols)
    # flatten the data, add an event id to identify clusters from different events
    evns = []
    for i, vec in enumerate(data[cols[0]]):
        evns += [i]*vec.size()
    for n, vals in data.items():
        # make sure ints are not converted to floats
        typename = vals[0].__class__.__name__.lower()
        # default
        dtype = np.float64
        if 'unsigned int' in typename or 'unsigned long' in typename:
            dtype = np.uint64
        elif 'int' in typename or 'long' in typename:
            dtype = np.int64
        # print(n, typename, dtype)
        # type safe creation
        data[n] = np.asarray([v for vec in vals for v in vec], dtype=dtype)
    # build data frame
    dfp = pd.DataFrame({c: pd.Series(v) for c, v in data.items()})
    dfp.loc[:, 'event'] = evns
    return dfp

# helper function to truncate color map (for a better view from the rainbow colormap)
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


# load root macros, input is an argument string
def load_root_macros(arg_macros):
    for path in arg_macros.split(','):
        path = path.strip()
        if os.path.exists(path):
            gROOT.Macro(path)
        else:
            print('\"{}\" does not exist, skip loading it.'.format(path))


# read mc particles from root file
def get_mcp_data(path, evnums=None, branch='MCParticles'):
    f = ROOT.TFile(path)
    events = f.events
    if evnums is None:
        evnums = np.arange(events.GetEntries())
    elif isinstance(evnums, int):
        evnums = [evnums]

    dbuf = np.zeros(shape=(5000*len(evnums), 6))
    idb = 0
    for iev in evnums:
        if iev >= events.GetEntries():
            print('Error: event {:d} is out of range (0 - {:d})'.format(iev, events.GetEntries() - 1))
            continue

        events.GetEntry(iev)
        # extract full mc particle data
        for part in getattr(events, branch):
            dbuf[idb] = (iev, part.momentum.x, part.momentum.y, part.momentum.z, part.PDG, part.simulatorStatus)
            idb += 1
    return pd.DataFrame(data=dbuf[:idb], columns=['event', 'px', 'py', 'pz', 'pid', 'status'])


# read mc particles from root file
def get_mcp_simple(path, evnums=None, branch='MCParticles'):
    f = ROOT.TFile(path)
    events = f.events
    if evnums is None:
        evnums = np.arange(events.GetEntries())
    elif isinstance(evnums, int):
        evnums = [evnums]

    dbuf = np.zeros(shape=(len(evnums), 6))
    idb = 0
    for iev in evnums:
        if iev >= events.GetEntries():
            print('Error: event {:d} is out of range (0 - {:d})'.format(iev, events.GetEntries() - 1))
            continue

        events.GetEntry(iev)
        # extract full mc particle data
        part = getattr(events, branch)[2]
        dbuf[idb] = (iev, part.momentum.x, part.momentum.y, part.momentum.z, part.PDG, part.simulatorStatus)
        idb += 1
    return pd.DataFrame(data=dbuf[:idb], columns=['event', 'px', 'py', 'pz', 'pid', 'status'])

#######################################
# read all mc particles from root file
#######################################
def get_all_mcp(path, evnums=None, branch='MCParticles'):
    f = ROOT.TFile(path)
    events = f.events
    if evnums is None:
        evnums = np.arange(events.GetEntries())
    elif isinstance(evnums, int):
        evnums = [evnums]

    dbuf = np.zeros(shape=(2000*len(evnums), 9))
    idb = 0
    for iev in evnums:
        if iev >= events.GetEntries():
            print('Error: event {:d} is out of range (0 - {:d})'.format(iev, events.GetEntries() - 1))
            continue

        events.GetEntry(iev)
        # extract mc particle data
        for ptl in getattr(events, branch):
            dbuf[idb] = (iev, ptl.momentum.x, ptl.momentum.y, ptl.momentum.z, ptl.PDG, ptl.simulatorStatus, ptl.endpoint.x, ptl.endpoint.y, ptl.endpoint.z)
            idb += 1

    return pd.DataFrame(data=dbuf[:idb], columns=['event', 'px', 'py', 'pz', 'pid', 'status', 'vex', 'vey', 'vez'])

# read hits data from root file
def get_hits_data(path, evnums=None, branch='RecoEcalBarrelImaginglHits'):
    f = ROOT.TFile(path)
    events = f.events
    if evnums is None:
        evnums = np.arange(events.GetEntries())
    elif isinstance(evnums, int):
        evnums = [evnums]

    dbuf = np.zeros(shape=(2000*len(evnums), 7))
    idb = 0
    for iev in evnums:
        if iev >= events.GetEntries():
            print('Error: event {:d} is out of range (0 - {:d})'.format(iev, events.GetEntries() - 1))
            continue

        events.GetEntry(iev)
        for ihit, hit in enumerate(getattr(events, branch)):
            dbuf[idb] = (iev, ihit, hit.layer, hit.position.x, hit.position.y,
                    hit.position.z, hit.energy*1000.)
            idb += 1

    return pd.DataFrame(data=dbuf[:idb], columns=['event', 'cluster', 'layer', 'x', 'y',
        'z', 'energy'])


# read layers data from root file
def get_layers_data(path, evnums=None, branch="EcalBarrelImagingClustersLayers"):
    f = ROOT.TFile(path)
    events = f.events
    if evnums is None:
        evnums = np.arange(events.GetEntries())
    elif isinstance(evnums, int):
        evnums = [evnums]

    dbuf = np.zeros(shape=(2000*len(evnums), 7))
    idb = 0
    for iev in evnums:
        if iev >= events.GetEntries():
            print('Error: event {:d} is out of range (0 - {:d})'.format(iev, events.GetEntries() - 1))
            continue

        events.GetEntry(iev)
        for icl, layer in enumerate(getattr(events, branch)):
            dbuf[idb] = (iev, icl, layer.layer,
                         layer.position.x, layer.position.y, layer.position.z,
                         layer.energy*1000.)
            idb += 1

    return pd.DataFrame(data=dbuf[:idb], columns=['event', 'cluster', 'layer', 'x', 'y',
        'z', 'energy'])


# read clusters data from root file
def get_clusters_data(path, evnums=None, branch='EcalBarrelImagingClustersReco'):
    f = ROOT.TFile(path)
    events = f.events
    if evnums is None:
        evnums = np.arange(events.GetEntries())
    elif isinstance(evnums, int):
        evnums = [evnums]

    dbuf = np.zeros(shape=(20*len(evnums), 6))
    idb = 0
    for iev in evnums:
        if iev >= events.GetEntries():
            print('Error: event {:d} is out of range (0 - {:d})'.format(iev, events.GetEntries() - 1))
            continue

        events.GetEntry(iev)
        for k, cl in enumerate(getattr(events, branch)):
            dbuf[idb] = (iev, k, cl.nhits, cl.energy*1000., cl.cl_theta, cl.cl_phi)
            idb += 1

    return pd.DataFrame(data=dbuf[:idb], columns=['event', 'cluster', 'nhits', 'energy', 'cl_theta', 'cl_phi'])


def compact_constants(path, names):
    if not os.path.exists(path):
        print('Cannot find compact file \"{}\".'.format(path))
        return []
    if not names:
        return []
    kernel = DDG4.Kernel()
    description = kernel.detectorDescription()
    kernel.loadGeometry("file:{}".format(path))
    try:
        vals = [description.constantAsDouble(n) for n in names]
    except:
        print('Fail to extract values from {}, return empty.'.format(names))
        vals = []
    kernel.terminate()
    return vals
