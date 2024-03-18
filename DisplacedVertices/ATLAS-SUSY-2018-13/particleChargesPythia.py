#!/usr/bin/env python3

from xml.etree import ElementTree as ET
import os 

# Load Pythia Particle Data and create a dictionary
pdataFile = os.path.join(os.path.dirname(__file__), "ParticleData.xml")
xml = ET.parse(pdataFile)
root_element = xml.getroot()

particleCharges = {}
colorCharges = {}
for child in root_element:
    particleCharges[int(child.attrib["id"])] = eval(child.attrib["chargeType"])
    colorCharges[int(child.attrib["id"])] = eval(child.attrib["colType"])
