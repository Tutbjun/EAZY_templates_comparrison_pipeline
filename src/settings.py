import yaml
import os
##################### IMPORT SETTINGS #####################
#from settings.yaml
with open(os.path.join(os.path.dirname(__file__),"settings.yaml"), 'r') as stream:
    settings = yaml.safe_load(stream)
for key in settings['CONFIG']:
    globals()[key] = settings['CONFIG'][key]
for key in settings['MOSAICZ']:
    globals()[key] = settings['MOSAICZ'][key]
del os, yaml, stream, key

#known keys
"""zmin = zmin
zMax = zMax
mosTiling = mosTiling
mm = eval(mm)
figwidth = eval(figwidth)
DPI = DPI
ro = ro
zCharacteristic = zCharacteristic
fontsize = fontsize"""
mm = eval(mm)
figwidth = eval(figwidth)