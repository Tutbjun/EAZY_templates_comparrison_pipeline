#remove all subdirectories with 'Kirkpatrick' in the name

import os
import shutil

removeType = 'Kirkpatrick' #string to be removed from subdirectory name
#path to the directory containing the subdirectories to be removed
for resamplefolder in os.listdir(os.path.dirname(__file__)):
    if resamplefolder == 'remove_set.py':
        continue
    subdirs = os.listdir(os.path.join(os.path.dirname(__file__), resamplefolder, 'eazy-output'))
    for subdir in subdirs:
        if removeType in subdir:
            shutil.rmtree(os.path.join(os.path.dirname(__file__), resamplefolder, 'eazy-output', subdir))