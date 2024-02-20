These data files contain mock integrated spectral energy distributions
(SEDs) for the simulations of isolated galaxies and galaxy mergers
presented in Lanz et al. (2013), along with auxiliary information regarding
the simulations. The file SIGS_sims_I_DISM_SEDs.sav contains the data for
the default ISM (DISM) model, whereas SIGS_sims_I_AISM_SEDs.sav contains
the data for the alternate ISM (AISM) model.

Contact Chris Hayward (current contact information at
www.cfa.harvard.edu/~chayward) if you have questions or problems.

The files are IDL save files, which can be read in IDL with the command

restore,'SIGS_sims_I_DISM_SEDs.sav'

which directly loads the data arrays into the IDL session. Type 'help' to
see the arrays that have been loaded.

The save files can also be read in Python with the commands

import scipy.io
idl_dict=scipy.io.readsav('SIGS_sims_I_DISM_SEDs.sav')

The readsav command loads the data into a case-insensitive dictionary with
item, attribute, and call access to variables. Type 'idl_dict.keys()' to see
the keys for the items that have been loaded. See
http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.readsav.html
for further details.

NOTE: Missing data (primarily because only a snapshot of the Gadget-3
snapshots were processed with Sunrise and because the simulations have
different numbers of snapshots) are indicated by NaN's.

The data contained in the files are the following:

lambda: array of wavelength values for the SEDs in meters.

l_lambda: integrated luminosity density values in units of W/m. The indices
	are [simulation, time snapshot, camera, wavelength].

runname: This array gives the simulation identifier for a given simulation
	index. The simulations 'M0', 'M1', 'M2', and 'M3' are isolated disk
	galaxies. See Table 2 of the publication for details. The other
	simulations are mergers; for these, the identifiers have the format
	<progenitor 1><progenitor 2><orbit>. Thus, 'M2M1e' is a merger of the
	'M2' and 'M1' disk galaxies on the 'e' orbit (which is the only orbit
	used in this work).

The following arrays have the dimensions [simulation, time snapshot]:

d_bh: separation of the central black holes of the galaxies (only applicable
	for the merger simulations) in kpc. This array can be used to define the
	time after coalescence that is used in the paper.

mdust: total dust mass in solar masses.

mgas: total gas mass in solar masses.

mstar: total stellar mass in solar masses.

sfr: total instantaneous star formation rate in solar masses per year.

time: time since the start of the simulation in Gyr.

Example (IDL): To plot the SED for the 'M3M3e' merger at snapshot 80 viewed
from camera 3, do

runidx = where(runname eq 'M3M3e')
plot,lamda,lambda*L_lambda[runidx,80,3,*],/xlog,/ylog

Then, to get information regarding this snapshot, do, e.g.,

print,time[runidx,80]
print,sfr[runidx,80]
