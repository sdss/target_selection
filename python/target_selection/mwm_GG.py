import sdssdb
import astropy
from os import path
from astropy.table import Table
from sdssdb.peewee.sdss5db.catalogdb import database, GaiaDR2Source, TwoMassPsc, GaiaDR2TmassBestNeighbour


#  TODO 
#  Split query for performance
#  Perform 2MASS phot quality checks in SQL  (optional, otherwise use python)
#  Output needs to eventually be input for targetDB
#  Use subqueries and perform outer join 
#  Make sure all relevant output columns are there


# 'GG_targetdb_eqivalent.fits'

names = ['Hlt7']
outfiles = ['GG_targetdb_eqivalent_{0}.fits'.format(names[0])]
outpaths = [path.join(path.expanduser('~'), outfile) for outfile in outfiles]
print(outpaths)

database.autorollback = True
database.set_profile('sdssadmin')

gaia_source = GaiaDR2Source.alias()
gaia_tmass = GaiaDR2TmassBestNeighbour.alias()
tmass = TwoMassPsc.alias()


tmass_mag_ranges = [(tmass.h_m < 7.0)]

query1 = tmass.select(gaia_source.source_id,
                    tmass.designation,
                    gaia_source.ra,
                    gaia_source.dec,
                    gaia_source.phot_g_mean_mag,
                    tmass.h_m, tmass.ph_qual,
                    gaia_tmass.angular_distance).join(gaia_tmass).join(gaia_source).where(tmass_mag_ranges[0] & (gaia_source.phot_g_mean_mag - tmass.h_m > 3.5))


targ_table = Table(list(query1.dicts()))
targ_table.write(outpaths[0])
database.close()
