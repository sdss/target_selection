# encoding: utf-8
#
# @Author: Tom Dwelly
# @Date: Jan 2020
# @Filename: bhm_aqmes.py
# @License: BSD 3-Clause
# @Copyright: Tom Dwelly
# Content:
#    Target selection for the AQMES components of BHM
#    AQMES targets are taken from SDSS QSO catalogue (DR16Q)

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import pkg_resources

import numpy as np
import numpy.lib.recfunctions as rf
from datetime import datetime
import pymangle as mangle

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

from .print_func import *
from .catalogdb_wrapper import retrieve_targets_from_catalogdb

__version__ = pkg_resources.require("sdss-target_selection")[0].version

def write_targets_to_fits(targets=None, filename=None):
    assert isinstance(targets, np.recarray)
    assert isinstance(filename, str)

    #convert to binaryTableHDU
    hdu = fits.BinTableHDU.from_columns(targets)

    # add some keywords
    hdu.header.append(fits.Card('CREATOR', __name__, 'Software name'))
    hdu.header.append(fits.Card('VERSION', __version__, 'Software version'))
    hdu.header.append(fits.Card('DATE', datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")[:-4]))

    # write to disk
    hdu.writeto(filename, overwrite=True)


def get_aqmes_initial_candidates(db_table='DR14QV44',
                                 peewee_selection_expr=["db.psfmag[3] > 16.0",
                                                        "db.psfmag[3] < 19.1"],
                                 col_dic={'col_id':'pk',
                                          'col_ra':'ra',
                                          'col_dec':'dec',
                                          'col_mag':'psfmag',}
                                 ):
    '''
    Carry out initial target selection for the BHM_AQMES components
    Using the DR1xQ catalogue as the basis
    '''
    assert db_table is not None
    assert col_dic is not None

    db_cols = list(col_dic.values())

    cands = retrieve_targets_from_catalogdb(tablename=db_table,
                                            constraints=peewee_selection_expr,
                                            fmt='recarray',
                                            columns=db_cols)
    print_comment("Num QSOs meeting PeeWee selection ({}): {}".format(
        peewee_selection_expr,
        len(cands)))

    return cands, col_dic




def select_aqmes_wide_targets(cands=None,
                              col_dic=None,
                              value_boost=10.0,
                              dec_thresh=+15.0):

    assert cands is not None
    mask_filename = pkg_resources.resource_filename(__name__,
                                                    'masks/AQMES_wboss_imaging_NGC_proc.ply')

    # filter on mask position
    try:
        mask = mangle.Mangle(mask_filename)
    except:
        raise Exception(f"Failed to find/open maskfile: {mask_filename}")

    in_mask = mask.contains(cands[col_dic['col_ra']], cands[col_dic['col_dec']])
    ntargets = np.sum(in_mask)

    print_comment(f"Num QSOs in spatial mask ({mask_filename}): {ntargets}")

    coldata_ra      = np.extract(in_mask, cands[col_dic['col_ra']]).astype(np.float64)
    coldata_dec     = np.extract(in_mask, cands[col_dic['col_dec']]).astype(np.float64)
    # split the sample by declination
    in_3x4 = np.where(coldata_dec>dec_thresh, 1, 0)

    print_comment (f"Num QSOs in AQMES_WIDE_3x4 region: {(in_3x4 == 1).sum()}")
    print_comment (f"Num QSOs in AQMES_WIDE_2x4 region: {(in_3x4 == 0).sum()}")

    coldata_cadence = np.where(in_3x4, 'bhm_aqmes_wide_3x4', 'bhm_aqmes_wide_2x4').astype('U24')
    coldata_type    = np.where(in_3x4, 'BHM_AQMES_WIDE3', 'BHM_AQMES_WIDE2').astype('U24')
    coldata_value   = np.where(in_3x4, 3.0*value_boost, 2.0*value_boost).astype(np.float32)
    coldata_priority = np.full(ntargets, 1.0, dtype=np.float32)

    targets = np.core.records.fromarrays(
        [coldata_ra,
         coldata_dec,
         coldata_type,
         coldata_cadence,
         coldata_value,
         coldata_priority,
        ],
         names=['ra',
               'dec',
               'type',
               'cadence',
               'value',
               'priority',
         ],
    )


    return targets



def select_aqmes_medium_targets(cands=None,
                                col_dic=None,
                                value_boost=10.0,
                                ntargets_req=3000,
                                dec_limits={'lower':-5.0,
                                            'upper':+45.0},):

    assert cands is not None
    assert ntargets_req > 0
    col_field_id = 'fieldid'
    col_field_ra = 'racen'
    col_field_dec = 'deccen'
    col_field_type = 'type'
    col_field_radius = 'radius'
    keep_cols=[col_field_id, col_field_ra, col_field_dec, col_field_type, col_field_radius]

    fields_filename = pkg_resources.resource_filename(__name__,
                                                      'masks/rsFields_boss_survey_sgc.fits')
    try:
        maskhdul = fits.open(fields_filename)
    except:
        raise Exception(f"Failed to find/open fieldlist file: {fields_filename}")

    d = maskhdul[1].data
    drop_cols = [c for c in d.columns.names if c not in keep_cols]
    if len(drop_cols) > 0:
        d = rf.drop_fields(d, drop_cols, usemask=False, asrecarray=True)

    # mask the fields based on Dec limits
    m = np.where((d[col_field_dec] >= dec_limits['lower']) & (d[col_field_dec] <= dec_limits['upper']), True, False)
    fields = np.extract(m, d)
    maskhdul.close()

    print_comment(f"Num fields meeting selection criteria ({dec_limits['lower']}<=Dec<={dec_limits['upper']}): {len(fields)}")

    # add some new columns to the recarry
    fields = rf.append_fields(fields,
                              ['nqso',
                               'nqsou',
                               'bhmcadence',
                              ],
                              [np.zeros(len(fields), dtype=np.int32),
                               np.zeros(len(fields), dtype=np.int32),
                               np.repeat('bhm_aqmes_med_12x4', len(fields)).astype(np.dtype('U24')),
                              ],
                              usemask=False,
                              asrecarray=True)


    # now choose the optimal fields
    last_nth_field = -1
    iter = 1  ## two or three iterations is generally enough for this to converge
    num_qsou = -1

    ### assume that the first field has the correct radius
    seplimit = u.Quantity(fields[col_field_radius][0], unit='deg')

    while iter <= 10 :
        # associate QSOs with fields
        # http://docs.astropy.org/en/stable/api/astropy.coordinates.search_around_sky.html
        coords_f = SkyCoord(fields[col_field_ra], fields[col_field_dec], frame="icrs", unit="deg")
        coords_q = SkyCoord(cands[col_dic['col_ra']], cands[col_dic['col_dec']], frame="icrs", unit="deg")
        idx_q, idx_f, d2d, d3d  = search_around_sky(coords_q, coords_f, seplimit)

        # count total qsos per field
        for f in np.arange(len(fields)) :
            fields['nqso'][f]=(idx_f == f).sum()

        # count unique qsos per field
        qsou = np.unique(idx_q)
        fields['nqsou'] = 0
        for q in qsou :
            fs = idx_f[(idx_q == q)]
            seps = d2d[(idx_q == q)]
            imin = np.argmin(seps)
            fields['nqsou'][fs[imin]] += 1

        # sort according to number of unique qsos
        fields = np.sort(fields, kind='quicksort', order='nqsou')[::-1]

        num_qso = 0
        num_qsou = 0
        nth_field = 0
        for f in fields :
            num_qso += f['nqso']
            num_qsou += f['nqsou']
            nth_field += 1
            if num_qsou >= ntargets_req :
                break

        # test for convergence
        if nth_field == last_nth_field :
            print_comment ("Converged")
            break

        fields = fields[:nth_field]
        print_comment("After iter {}: The top {} fields contain {} unique QSOs ({} QSO-field matches)".format(
            iter, nth_field, num_qsou, num_qso))
        iter += 1
        last_nth_field = nth_field


    ## Now prepare the target catalogue data
    ## make a mask of input targets that fall in the fields
    mask = np.zeros(len(cands), dtype=int)
    mask[qsou] = 1

    coldata_ra      = np.extract(mask==1, cands[col_dic['col_ra']]).astype(np.float64)
    coldata_dec     = np.extract(mask==1, cands[col_dic['col_dec']]).astype(np.float64)
    ntargets = len(coldata_ra)
    coldata_cadence  = np.repeat('bhm_aqmes_med_12x4', ntargets).astype('U24')
    coldata_type     = np.repeat('BHM_AQMES_MEDIUM', ntargets).astype('U24')
    coldata_value    = np.repeat(12.0*value_boost, ntargets).astype(np.float32)
    coldata_priority = np.repeat(1, ntargets).astype(np.float32)

    targets = np.core.records.fromarrays([coldata_ra,
                                          coldata_dec,
                                          coldata_type,
                                          coldata_cadence,
                                          coldata_value,
                                          coldata_priority],
                                         names=['ra',
                                                'dec',
                                                'type',
                                                'cadence',
                                                'value',
                                                'priority',
                                         ],)

    return targets, fields
