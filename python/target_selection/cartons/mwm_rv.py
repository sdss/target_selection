#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-06-10
# @Filename: mwm_rv.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
import peewee

# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             CatalogToTIC_v8,
                                             TIC_v8,
                                             SDSS_APOGEE_AllStarMerge_r13,
                                             TwoMassPSC)
# from . import BaseCarton
from target_selection.cartons import BaseCarton


class MWM_YSO_RV_Long_RM_Carton(BaseCarton):
    """ 2.2.1. Long Baseline (Legacy Targets)

Shorthand name: mwm_rv_long (Not an actual target class.
Defined as a parent sample to select from for 2.2.1.x sections below)

Simplified Description of selection criteria:
APOGEE-1&-2 "main survey" targets with at least 3 previous epochs
brighter than H of 12.8 that have a Gaia parallax.

Wiki page: Binaries and Substellar Companions

Additional source catalogs needed: Just need sdss_apogeeallstarmerge_r13

Additional cross-matching needed: (None)

Return columns: apogee_id, nvisits, ra, dec, pmra, pmdec, h, baseline, fields

cadence options for these targets (list all options,
even though no single target will receive more than one):
(See entries for this in 2.2.1.x below)

Pseudo SQL (optional): SELECT apogee_id, nvisits, ra, dec,
pmra, pmdec, h, baseline, fields FROM sdss_apogeeallstarmerge_r13
WHERE h<12.8 AND nvisits>=3 AND dist_src == 'gaia' AND
[targflags contains 'APOGEE_SHORT' OR 'APOGEE_INTERMEDIATE' OR
'APOGEE_LONG' OR '*APOGEE2_*BIN_’] AS mwm_rv_long

Implementation:

Non-SQL implementation:

lead contact: Nicholas Troup

Target selection final for v0?: No
2.2.1.1. Long Baseline Targets for RM Plates

Shorthand name: mwm_rv_long_rm

Simplified Description of selection criteria:
Select from long-baseline targets (above) within
the 3 BHM-RM fields observable from APO
(sexagesimal coordinates given in Pseudo SQL below)

Wiki page: Binaries and Substellar Companions

Additional source catalogs needed: Select from mwm_rv_long (2.2.1)

Additional cross-matching needed: (None)

Return columns: apogee_id, nvisits, ra, dec, pmra, pmdec, h,
baseline, fields (Same as 2.2.1 above)

cadence options for these targets (list all options,
even though no single target will receive more than one):
mwm_rv_<nn>x8 (where <nn> = number of visits to RM fields)

Pseudo SQL (optional): SELECT apogee_id, nvisits, ra, dec,
pmra, pmdec, h, baseline, fields FROM mwm_rv_long
WHERE [(ra,dec) < 3 degrees on sky from
(14:14:49 +53:05:00) OR (10:00:00 +02:12:00) OR (02:23:30 -04:15:00)]

Implementation:

Non-SQL implementation:

lead contact: Nicholas Troup

Target selection final for v0?: No



We associate each APOGEE target with a 2MASS source since the
vast majority of APOGEE targets were selected from the 2MASS PSC.
APOGEE_ID is essentially the same as
the “designation” column of the 2MASS PSC;
We just have to strip off the “2M” in the APOGEE_ID.

For example:
sdss5db=# select designation from  catalogdb.twomass_psc limit 2;
   designation
------------------
 12034281-5740002
 12032366-5738380
(2 rows)

sdss5db=# select apogee_id from  catalogdb.sdss_apogeeallstarmerge_r13 limit 2;
     apogee_id
--------------------
 2M14044120-1550575
 2M14033676-1554164
(2 rows)

sdss5db=# select trim( leading '2M' from apogee_id) from
 catalogdb.sdss_apogeeallstarmerge_r13 limit 2;
      ltrim
------------------
 14044120-1550575
 14033676-1554164
(2 rows)

    """
    name = 'mwm_rv_long_rm'
    category = 'science'
    cadence = None
    program = 'program_mwm_rv_long_rm'

# peewee Model name ---> postgres table name
# SDSS_APOGEE_AllStarMerge_r13(CatalogdbModel)--->'sdss_apogeeallstarmerge_r13'

    def build_query(self, version_id, query_region=None):
        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid)
                 .join(TIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .join(SDSS_APOGEE_AllStarMerge_r13,
                       on=(TwoMassPSC.designation ==
                           peewee.fn.ltrim(SDSS_APOGEE_AllStarMerge_r13.apogee_id,
                                           '2M')))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        SDSS_APOGEE_AllStarMerge_r13.h < 12.8,
                        SDSS_APOGEE_AllStarMerge_r13.nvisits > 3,
                        SDSS_APOGEE_AllStarMerge_r13.dist_src == 'gaia',
                        (SDSS_APOGEE_AllStarMerge_r13.targflags % '%APOGEE_SHORT%') |
                        (SDSS_APOGEE_AllStarMerge_r13.targflags % '%APOGEE_INTERMEDIATE%') |
                        (SDSS_APOGEE_AllStarMerge_r13.targflags % '%APOGEE_LONG%') |
                        (SDSS_APOGEE_AllStarMerge_r13.targflags % '%APOGEE2_%BIN_%')))
        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                           Catalog.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))
        return query
