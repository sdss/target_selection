#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2023-04-24
# @Filename: mwm_monitor_apogee.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (Catalog, CatalogToGaia_DR3,
                                             CatalogToTwoMassPSC, Gaia_DR3,
                                             SDSS_DR17_APOGEE_Allstarmerge,
                                             TwoMassPSC)

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py
#
# peewee model : postgres table name
# Catalog : catalogdb.catalog
# Gaia_DR3 : catalogdb.gaia_dr3_source
# SDSS_DR17_APOGEE_Allstarmerge : catalogdb.sdss_dr17_apogee_allstarmerge
# TwoMassPSC : catalogdb.twomass_psc
#
#
# We associate each APOGEE target with a 2MASS source since the
# vast majority of APOGEE targets were selected from the 2MASS PSC.
# APOGEE_ID is the same as
# the “designation” column of the 2MASS PSC;
#
# For example:
# sdss5db=# select designation from  catalogdb.twomass_psc limit 2;
#    designation
# ------------------
#  12034281-5740002
#  12032366-5738380
# (2 rows)
#
# For sdss_dr17_apogee_allstarmerge we
# have to strip off the “2M” in the APOGEE_ID.
# However, some apogee_id values do not have 2M in the front.
# (old: for sdss_apogeeallstarmerge_r13 we also had to strip off the
#  “2M” in the APOGEE_ID.)
# For example:
# sdss5db=# select designation from  catalogdb.twomass_psc limit 2;
#    designation
# ------------------
#  12034281-5740002
#  12032366-5738380
# (2 rows)
#
# sdss5db=# select apogee_id from catalogdb.SDSS_DR17_APOGEE_Allstarmerge
# where apogee_id like '2M%';
#      apogee_id
# ---------------------
# 2M00000002+7417074
# 2M00000019-1924498
# etc.
# sdss5db=# select apogee_id from catalogdb.SDSS_DR17_APOGEE_Allstarmerge
# where apogee_id not like '2M%';
#    apogee_id
# --------------------
# 19140272-1554055
# 19155129-1617591
#
#

class MWM_monitor_apogee_n188_long_Carton(BaseCarton):
    """ 5.1.25. mwm_monitor_apogee_*

Shorthand name: mwm_monitor_apogee

Existing carton code: None

Simplified Description of selection criteria:
These are targets selected from globular and open clusters that
have been observed many times over many years with the apogee instrument.

mwm_monitor_apogee_n188_long:
SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%N188%'
AND baseline >= 3000 AND nvisits >= 12  high priority

mwm_monitor_apogee_n188_short:
SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%N188%'
AND baseline < 3000 AND baseline >= 1800 AND nvisits >= 12  lower priority

mwm_monitor_apogee_m67_long:
SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%M67%'
AND baseline >= 2000 AND nvisits >= 12  high priority

mwm_monitor_apogee_m67_short:
SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%M67%'
AND baseline < 2000 AND baseline >= 1800 AND nvisits >= 12  lower priority

mwm_monitor_apogee_m15_long:
SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%M15%'
AND baseline >= 1300 AND nvisits >= 12  high priority

mwm_monitor_apogee_m15_short:
SELECT apogee_id,ra,dec,nvisits,baseline,j,j_err,h,h_err,k,k_err,fields
FROM catalogdb.sdss_dr17_apogee_allstarmerge WHERE fields LIKE '%M15%'
AND baseline < 1300 AND baseline >= 900 AND nvisits >= 12  lower priority

Gaia DR2 parameters to be converted to Gaia DR3: Nothing to be converted.

Return columns: RA, DEC, proper motion, hmag, 2mass name.

Metadata:

Priority:
High priority=1300
Lower Priority =1302

Cadence: bright_1x4

Instrument: APOGEE

can_offset=True

Lead contact: Nathan De Lee
    """

    name = 'mwm_monitor_apogee_n188_long'
    category = 'science'
    instrument = 'APOGEE'
    cadence = 'bright_1x4'
    program = 'mwm_legacy'
    mapper = 'MWM'
    priority = 1300
    can_offset = True

    # In the below query, we use replace() instead of ltrim() since
    # ltrim('2M20', '2M') will also trim the second 2.

    def build_query(self, version_id, query_region=None):

        query = (Catalog
                 .select(CatalogToGaia_DR3.catalogid,
                         Gaia_DR3.source_id,
                         Gaia_DR3.ra.alias('gaia_dr3_ra'),
                         Gaia_DR3.dec.alias('gaia_dr3_dec'),
                         Gaia_DR3.phot_g_mean_mag,
                         Gaia_DR3.phot_bp_mean_mag,
                         Gaia_DR3.phot_rp_mean_mag,
                         TwoMassPSC.pts_key,
                         TwoMassPSC.designation)
                 .join(CatalogToGaia_DR3,
                       on=(Catalog.catalogid == CatalogToGaia_DR3.catalogid))
                 .join(Gaia_DR3,
                       on=(CatalogToGaia_DR3.target_id == Gaia_DR3.source_id))
                 .switch(CatalogToGaia_DR3)
                 .join(CatalogToTwoMassPSC,
                       on=(CatalogToGaia_DR3.catalogid == CatalogToTwoMassPSC.catalogid))
                 .join(TwoMassPSC,
                       on=(CatalogToTwoMassPSC.target_id == TwoMassPSC.pts_key))
                 .join(SDSS_DR17_APOGEE_Allstarmerge,
                       on=(TwoMassPSC.designation ==
                           peewee.fn.replace(SDSS_DR17_APOGEE_Allstarmerge.apogee_id, '2M', '')))
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        Gaia_DR3.phot_g_mean_mag.between(8, 18)))

        # Gaia_DR3 peewee model class corresponds to
        # table catalogdb.gaia_dr3_source.

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
