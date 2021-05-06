#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-17
# @Filename: bhm_spiders_clusters.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# isort: skip_file

import peewee
from peewee import JOIN
from peewee import fn


from target_selection.cartons.base import BaseCarton
from target_selection.mag_flux import AB2nMgy, AB2Jy

# general catalogdb imports
from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    EROSITASupersetClusters,
)

# imports of existing spectro catalogues
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogToSDSS_DR16_SpecObj,
    SDSS_DR16_SpecObj,
    CatalogToBHM_eFEDS_Veto,
    BHM_eFEDS_Veto,
    SDSSV_BOSS_SPALL,
    SDSSV_BOSS_Conflist,
    SDSSV_Plateholes,
    SDSSV_Plateholes_Meta,
)

# additional imports required by bhm_spiders_clusters_lsdr8
from sdssdb.peewee.sdss5db.catalogdb import (
    CatalogToLegacy_Survey_DR8,
    Legacy_Survey_DR8,
)

# additional imports required by bhm_spiders_clusters_ps1dr2
from sdssdb.peewee.sdss5db.catalogdb import (
    Panstarrs1,
    CatalogToPanstarrs1,    # only exists after v0.5 cross-match
)

# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms

# ############################################
# ############################################
# ############################################
# ############################################
# # This file provides the following BHM cartons in v0.5:
#
#   bhm_spiders_clusters_efeds_stragglers
#   bhm_spiders_clusters_lsdr8
#   bhm_spiders_clusters_ps1dr2
#
# ############################################
# ############################################
# ############################################
# ############################################

# Notes on how many targets to expect:
# => SELECT ero_version,xmatch_method,xmatch_version,opt_cat,count(*)
#      FROM erosita_superset_clusters GROUP BY ero_version,xmatch_method,xmatch_version,opt_cat;
#        ero_version        |      xmatch_method |      xmatch_version      |   opt_cat    | count
# --------------------------+--------------------+--------------------------+--------------+--------
#  eFEDS_c001_V18C_legacy_d | EROMAPPER_LS_DR8   | grzw1_v0.3_2020-12-04    | lsdr8        |  73768
#  em01_c946_201008_poscorr | EROMAPPER_LS_DR8   | grzw1_v0.3_2020-12-04    | lsdr8        | 160465
#  em01_c946_201008_poscorr | EROMAPPER_PS1_DR2  | eromapper_2020-10-23     | ps1dr2       | 140808

#
#
# END PREAMBLE
# ##################################################################################


class BhmSpidersClustersLsdr8Carton(BaseCarton):

    name = 'bhm_spiders_clusters_lsdr8'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        ls = Legacy_Survey_DR8.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()
        s2020 = BHM_eFEDS_Veto.alias()
        sV = SDSSV_BOSS_SPALL.alias()

        xx = EROSITASupersetClusters.alias()
        x = (
            xx
            .select(
                fn.rank().over(partition_by=[xx.ero_detuid],
                               order_by=[xx.xmatch_metric.desc()]).alias('x_rank'),
                xx.ero_detuid.alias('ero_detuid'),
                xx.ls_id.alias('ls_id'),
                xx.target_has_spec.alias('target_has_spec'),
            )
            .where(
                (xx.ero_version == self.parameters['ero_version']),
                (xx.xmatch_method == self.parameters['xmatch_method']),
                (xx.xmatch_version == self.parameters['xmatch_version']),
                (xx.opt_cat == self.parameters['opt_cat']),
                (xx.xmatch_metric > self.parameters['xmatch_metric_min']),
                (xx.ero_det_like > self.parameters['det_like_min']),
            )
            .alias('x')
        )

        instrument = peewee.Value(self.instrument)
        inertial = peewee.Value(self.inertial).cast('bool')

        fibertotflux_r_max = AB2nMgy(self.parameters['fibertotmag_r_min'])
        fibertotflux_r_min = AB2nMgy(self.parameters['fibertotmag_r_max'])
        fibertotflux_z_max = AB2nMgy(self.parameters['fibertotmag_z_min'])
        fibertotflux_z_min = AB2nMgy(self.parameters['fibertotmag_z_max'])

        fibertotflux_r_min_for_cadence1 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence1'])
        fibertotflux_z_min_for_cadence1 = AB2nMgy(self.parameters['fibertotmag_z_for_cadence1'])
        fibertotflux_r_min_for_cadence2 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence2'])
        gaia_g_max_for_cadence1 = self.parameters['gaia_g_max_for_cadence1']
        gaia_rp_max_for_cadence1 = self.parameters['gaia_rp_max_for_cadence1']

        flux30 = AB2nMgy(30.00)
        # match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0

        # #########################################################################
        # prepare the spectroscopy catalogues

        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        # SDSS DR16
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        ss16 = SDSS_DR16_SpecObj.alias()
        s16 = (
            ss16.select(
                ss16.specobjid.alias('specobjid'),
            )
            .where(
                ss16.snmedian >= spec_sn_thresh,
                ss16.zwarning == 0,
                ss16.zerr <= spec_z_err_thresh,
                ss16.zerr > 0.0,
                ss16.scienceprimary > 0,
            )
            .alias('s16')
        )

        # SDSS-IV/eFEDS March2020
        c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
        ss2020 = BHM_eFEDS_Veto.alias()
        s2020 = (
            ss2020.select(
                ss2020.pk.alias('pk'),
            )
            .where(
                ss2020.sn_median_all >= spec_sn_thresh,
                ss2020.zwarning == 0,
                ss2020.z_err <= spec_z_err_thresh,
                ss2020.z_err > 0.0,
            )
            .alias('s2020')
        )

        # SDSS-V spAll
        ssV = SDSSV_BOSS_SPALL.alias()
        sV = (
            ssV.select(
                ssV.specobjid.alias('specobjid'),
                ssV.plug_ra.alias('plug_ra'),
                ssV.plug_dec.alias('plug_dec'),
            )
            .where(
                ssV.sn_median_all >= spec_sn_thresh,
                ssV.zwarning == 0,
                ssV.z_err <= spec_z_err_thresh,
                ssV.z_err > 0.0,
                ssV.specprimary > 0,
            )
            .alias('sV')
        )

        # SDSS-V plateholes - only consider plateholes that
        # were drilled+shipped but that were not yet observed
        ssph = SDSSV_Plateholes.alias()
        ssphm = SDSSV_Plateholes_Meta.alias()
        ssconf = SDSSV_BOSS_Conflist.alias()
        sph = (
            ssph.select(
                ssph.pkey.alias('pkey'),
                ssph.target_ra.alias('target_ra'),
                ssph.target_dec.alias('target_dec'),
            )
            .join(
                ssphm,
                on=(ssph.yanny_uid == ssphm.yanny_uid)
            )
            .join(
                ssconf, JOIN.LEFT_OUTER,
                on=(ssphm.plateid == ssconf.plate)
            )
            .where(
                (ssph.holetype == 'BOSS_SHARED'),
                (ssph.sourcetype == 'SCI') | (ssph.sourcetype == 'STA'),
                ssphm.isvalid > 0,
                ssconf.plate.is_null(),
            )
            .alias('sph')
        )

        # priority is determined by target rank within cluster
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:

        priority = peewee.Case(
            None,
            (
                (
                    x.c.x_rank == 1,
                    self.parameters['priority_floor_bcg']
                ),
                (
                    x.c.x_rank > 1,
                    self.parameters['priority_floor_member'] +
                    fn.least(self.parameters['priority_levels'] - 2,
                             x.c.x_rank - 2)
                ),
            ),
            None)

        value = peewee.Case(
            None,
            (
                (x.c.x_rank == 1, self.parameters['value_bcg']),
                (x.c.x_rank > 1, self.parameters['value_member']),
            ),
            None)

        # choose cadence based on fiber magnitude in r-band
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence4 = 'unknown_cadence'  # catch failures
        cadence = peewee.Case(
            None,
            (
                (
                    ((ls.fibertotflux_r > fibertotflux_r_min_for_cadence1) |
                     (ls.fibertotflux_z > fibertotflux_z_min_for_cadence1) |
                     (ls.gaia_phot_g_mean_mag.between(0.1, gaia_g_max_for_cadence1)) |
                     (ls.gaia_phot_rp_mean_mag.between(0.1, gaia_rp_max_for_cadence1))),
                    cadence1),
                (ls.fibertotflux_r > fibertotflux_r_min_for_cadence2, cadence2),
                (ls.fibertotflux_r <= fibertotflux_r_min_for_cadence2, cadence3),
            ),
            cadence4)

        # We want to switch between psfmags and fibertotmags depending on
        # ls.type parameter (PSF or extended)
        # For 'PSF' targets, we use psfmags, but for extended sources use fiber2mags
        opt_prov = peewee.Case(
            ls.type,
            (('PSF', 'ls_psfmag'),),
            'ls_fibertotmag')

        magnitude_g = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_g))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_g))).cast('float'))

        magnitude_r = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_r))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_r))).cast('float'))

        magnitude_z = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_z))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_z))).cast('float'))

        magnitude_i = peewee.Case(
            ls.type,
            (('PSF',
              (22.5 - 2.5 * fn.log10(
                  fn.greatest(flux30, 0.5 * (ls.flux_r + ls.flux_z)))).cast('float')),),
            (22.5 - 2.5 * fn.log10(
                fn.greatest(flux30, 0.5 * (ls.fibertotflux_r + ls.fibertotflux_z)))).cast('float'))

        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        query = (
            c.select(
                c.catalogid.alias('catalogid'),
                ls.ls_id.alias('ls_id'),
                x.c.ero_detuid.cast('text').alias('ero_detuid'),
                c.ra.alias('ra'),
                c.dec.alias('dec'),
                priority.alias('priority'),
                value.cast('float').alias('value'),
                cadence.alias('cadence'),
                instrument.alias('instrument'),
                magnitude_g.alias('g'),
                magnitude_r.alias('r'),
                magnitude_i.alias('i'),
                magnitude_z.alias('z'),
                opt_prov.alias('optical_prov'),
                inertial.alias('inertial'),
            )
            .join(c2ls)
            .join(ls)
            .join(x, on=(ls.ls_id == x.c.ls_id))
            # start joining the spectroscopy
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid) &
                    (c2s16.version_id == version_id)
                )
            )
            .switch(c)
            .join(c2s2020, JOIN.LEFT_OUTER)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    (c2s2020.target_id == s2020.c.pk) &
                    (c2s2020.version_id == version_id)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            .join(
                sph, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            # finished joining the spectroscopy
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                c2ls.best >> True
            )
            .where(
                s16.c.specobjid.is_null(True),  # all of these must be satisfied
                s2020.c.pk.is_null(True),
                sV.c.specobjid.is_null(True),
                sph.c.pkey.is_null(True),
            )
            .where(
                (
                    (ls.fibertotflux_r.between(fibertotflux_r_min, fibertotflux_r_max)) |
                    (ls.fibertotflux_z.between(fibertotflux_z_min, fibertotflux_z_max))
                ),
                (x.c.target_has_spec == 0),
                # gaia safety checks to avoid bad ls photometry
                ~(ls.gaia_phot_g_mean_mag.between(0.1, self.parameters['gaia_g_mag_limit'])),
                ~(ls.gaia_phot_rp_mean_mag.between(0.1, self.parameters['gaia_rp_mag_limit'])),
            )
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query

#
# END BhmSpidersClustersLsdr8Carton
# ##################################################################################


class BhmSpidersClustersEfedsStragglersCarton(BaseCarton):

    name = 'bhm_spiders_clusters_efeds_stragglers'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        ls = Legacy_Survey_DR8.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()
        # s2020 = BHM_eFEDS_Veto.alias()
        # sV = SDSSV_BOSS_SPALL.alias()
        # ph = SDSSV_Plateholes.alias()
        # phm = SDSSV_Plateholes_Meta.alias()

        xx = EROSITASupersetClusters.alias()
        x = (
            xx
            .select(
                fn.rank().over(partition_by=[xx.ero_detuid],
                               order_by=[xx.xmatch_metric.desc()]).alias('x_rank'),
                xx.ero_detuid.alias('ero_detuid'),
                xx.ls_id.alias('ls_id'),
                xx.target_has_spec.alias('target_has_spec'),
            )
            .where(
                (xx.ero_version == self.parameters['ero_version']),
                (xx.xmatch_method == self.parameters['xmatch_method']),
                (xx.xmatch_version == self.parameters['xmatch_version']),
                (xx.opt_cat == self.parameters['opt_cat']),
                (xx.xmatch_metric > self.parameters['xmatch_metric_min']),
                (xx.ero_det_like > self.parameters['det_like_min']),
            )
            .alias('x')
        )

        instrument = peewee.Value(self.instrument)
        inertial = peewee.Value(self.inertial).cast('bool')

        fibertotflux_r_max = AB2nMgy(self.parameters['fibertotmag_r_min'])
        fibertotflux_r_min = AB2nMgy(self.parameters['fibertotmag_r_max'])
        fibertotflux_z_max = AB2nMgy(self.parameters['fibertotmag_z_min'])
        fibertotflux_z_min = AB2nMgy(self.parameters['fibertotmag_z_max'])

        fibertotflux_r_min_for_cadence1 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence1'])
        fibertotflux_z_min_for_cadence1 = AB2nMgy(self.parameters['fibertotmag_z_for_cadence1'])
        fibertotflux_r_min_for_cadence2 = AB2nMgy(self.parameters['fibertotmag_r_for_cadence2'])
        gaia_g_max_for_cadence1 = self.parameters['gaia_g_max_for_cadence1']
        gaia_rp_max_for_cadence1 = self.parameters['gaia_rp_max_for_cadence1']

        flux30 = AB2nMgy(30.00)

        # #########################################################################
        # prepare the spectroscopy catalogues

        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        # SDSS DR16
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        ss16 = SDSS_DR16_SpecObj.alias()
        s16 = (
            ss16.select(
                ss16.specobjid.alias('specobjid'),
            )
            .where(
                ss16.snmedian >= spec_sn_thresh,
                ss16.zwarning == 0,
                ss16.zerr <= spec_z_err_thresh,
                ss16.zerr > 0.0,
                ss16.scienceprimary > 0,
            )
            .alias('s16')
        )

        # SDSS-IV/eFEDS March2020
        c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
        ss2020 = BHM_eFEDS_Veto.alias()
        s2020 = (
            ss2020.select(
                ss2020.pk.alias('pk'),
            )
            .where(
                ss2020.sn_median_all >= spec_sn_thresh,
                ss2020.zwarning == 0,
                ss2020.z_err <= spec_z_err_thresh,
                ss2020.z_err > 0.0,
            )
            .alias('s2020')
        )

        # SDSS-V spAll
        ssV = SDSSV_BOSS_SPALL.alias()
        sV = (
            ssV.select(
                ssV.specobjid.alias('specobjid'),
                ssV.plug_ra.alias('plug_ra'),
                ssV.plug_dec.alias('plug_dec'),
            )
            .where(
                ssV.sn_median_all >= spec_sn_thresh,
                ssV.zwarning == 0,
                ssV.z_err <= spec_z_err_thresh,
                ssV.z_err > 0.0,
                ssV.specprimary > 0,
            )
            .alias('sV')
        )

        # All eFEDS plates have been observed so ignore plateholes now
        # # SDSS-V plateholes - only consider plateholes that
        # # were drilled+shipped but that were not yet observed
        # ssph = SDSSV_Plateholes.alias()
        # ssphm = SDSSV_Plateholes_Meta.alias()
        # ssconf = SDSSV_BOSS_Conflist.alias()
        # sph = (
        #     ssph.select(
        #         ssph.pkey.alias('pkey'),
        #         ssph.target_ra.alias('target_ra'),
        #         ssph.target_dec.alias('target_dec'),
        #     )
        #     .join(
        #         ssphm,
        #         on=(ssph.yanny_uid == ssphm.yanny_uid)
        #     )
        #     .join(
        #         ssconf, JOIN.LEFT_OUTER,
        #         on=(ssphm.plateid == ssconf.plate)
        #     )
        #     .where(
        #         (ssph.holetype == 'BOSS_SHARED'),
        #         (ssph.sourcetype == 'SCI') | (ssph.sourcetype == 'STA'),
        #         ssphm.isvalid > 0,
        #         ssconf.plate.is_null(),
        #     )
        #     .alias('sph')
        # )

        # priority is determined by target rank within cluster
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:

        priority = peewee.Case(
            None,
            (
                (
                    x.c.x_rank == 1,
                    self.parameters['priority_floor_bcg']
                ),
                (
                    x.c.x_rank > 1,
                    self.parameters['priority_floor_member'] +
                    fn.least(self.parameters['priority_levels'] - 2,
                             x.c.x_rank - 2)
                ),
            ),
            None)

        value = peewee.Case(
            None,
            (
                (x.c.x_rank == 1, self.parameters['value_bcg']),
                (x.c.x_rank > 1, self.parameters['value_member']),
            ),
            None)

        # choose cadence based on fiber magnitude in r-band
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence4 = 'unknown_cadence'  # catch failures
        cadence = peewee.Case(
            None,
            (
                (
                    ((ls.fibertotflux_r > fibertotflux_r_min_for_cadence1) |
                     (ls.fibertotflux_z > fibertotflux_z_min_for_cadence1) |
                     (ls.gaia_phot_g_mean_mag.between(0.1, gaia_g_max_for_cadence1)) |
                     (ls.gaia_phot_rp_mean_mag.between(0.1, gaia_rp_max_for_cadence1))),
                    cadence1),
                (ls.fibertotflux_r > fibertotflux_r_min_for_cadence2, cadence2),
                (ls.fibertotflux_r <= fibertotflux_r_min_for_cadence2, cadence3),
            ),
            cadence4)

        # We want to switch between psfmags and fibertotmags depending on
        # ls.type parameter (PSF or extended)
        # For 'PSF' targets, we use psfmags, but for extended sources use fiber2mags
        opt_prov = peewee.Case(
            ls.type,
            (('PSF', 'ls_psfmag'),),
            'ls_fibertotmag')

        magnitude_g = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_g))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_g))).cast('float'))

        magnitude_r = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_r))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_r))).cast('float'))

        magnitude_z = peewee.Case(
            ls.type,
            (('PSF', (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.flux_z))).cast('float')),),
            (22.5 - 2.5 * fn.log10(fn.greatest(flux30, ls.fibertotflux_z))).cast('float'))

        magnitude_i = peewee.Case(
            ls.type,
            (('PSF',
              (22.5 - 2.5 * fn.log10(
                  fn.greatest(flux30, 0.5 * (ls.flux_r + ls.flux_z)))).cast('float')),),
            (22.5 - 2.5 * fn.log10(
                fn.greatest(flux30, 0.5 * (ls.fibertotflux_r + ls.fibertotflux_z)))).cast('float'))

        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        query = (
            c.select(
                c.catalogid.alias('catalogid'),
                ls.ls_id.alias('ls_id'),
                x.c.ero_detuid.cast('text').alias('ero_detuid'),
                c.ra.alias('ra'),
                c.dec.alias('dec'),
                priority.alias('priority'),
                value.cast('float').alias('value'),
                cadence.alias('cadence'),
                instrument.alias('instrument'),
                magnitude_g.alias('g'),
                magnitude_r.alias('r'),
                magnitude_i.alias('i'),
                magnitude_z.alias('z'),
                opt_prov.alias('optical_prov'),
                inertial.alias('inertial'),
            )
            .join(c2ls)
            .join(ls)
            .join(x, on=(ls.ls_id == x.c.ls_id))
            # start joining the spectroscopy
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid) &
                    (c2s16.version_id == version_id)
                )
            )
            .switch(c)
            .join(c2s2020, JOIN.LEFT_OUTER)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    (c2s2020.target_id == s2020.c.pk) &
                    (c2s2020.version_id == version_id)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            # .join(
            #     sph, JOIN.LEFT_OUTER,
            #     on=(
            #         fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
            #                     c.ra, c.dec,
            #                     match_radius_spectro)
            #     )
            # )
            # finished joining the spectroscopy
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                c2ls.best >> True
            )
            .where(
                s16.c.specobjid.is_null(True),  # all of these must be satisfied
                s2020.c.pk.is_null(True),
                sV.c.specobjid.is_null(True),
                # sph.c.pkey.is_null(True),
            )
            .where(
                (
                    (ls.fibertotflux_r.between(fibertotflux_r_min, fibertotflux_r_max)) |
                    (ls.fibertotflux_z.between(fibertotflux_z_min, fibertotflux_z_max))
                ),
                (x.c.target_has_spec == 0),
                # gaia safety checks to avoid bad ls photometry
                ~(ls.gaia_phot_g_mean_mag.between(0.1, self.parameters['gaia_g_mag_limit'])),
                ~(ls.gaia_phot_rp_mean_mag.between(0.1, self.parameters['gaia_rp_mag_limit'])),
            )
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query

#
# END BhmSpidersClustersEfedsStragglersCarton
# ##################################################################################


class BhmSpidersClustersPs1dr2Carton(BaseCarton):

    name = 'bhm_spiders_clusters_ps1dr2'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    instrument = 'BOSS'

    def build_query(self, version_id, query_region=None):

        c = Catalog.alias()
        ps = Panstarrs1.alias()
        c2ps = CatalogToPanstarrs1.alias()   # only exists after v0.5 cross-match
        # s2020 = BHM_eFEDS_Veto.alias()
        # sV = SDSSV_BOSS_SPALL.alias()

        xx = EROSITASupersetClusters.alias()
        x = (
            xx
            .select(
                fn.rank().over(partition_by=[xx.ero_detuid],
                               order_by=[xx.xmatch_metric.desc()]).alias('x_rank'),
                xx.ero_detuid.alias('ero_detuid'),
                xx.ps1_dr2_id.alias('ps1_dr2_id'),
                xx.target_has_spec.alias('target_has_spec'),
            )
            .where(
                (xx.ero_version == self.parameters['ero_version']),
                (xx.xmatch_method == self.parameters['xmatch_method']),
                (xx.xmatch_version == self.parameters['xmatch_version']),
                (xx.opt_cat == self.parameters['opt_cat']),
                (xx.xmatch_metric > self.parameters['xmatch_metric_min']),
                (xx.ero_det_like > self.parameters['det_like_min']),
            )
            .alias('x')
        )

        instrument = peewee.Value(self.instrument)
        inertial = peewee.Value(self.inertial).cast('bool')

        r_psf_flux_max = AB2Jy(self.parameters['r_psf_mag_min'])
        i_psf_flux_max = AB2Jy(self.parameters['i_psf_mag_min'])
        z_psf_flux_max = AB2Jy(self.parameters['z_psf_mag_min'])
        r_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['r_psf_mag_max_for_cadence1'])
        i_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['i_psf_mag_max_for_cadence1'])
        z_psf_flux_min_for_cadence1 = AB2Jy(self.parameters['z_psf_mag_max_for_cadence1'])
        r_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['r_psf_mag_max_for_cadence2'])
        i_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['i_psf_mag_max_for_cadence2'])
        z_psf_flux_min_for_cadence2 = AB2Jy(self.parameters['z_psf_mag_max_for_cadence2'])

        # match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0

        # #########################################################################
        # prepare the spectroscopy catalogues

        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        # SDSS DR16
        c2s16 = CatalogToSDSS_DR16_SpecObj.alias()
        ss16 = SDSS_DR16_SpecObj.alias()
        s16 = (
            ss16.select(
                ss16.specobjid.alias('specobjid'),
            )
            .where(
                ss16.snmedian >= spec_sn_thresh,
                ss16.zwarning == 0,
                ss16.zerr <= spec_z_err_thresh,
                ss16.zerr > 0.0,
                ss16.scienceprimary > 0,
            )
            .alias('s16')
        )

        # SDSS-IV/eFEDS March2020
        c2s2020 = CatalogToBHM_eFEDS_Veto.alias()
        ss2020 = BHM_eFEDS_Veto.alias()
        s2020 = (
            ss2020.select(
                ss2020.pk.alias('pk'),
            )
            .where(
                ss2020.sn_median_all >= spec_sn_thresh,
                ss2020.zwarning == 0,
                ss2020.z_err <= spec_z_err_thresh,
                ss2020.z_err > 0.0,
            )
            .alias('s2020')
        )

        # SDSS-V spAll
        ssV = SDSSV_BOSS_SPALL.alias()
        sV = (
            ssV.select(
                ssV.specobjid.alias('specobjid'),
                ssV.plug_ra.alias('plug_ra'),
                ssV.plug_dec.alias('plug_dec'),
            )
            .where(
                ssV.sn_median_all >= spec_sn_thresh,
                ssV.zwarning == 0,
                ssV.z_err <= spec_z_err_thresh,
                ssV.z_err > 0.0,
                ssV.specprimary > 0,
            )
            .alias('sV')
        )

        # SDSS-V plateholes - only consider plateholes that
        # were drilled+shipped but that were not yet observed
        ssph = SDSSV_Plateholes.alias()
        ssphm = SDSSV_Plateholes_Meta.alias()
        ssconf = SDSSV_BOSS_Conflist.alias()
        sph = (
            ssph.select(
                ssph.pkey.alias('pkey'),
                ssph.target_ra.alias('target_ra'),
                ssph.target_dec.alias('target_dec'),
            )
            .join(
                ssphm,
                on=(ssph.yanny_uid == ssphm.yanny_uid)
            )
            .join(
                ssconf, JOIN.LEFT_OUTER,
                on=(ssphm.plateid == ssconf.plate)
            )
            .where(
                (ssph.holetype == 'BOSS_SHARED'),
                (ssph.sourcetype == 'SCI') | (ssph.sourcetype == 'STA'),
                ssphm.isvalid > 0,
                ssconf.plate.is_null(),
            )
            .alias('sph')
        )

        # priority is determined by target rank within cluster
        # start with a priority floor value (per carton)
        # then increment if any conditions are met:

        priority = peewee.Case(
            None,
            (
                (
                    x.c.x_rank == 1,
                    self.parameters['priority_floor_bcg']
                ),
                (
                    x.c.x_rank > 1,
                    self.parameters['priority_floor_member'] +
                    fn.least(self.parameters['priority_levels'] - 2,
                             x.c.x_rank - 2)
                ),
            ),
            None)

        value = peewee.Case(
            None,
            (
                (x.c.x_rank == 1, self.parameters['value_bcg']),
                (x.c.x_rank > 1, self.parameters['value_member']),
            ),
            None)

        # choose cadence based on psf_flux magnitude in panstarrs1 g,r,i-bands
        cadence1 = self.parameters['cadence1']
        cadence2 = self.parameters['cadence2']
        cadence3 = self.parameters['cadence3']
        cadence4 = 'unknown_cadence'
        cadence = peewee.Case(
            None,
            (
                ((ps.r_stk_psf_flux > r_psf_flux_min_for_cadence1) |
                 (ps.i_stk_psf_flux > i_psf_flux_min_for_cadence1) |
                 (ps.z_stk_psf_flux > z_psf_flux_min_for_cadence1), cadence1),
                ((ps.r_stk_psf_flux > r_psf_flux_min_for_cadence2) |
                 (ps.i_stk_psf_flux > i_psf_flux_min_for_cadence2) |
                 (ps.z_stk_psf_flux > z_psf_flux_min_for_cadence2), cadence2),
                ((ps.r_stk_psf_flux <= r_psf_flux_min_for_cadence2) &
                 (ps.i_stk_psf_flux <= i_psf_flux_min_for_cadence2) &
                 (ps.z_stk_psf_flux <= z_psf_flux_min_for_cadence2), cadence3),
            ),
            cadence4)
        # We want to switch between psfmags and fibertotmags depending on
        # ps.flags EXT+EXT_ALT (i.e. extended sources)
        # For non-extended targets, we use psfmags, but for extended sources use apermag
        flux30 = AB2Jy(30.00)
        ps1_ext_flags = 8388608 + 16777216
        ps1_good_stack_flag = 134217728
        opt_prov = peewee.Case(
            ps.flags.bin_and(ps1_ext_flags),
            ((0, 'ps_psfmag'),),
            'ps_apermag')

        magnitude_g = peewee.Case(
            ps.flags.bin_and(ps1_ext_flags),
            ((0, (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.g_stk_psf_flux))).cast('float')),),
            (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.g_stk_aper_flux))).cast('float'))

        magnitude_r = peewee.Case(
            ps.flags.bin_and(ps1_ext_flags),
            ((0, (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.r_stk_psf_flux))).cast('float')),),
            (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.r_stk_aper_flux))).cast('float'))

        magnitude_i = peewee.Case(
            ps.flags.bin_and(ps1_ext_flags),
            ((0, (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.i_stk_psf_flux))).cast('float')),),
            (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.i_stk_aper_flux))).cast('float'))

        magnitude_z = peewee.Case(
            ps.flags.bin_and(ps1_ext_flags),
            ((0, (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.z_stk_psf_flux))).cast('float')),),
            (8.9 - 2.5 * fn.log10(fn.greatest(flux30, ps.z_stk_aper_flux))).cast('float'))

        # these control matching to spectroscopy
        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0
        spec_sn_thresh = self.parameters['spec_sn_thresh']
        spec_z_err_thresh = self.parameters['spec_z_err_thresh']

        query = (
            c.select(
                c.catalogid.alias('catalogid'),
                ps.catid_objid.alias('ps1_catid_objid'),
                x.c.ero_detuid.cast('text').alias('ero_detuid'),
                c.ra.alias('ra'),
                c.dec.alias('dec'),
                priority.alias('priority'),
                value.cast('float').alias('value'),
                cadence.alias('cadence'),
                instrument.alias('instrument'),
                magnitude_g.alias('g'),
                magnitude_r.alias('r'),
                magnitude_i.alias('i'),
                magnitude_z.alias('z'),
                opt_prov.alias('optical_prov'),
                (ps.flags.bin_and(ps1_good_stack_flag) > 0)
                .cast('bool').alias('ps1_good_stack_flag'),
                inertial.alias('inertial'),
            )
            .join(c2ps)
            .join(ps)
            .join(x, on=(ps.catid_objid == x.c.ps1_dr2_id))
            # start joining the spectroscopy
            .switch(c)
            .join(c2s16, JOIN.LEFT_OUTER)
            .join(
                s16, JOIN.LEFT_OUTER,
                on=(
                    (c2s16.target_id == s16.c.specobjid) &
                    (c2s16.version_id == version_id)
                )
            )
            .switch(c)
            .join(c2s2020, JOIN.LEFT_OUTER)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    (c2s2020.target_id == s2020.c.pk) &
                    (c2s2020.version_id == version_id)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.c.plug_ra, sV.c.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            .join(
                sph, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sph.c.target_ra, sph.c.target_dec,
                                c.ra, c.dec,
                                match_radius_spectro)
                )
            )
            # finished joining the spectroscopy
            .where(
                c.version_id == version_id,
                c2ps.version_id == version_id,
                c2ps.best >> True
            )
            .where(
                s16.c.specobjid.is_null(True),  # all of these must be satisfied
                s2020.c.pk.is_null(True),
                sV.c.specobjid.is_null(True),
                sph.c.pkey.is_null(True),
            )
            .where(
                (x.c.target_has_spec == 0),
                (ps.r_stk_psf_flux < r_psf_flux_max),
                (ps.i_stk_psf_flux < i_psf_flux_max),
                (ps.z_stk_psf_flux < z_psf_flux_max),
                (ps.r_stk_psf_flux != 'NaN'),   # TODO check this is correct test via peewee
                (ps.i_stk_psf_flux != 'NaN'),
                (ps.z_stk_psf_flux != 'NaN'),
                # TODO - check panstarrs photometry quality ??
                # (ps.flags.bin_and(ps1_good_stack_flag) > 0),
                # TODO gaia safety checks to avoid bad ls photometry???
            )
            .order_by(x.c.ps1_dr2_id, x.c.x_rank.asc())
            .distinct([x.c.ps1_dr2_id, ])   # avoid duplicate entries
        )

        if query_region:
            query = query.where(peewee.fn.q3c_radial_query(c.ra, c.dec,
                                                           query_region[0],
                                                           query_region[1],
                                                           query_region[2]))

        return query

#
# END BhmSpidersClustersPs1dr2Carton
# ##################################################################################
