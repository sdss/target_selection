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
from target_selection.mag_flux import AB2nMgy  # , AB2Jy

# general catalogdb imports
from sdssdb.peewee.sdss5db.catalogdb import (
    Catalog,
    EROSITASupersetClusters,
)

# imports of existing spectro catalogues
from sdssdb.peewee.sdss5db.catalogdb import (
    # CatalogToSDSS_DR16_SpecObj,
    # SDSS_DR16_SpecObj,
    BHM_eFEDS_Veto,
    SDSSV_BOSS_SPALL,
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
    # CatalogToPanstarrs1,    # only exists after v0.5 cross-match
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

    '''

    DROP TABLE sandbox.temp_td_bhm_spiders_clusters_lsdr8;

    SELECT
         c.catalogid,
         ls.ls_id,
         x.ero_detuid,
         c.ra,
         c.dec,
         (22.5-2.5*log10(greatest(1e-3,ls.fiberflux_g))) AS fibermag_g,
         (22.5-2.5*log10(greatest(1e-3,ls.fiberflux_r))) AS fibermag_r,
         (22.5-2.5*log10(greatest(1e-3,ls.fiberflux_z))) AS fibermag_z,
         (22.5-2.5*log10(greatest(1e-3,ls.fibertotflux_g))) AS fibertotmag_g,
         (22.5-2.5*log10(greatest(1e-3,ls.fibertotflux_r))) AS fibertotmag_r,
         (22.5-2.5*log10(greatest(1e-3,ls.fibertotflux_z))) AS fibertotmag_z,
         (22.5-2.5*log10(greatest(1e-3,ls.flux_g))) AS mag_g,
         (22.5-2.5*log10(greatest(1e-3,ls.flux_r))) AS mag_r,
         (22.5-2.5*log10(greatest(1e-3,ls.flux_z))) AS mag_z,
         x_rank,
         (CASE WHEN x_rank = 1 THEN 1500
               ELSE 1550+LEAST(29,x_rank-2) END) AS priority,
         (CASE WHEN x_rank = 1 THEN 5.0 ELSE 1.0 END) AS value,
         x.target_priority AS x_priority,
         (CASE WHEN (     ls.fibertotflux_r > 251.189
                       OR ls.fibertotflux_z > 251.189
                       OR ls.gaia_phot_g_mean_mag BETWEEN 0.1 AND 16.0
                       OR ls.gaia_phot_rp_mean_mag BETWEEN 0.1 AND 16.0
                     ) THEN 'bright_2x1'
               WHEN ls.fibertotflux_r > 39.8107 THEN 'dark_1x2'
               WHEN ls.fibertotflux_r <= 39.8107 THEN 'dark_1x4'
               ELSE 'unknown_cadence' END) AS cadence,
         ls.gaia_phot_g_mean_mag AS gaia_g,
         ls.gaia_phot_rp_mean_mag AS gaia_rp,
         (CASE WHEN ls.type = 'PSF' THEN 'ls_psfmag' ELSE 'ls_fibertotmag' END) AS opt_prov,
         x.xmatch_metric,
         x.xmatch_flags,
         x.ero_det_like
    INTO sandbox.temp_td_bhm_spiders_clusters_lsdr8
    FROM catalogdb.catalog AS c
    JOIN catalog_to_legacy_survey_dr8 AS c2ls
         ON c.catalogid = c2ls.catalogid
    JOIN legacy_survey_dr8 AS ls
         ON c2ls.target_id = ls.ls_id
    JOIN ( SELECT
              xx.*, (RANK() OVER (partition BY xx.ero_detuid ORDER BY xx.xmatch_metric DESC)) as x_rank
            FROM erosita_superset_clusters as xx
            WHERE xx.ero_version = 'em01_c946_201008_poscorr'
              AND xx.xmatch_method = 'EROMAPPER_LS_DR8'
              AND xx.xmatch_version = 'grzw1_v0.3_2020-12-04'
              AND xx.opt_cat = 'lsdr8'
              AND xx.ero_det_like > 8.0
         ) AS x
         ON x.ls_id = ls.ls_id
    LEFT OUTER JOIN bhm_efeds_veto AS s2020
          ON ( q3c_join(s2020.plug_ra,s2020.plug_dec,c.ra,c.dec,1.0/3600.)
               AND s2020.zwarning = 0
               AND s2020.sn_median_all > 2.0
               AND s2020.z_err < 0.01
               AND s2020.z_err > 0.0  )
    LEFT OUTER JOIN sdssv_boss_spall AS sV
          ON ( q3c_join(sV.plug_ra,sV.plug_dec,c.ra,c.dec,1.0/3600.)
               AND sV.zwarning = 0
               AND sV.sn_median_all > 2.0
               AND sV.z_err < 0.01
               AND sV.z > 0.0  )
    WHERE
            x.target_has_spec = 0
        AND ((ls.fibertotflux_r BETWEEN 3.98107 AND 3981.07 ) OR
             (ls.fibertotflux_z BETWEEN 10.0 AND 3981.07))
        AND (ls.gaia_phot_g_mean_mag NOT BETWEEN 0.1 AND 13.5 )
        AND (ls.gaia_phot_rp_mean_mag NOT BETWEEN 0.1 AND 13.5 )
        AND c.version_id = 21
        AND c2ls.version_id = 21
        AND c2ls.best IS TRUE
        AND s2020.pk  IS NULL
        AND sV.specobjid IS NULL
    ;



    select cadence,count(*) from sandbox.temp_td_bhm_spiders_clusters_lsdr8 group by cadence;
    select count(*) from sandbox.temp_td_bhm_spiders_clusters_lsdr8 where fibertotmag_r < 17;
    select priority,count(*)
      from sandbox.temp_td_bhm_spiders_clusters_lsdr8
      group by priority order by priority;

    select cadence,count(*) from sandbox.temp_bhm_spiders_clusters_lsdr8 group by cadence;
    select priority,count(*)
      from sandbox.temp_bhm_spiders_clusters_lsdr8
      group by priority order by priority;

    \copy (SELECT * FROM sandbox.temp_td_bhm_spiders_clusters_lsdr8) TO '/home/dwelly/scratch/bhm_spiders_clusters_lsdr8.csv' with csv header

    stilts tpipe \
    in=/home/dwelly/scratch/bhm_spiders_clusters_lsdr8.csv \
    out=/home/dwelly/scratch/bhm_spiders_clusters_lsdr8.fits ofmt=fits-basic

    '''

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
                (xx.xmatch_metric >= self.parameters['xmatch_metric_min']),
                (xx.ero_det_like > self.parameters['det_like_min']),
            )
            .alias('x')
        )

        instrument = peewee.Value(self.instrument)

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
        match_radius_spectro = self.parameters['spec_join_radius'] / 3600.0

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
                    self.parameters['priority_floor_member'] + fn.least(29, x.c.x_rank - 2)
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
                x.c.ero_detuid.alias('ero_detuid'),
                c.ra.alias('ra'),
                c.dec.alias('dec'),
                priority.alias('priority'),
                value.alias('value'),
                cadence.alias('cadence'),
                instrument.alias('instrument'),
                magnitude_g.alias('g'),
                magnitude_r.alias('r'),
                magnitude_i.alias('i'),
                magnitude_z.alias('z'),
                opt_prov.alias('opt_prov'),
            )
            .join(c2ls)
            .join(ls)
            .join(x, on=(ls.ls_id == x.c.ls_id))
            .switch(c)
            .join(
                s2020, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(s2020.plug_ra, s2020.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (s2020.sn_median_all >= spec_sn_thresh) &
                    (s2020.zwarning == 0) &
                    (s2020.z_err <= spec_z_err_thresh) &
                    (s2020.z_err > 0.0)
                )
            )
            .join(
                sV, JOIN.LEFT_OUTER,
                on=(
                    fn.q3c_join(sV.plug_ra, sV.plug_dec,
                                c.ra, c.dec,
                                match_radius_spectro) &
                    (sV.sn_median_all >= spec_sn_thresh) &
                    (sV.zwarning == 0) &
                    (sV.z_err <= spec_z_err_thresh) &
                    (sV.z_err > 0.0)
                )
            )
            .where(
                c.version_id == version_id,
                c2ls.version_id == version_id,
                c2ls.best >> True
            )
            .where(
                (
                    (ls.fibertotflux_r.between(fibertotflux_r_min, fibertotflux_r_max)) |
                    (ls.fibertotflux_z.between(fibertotflux_z_min, fibertotflux_z_max))
                ),
                (x.c.target_has_spec == 0),
                (s2020.pk.is_null(True)),
                (sV.specobjid.is_null(True)),
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



#
#
#
# '''
# Pseudo SQL (optional):
#
#     all Cartons selections use the following pseudo-SQL:
#         SELECT * from bhm_spiders_clusters_superset AS x
#         INNER JOIN legacy_survey_dr8 AS ls ON  x.ls_id = ls.ls_id
#         LEFT JOIN bhm_efeds_veto AS v ON q3c_join(x.opt_ra,x.opt_dec,v.plug_ra,v.plug_dec,1.0)
#         LEFT JOIN sdss_specobj_dr16 AS sp ON q3c_join(x.opt_ra,x.opt_dec,so.plug_ra,so.plug_dec,1.0)
#         WHERE x.ero_version = "efeds_c940_clus"
#         AND WHERE x.ero_det_like > X.X
#         AND WHERE (ls.fiberflux_r > D.D OR ls.fiberflux_z > E.E)    # faint limits
#         AND WHERE ls.fiberflux_r < CCC.C                                               # bright limit
#         AND WHERE (v.plug_ra = NULL  OR v.sn_median_all < 1.x  OR v.zwarning > 0  OR v.z_err > 2e-3  OR v.z_err <= 0.0)
#         AND WHERE (so.plug_ra = NULL OR so.sn_median_all < 1.x OR so.zwarning > 0 OR so.z_err > 2e-3 OR so.z_err <= 0.0 OR so.z_err >so.z )
#     Carton-specific pseudo-SQL based on targetting bits in superset catalogue:
#         bhm_spiders_clusters_efeds_sdss_redmapper
#             AND WHERE x.xmatch_flags & (2^2 + 2^6 + 2^10 + 2^12 + 2^17 + 2^19) !=0
#         bhm_spiders_clusters_efeds_hsc_redmapper
#             AND WHERE x.xmatch_flags & (2^5 + 2^11) !=0
#         bhm_spiders_clusters_efeds_ls_redmapper
#             AND WHERE x.xmatch_flags & (2^1 + 2^3 + 2^8 + 2^15) !=0
#         bhm_spiders_clusters_efeds_erosita
#             AND WHERE x.xmatch_flags & (2^2 + 2^6 + 2^10 + 2^12 + 2^17    +    2^19 + 2^5 + 2^11  +    2^1 + 2^3 + 2^8 + 2^15)  = 0
# '''
#
#
# class BhmSpidersClusEfedsCarton(BaseCarton):
#     ''' Parent class that provides the underlying selections for all SPIDERS Clusters eFEDS cartons'''
#
#     name = 'bhm_spiders_clusters-efeds-base'
#     base_name = 'bhm_spiders_clusters-efeds-base'
#     category = 'science'
#     mapper = 'BHM'
#     program = 'bhm_spiders'
#     tile = False
#     priority = None
#     cadence = 'bhm_spiders_1x8'
#     bitmask = None
#     mask_sense = True
#
#     def build_query(self, version_id, query_region=None):
#         c = Catalog.alias()
#         x = BHM_Spiders_Clusters_Superset.alias()
#         ls = Legacy_Survey_DR8.alias()
#         c2ls = CatalogToLegacy_Survey_DR8.alias()
#         v = BHM_eFEDS_Veto.alias()
#         s = SDSS_DR16_SpecObj.alias()
#
#         base_parameters = self.config['parameters'].get(self.base_name, None)
#         assert base_parameters is not None
#
#
#         fiberflux_r_max = AB2nMgy(base_parameters['fibermag_r_min'])
#         fiberflux_r_min = AB2nMgy(base_parameters['fibermag_r_max'])
#         fiberflux_z_min = AB2nMgy(base_parameters['fibermag_z_max'])
#
#         flux30 = AB2nMgy(30.0)
#
#         value = peewee.Value(base_parameters.get('value', 1.0)).cast('float').alias('value')
#         match_radius_spectro = base_parameters['spec_join_radius']/3600.0
#
#         p_f = base_parameters['priority_floor']
#         priority = peewee.Case(None,
#                                (
#                                    ((x.target_priority == 0), p_f+0), # BCGs
#                                ),
#                                p_f + x.target_priority + 10)     # Member galaxies
#
#         # legacysurvey mags - derived from fiberfluxes - with limits to avoid divide by zero errors
#         # notes on convertion from ls_fibermag to sdss_fiber2mag:
#         # https://wiki.mpe.mpg.de/eRosita/EroAGN_eFEDS/SDSSIVSpecialPlates#Estimating_SDSS_fiber2mag_.2A_from_legacysurvey_photometry
#         # for mostly point-like sources, the average offset is sdss_fiber2mag_[griz] = ls_fibermag_[griz] + 0.55 mag
#         # for mostly galaxies, there is a magnitude and band-dependent shift
#         # sdss_fiber2mag_g = ls_fibermag_g + 0.46 mag   flux_ratio = ~0.65
#         # sdss_fiber2mag_r = ls_fibermag_r + 0.55 mag   flux_ratio = ~0.60
#         # sdss_fiber2mag_i = ls_fibermag_i + 0.44 mag   flux_ratio = ~0.67
#         # sdss_fiber2mag_z = ls_fibermag_z + 0.39 mag   flux_ratio = ~0.70
#         flux_ratio = {'g' : 0.65, 'r' : 0.60, 'i' : 0.67, 'z' : 0.70 }
#         # Then add the correction from sdss_fiber2mag to sdss_psfmag
#
#         # Notes on converting from sdss_fiber2mag to sdss_psfmag
#         # https://wiki.sdss.org/display/OPS/Contents+of+targetdb.magnitude#Contentsoftargetdb.magnitude-WhatmagnitudestoputintotheplPlugMapfilesforBOSSplatetargets?
#
#         magnitude_g = (psfmag_minus_fiber2mag('g') +
#                        22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_g*flux_ratio['g']))).cast('float')
#         magnitude_r = (psfmag_minus_fiber2mag('r') +
#                        22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_r*flux_ratio['r']))).cast('float')
#         magnitude_z = (psfmag_minus_fiber2mag('z') +
#                        22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_z*flux_ratio['z']))).cast('float')
#         # the simplest possible interpolation - TODO do this better
#         magnitude_i = (psfmag_minus_fiber2mag('i') +
#                        22.5-2.5*fn.log10(fn.greatest(flux30,
#                                                      0.5*(ls.fiberflux_r+
#                                                           ls.fiberflux_z)*flux_ratio['i']))).cast('float')
#
#
#
#         pmra = peewee.Value(0.0).cast('float').alias('pmra')
#         pmdec = peewee.Value(0.0).cast('float').alias('pmdec')
#         parallax = peewee.Value(0.0).cast('float').alias('parallax')
#
#         query = (
#             c
#             .select(c.catalogid,
#                     priority.alias('priority'),
#                     value,
#                     pmra,
#                     pmdec,
#                     parallax,
#                     magnitude_g.alias("g"),
#                     magnitude_r.alias("r"),
#                     magnitude_i.alias("i"),
#                     magnitude_z.alias("z"),
#             )
#             .join(c2ls)
#             .join(ls)
#             .join(x)
#             .join(v, JOIN.LEFT_OUTER,
#                   on=(fn.q3c_join(c.ra,c.dec,
#                                  v.plug_ra,v.plug_dec,
#                                  match_radius_spectro) &
#                       (v.sn_median[3] >= base_parameters['spec_sn_thresh']) &
#                       (v.zwarning == 0) &
#                       (v.z_err <= base_parameters['spec_z_err_thresh']) &
#                       (v.z_err > 0.0)
#                       )
#                   )
#             .join(s, JOIN.LEFT_OUTER,
#                   on=(fn.q3c_join(c.ra,c.dec,
#                                   s.ra,s.dec,
#                                   match_radius_spectro) &
#                       (s.snmedian_i >= base_parameters['spec_sn_thresh']) &
#                       (s.zwarning == 0) &
#                       (s.zerr <= base_parameters['spec_z_err_thresh']) &
#                       (s.zerr > 0.0) &
#                       (s.scienceprimary > 0)
#                   )
#             )
#             .where(c.version_id == version_id,
#                    c2ls.version_id == version_id,
#                    c2ls.best == True)
#             .where(
#                 (x.ero_version == base_parameters['ero_version'] ),
#                 (v.pk.is_null()),
#                 (s.specobjid.is_null()),
#                 (ls.fibertotflux_r < fiberflux_r_max),
#                 (
#                     (ls.fiberflux_r >= fiberflux_r_min) |
#                     (ls.fiberflux_z >= fiberflux_z_min)
#                 ),
#                 (x.ero_det_like > base_parameters['det_like_min']),
#                 (x.target_has_spec == 0),
#             )
#             .distinct([ls.ls_id])   # avoid duplicates - trust the ls_id
#         )
#         if self.mask_sense is True:
#             query = query.where(x.xmatch_flags.bin_and(self.bitmask) != 0)
#         else:
#             query = query.where(x.xmatch_flags.bin_and(self.bitmask) == 0)
#
#         return query
#
#
# #-------BHM SPIDERS eFEDS Clusters SDSS RedMapper------ #
#
# class BhmSpidersClusEfedsSdssRedmapperCarton(BhmSpidersClusEfedsCarton):
#     '''
#     Carton-specific pseudo-SQL based on targetting bits in superset catalogue:
#             AND WHERE x.xmatch_flags & (2^2 + 2^6 + 2^10 + 2^12 + 2^17 + 2^19) !=0
#     '''
#
#     name = 'bhm_spiders_clusters-efeds-sdss-redmapper'
#     bitmask = 2**2 + 2**6 + 2**10 + 2**12 + 2**17 + 2**19
#
#
# #-------BHM SPIDERS eFEDS Clusters HSC RedMapper------ #
#
# class BhmSpidersClusEfedsHscRedmapperCarton(BhmSpidersClusEfedsCarton):
#     '''
#     Carton-specific pseudo-SQL based on targetting bits in superset catalogue:
#             AND WHERE x.xmatch_flags & (2^5 + 2^11) !=0
#     '''
#
#     name = 'bhm_spiders_clusters-efeds-hsc-redmapper'
#     bitmask = 2**5 + 2**11
#
#
# #-------BHM SPIDERS eFEDS Clusters LS RedMapper------ #
#
# class BhmSpidersClusEfedsLsRedmapperCarton(BhmSpidersClusEfedsCarton):
#     '''
#     Carton-specific pseudo-SQL based on targetting bits in superset catalogue:
#             AND WHERE x.xmatch_flags & (2^1 + 2^3 + 2^8 + 2^15) !=0
#     '''
#
#     name = 'bhm_spiders_clusters-efeds-ls-redmapper'
#     bitmask = 2**1 + 2**3 + 2**8 + 2**15
#
#
# #-------BHM SPIDERS eFEDS Clusters eROSITA------ #
#
# class BhmSpidersClusEfedsErositaCarton(BhmSpidersClusEfedsCarton):
#     '''
#     Carton-specific pseudo-SQL based on targetting bits in superset catalogue:
#     AND WHERE x.xmatch_flags & (2^2 + 2^6 + 2^10 + 2^12 + 2^17    +
#                                 2^19 + 2^5 + 2^11  +
#                                 2^1 + 2^3 + 2^8 + 2^15)
#                                 = 0
#     '''
#
#     name = 'bhm_spiders_clusters-efeds-erosita'
#     bitmask =  (2**2 + 2**6 + 2**10 + 2**12 + 2**17 + 2**19 +
#                 2**5 + 2**11 +
#                 2**1 + 2**3 + 2**8 + 2**15)
#     mask_sense = False
#
#
#
#
#
# '''
#
#  target_selection --profile tunnel_operations_sdss --verbose run --include bhm_spiders_clusters_efeds_sdss_redmapper,bhm_spiders_clusters_efeds_hsc_redmapper,bhm_spiders_clusters_efeds_ls_redmapper,bhm_spiders_clusters_efeds_erosita --keep --overwrite '0.1.0' --write-table
#
#
# Exporting from the temp table
#
# \copy (SELECT * FROM sandbox.temp_bhm_spiders_clusters_efeds_sdss_redmapper)  TO '/home/tdwelly/scratch/targetdb/bhm_spiders_clusters_efeds_sdss_redmapper.csv' with csv header
# \copy (SELECT * FROM sandbox.temp_bhm_spiders_clusters_efeds_hsc_redmapper)  TO '/home/tdwelly/scratch/targetdb/bhm_spiders_clusters_efeds_hsc_redmapper.csv' with csv header
# \copy (SELECT * FROM sandbox.temp_bhm_spiders_clusters_efeds_ls_redmapper)  TO '/home/tdwelly/scratch/targetdb/bhm_spiders_clusters_efeds_ls_redmapper.csv' with csv header
# \copy (SELECT * FROM sandbox.temp_bhm_spiders_clusters_efeds_erosita)  TO '/home/tdwelly/scratch/targetdb/bhm_spiders_clusters_efeds_erosita.csv' with csv header
#
#
# for F in bhm_spiders_clusters_*.csv; do   stilts tpipe in=${F} out="${F%.*}.fits" ifmt=csv ofmt=fits-basic; done
#
# '''
#
