#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-06-17
# @Filename: bhm_spiders_clusters.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

# flake8: noqa
# isort: skip_file

import peewee
from peewee import JOIN
from peewee import fn
import sdssdb

from target_selection.cartons.base import BaseCarton
from target_selection.mag_flux import psfmag_minus_fiber2mag, AB2nMgy


from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             CatalogToLegacy_Survey_DR8,
                                             Legacy_Survey_DR8,
                                             BHM_Spiders_Clusters_Superset,
                                             BHM_eFEDS_Veto,
                                             SDSS_DR16_SpecObj,
                                             )


## This provides the following BHM cartons:

# bhm_spiders_clusters-efeds-sdss-redmapper
# bhm_spiders_clusters-efeds-hsc-redmapper
# bhm_spiders_clusters-efeds-ls-redmapper
# bhm_spiders_clusters-efeds-erosita
# bhm_spiders_clusters-wide
# bhm_spiders_clusters-deep

# Details: Start here
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms


'''
Pseudo SQL (optional):

    all Cartons selections use the following pseudo-SQL:
        SELECT * from bhm_spiders_clusters_superset AS x
        INNER JOIN legacy_survey_dr8 AS ls ON  x.ls_id = ls.ls_id
        LEFT JOIN bhm_efeds_veto AS v ON q3c_join(x.opt_ra,x.opt_dec,v.plug_ra,v.plug_dec,1.0)
        LEFT JOIN sdss_specobj_dr16 AS sp ON q3c_join(x.opt_ra,x.opt_dec,so.plug_ra,so.plug_dec,1.0)
        WHERE x.ero_version = "efeds_c940_clus"
        AND WHERE x.ero_det_like > X.X
        AND WHERE (ls.fiberflux_r > D.D OR ls.fiberflux_z > E.E)    # faint limits
        AND WHERE ls.fiberflux_r < CCC.C                                               # bright limit
        AND WHERE (v.plug_ra = NULL  OR v.sn_median_all < 1.x  OR v.zwarning > 0  OR v.z_err > 2e-3  OR v.z_err <= 0.0)
        AND WHERE (so.plug_ra = NULL OR so.sn_median_all < 1.x OR so.zwarning > 0 OR so.z_err > 2e-3 OR so.z_err <= 0.0 OR so.z_err >so.z )
    Carton-specific pseudo-SQL based on targetting bits in superset catalogue:
        bhm_spiders_clusters_efeds_sdss_redmapper
            AND WHERE x.xmatch_flags & (2^2 + 2^6 + 2^10 + 2^12 + 2^17 + 2^19) !=0
        bhm_spiders_clusters_efeds_hsc_redmapper
            AND WHERE x.xmatch_flags & (2^5 + 2^11) !=0
        bhm_spiders_clusters_efeds_ls_redmapper
            AND WHERE x.xmatch_flags & (2^1 + 2^3 + 2^8 + 2^15) !=0
        bhm_spiders_clusters_efeds_erosita
            AND WHERE x.xmatch_flags & (2^2 + 2^6 + 2^10 + 2^12 + 2^17    +    2^19 + 2^5 + 2^11  +    2^1 + 2^3 + 2^8 + 2^15)  = 0
'''


class BhmSpidersClusEfedsCarton(BaseCarton):
    ''' Parent class that provides the underlying selections for all SPIDERS Clusters eFEDS cartons'''

    name = 'bhm_spiders_clusters-efeds-base'
    base_name = 'bhm_spiders_clusters-efeds-base'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_spiders'
    tile = False
    priority = None
    cadence = 'bhm_spiders_1x8'
    bitmask = None
    mask_sense = True

    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        x = BHM_Spiders_Clusters_Superset.alias()
        ls = Legacy_Survey_DR8.alias()
        c2ls = CatalogToLegacy_Survey_DR8.alias()
        v = BHM_eFEDS_Veto.alias()
        s = SDSS_DR16_SpecObj.alias()

        base_parameters = self.config['parameters'].get(self.base_name, None)
        assert base_parameters is not None


        fiberflux_r_max = AB2nMgy(base_parameters['fibermag_r_min'])
        fiberflux_r_min = AB2nMgy(base_parameters['fibermag_r_max'])
        fiberflux_z_min = AB2nMgy(base_parameters['fibermag_z_max'])

        flux30 = AB2nMgy(30.0)

        value = peewee.Value(base_parameters.get('value', 1.0)).cast('float').alias('value')
        match_radius_spectro = base_parameters['spec_join_radius']/3600.0

        p_f = base_parameters['priority_floor']
        priority = peewee.Case(None,
                               (
                                   ((x.target_priority == 0), p_f+0), # BCGs
                               ),
                               p_f + x.target_priority + 10)     # Member galaxies

        # legacysurvey mags - derived from fiberfluxes - with limits to avoid divide by zero errors
        # notes on convertion from ls_fibermag to sdss_fiber2mag:
        # https://wiki.mpe.mpg.de/eRosita/EroAGN_eFEDS/SDSSIVSpecialPlates#Estimating_SDSS_fiber2mag_.2A_from_legacysurvey_photometry
        # for mostly point-like sources, the average offset is sdss_fiber2mag_[griz] = ls_fibermag_[griz] + 0.55 mag
        # for mostly galaxies, there is a magnitude and band-dependent shift
        # sdss_fiber2mag_g = ls_fibermag_g + 0.46 mag   flux_ratio = ~0.65
        # sdss_fiber2mag_r = ls_fibermag_r + 0.55 mag   flux_ratio = ~0.60
        # sdss_fiber2mag_i = ls_fibermag_i + 0.44 mag   flux_ratio = ~0.67
        # sdss_fiber2mag_z = ls_fibermag_z + 0.39 mag   flux_ratio = ~0.70
        flux_ratio = {'g' : 0.65, 'r' : 0.60, 'i' : 0.67, 'z' : 0.70 }
        # Then add the correction from sdss_fiber2mag to sdss_psfmag

        # Notes on converting from sdss_fiber2mag to sdss_psfmag
        # https://wiki.sdss.org/display/OPS/Contents+of+targetdb.magnitude#Contentsoftargetdb.magnitude-WhatmagnitudestoputintotheplPlugMapfilesforBOSSplatetargets?

        magnitude_g = (psfmag_minus_fiber2mag('g') +
                       22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_g*flux_ratio['g']))).cast('float')
        magnitude_r = (psfmag_minus_fiber2mag('r') +
                       22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_r*flux_ratio['r']))).cast('float')
        magnitude_z = (psfmag_minus_fiber2mag('z') +
                       22.5-2.5*fn.log10(fn.greatest(flux30,ls.fiberflux_z*flux_ratio['z']))).cast('float')
        # the simplest possible interpolation - TODO do this better
        magnitude_i = (psfmag_minus_fiber2mag('i') +
                       22.5-2.5*fn.log10(fn.greatest(flux30,
                                                     0.5*(ls.fiberflux_r+
                                                          ls.fiberflux_z)*flux_ratio['i']))).cast('float')



        pmra = peewee.Value(0.0).cast('float').alias('pmra')
        pmdec = peewee.Value(0.0).cast('float').alias('pmdec')
        parallax = peewee.Value(0.0).cast('float').alias('parallax')

        query = (
            c
            .select(c.catalogid,
                    priority.alias('priority'),
                    value,
                    pmra,
                    pmdec,
                    parallax,
                    magnitude_g.alias("g"),
                    magnitude_r.alias("r"),
                    magnitude_i.alias("i"),
                    magnitude_z.alias("z"),
            )
            .join(c2ls)
            .join(ls)
            .join(x)
            .join(v, JOIN.LEFT_OUTER,
                  on=(fn.q3c_join(c.ra,c.dec,
                                 v.plug_ra,v.plug_dec,
                                 match_radius_spectro) &
                      (v.sn_median[3] >= base_parameters['spec_sn_thresh']) &
                      (v.zwarning == 0) &
                      (v.z_err <= base_parameters['spec_z_err_thresh']) &
                      (v.z_err > 0.0)
                      )
                  )
            .join(s, JOIN.LEFT_OUTER,
                  on=(fn.q3c_join(c.ra,c.dec,
                                  s.ra,s.dec,
                                  match_radius_spectro) &
                      (s.snmedian_i >= base_parameters['spec_sn_thresh']) &
                      (s.zwarning == 0) &
                      (s.zerr <= base_parameters['spec_z_err_thresh']) &
                      (s.zerr > 0.0) &
                      (s.scienceprimary > 0)
                  )
            )
            .where(c.version_id == version_id,
                   c2ls.version_id == version_id,
                   c2ls.best == True)
            .where(
                (x.ero_version == base_parameters['ero_version'] ),
                (v.pk.is_null()),
                (s.specobjid.is_null()),
                (ls.fibertotflux_r < fiberflux_r_max),
                (
                    (ls.fiberflux_r >= fiberflux_r_min) |
                    (ls.fiberflux_z >= fiberflux_z_min)
                ),
                (x.ero_det_like > base_parameters['det_like_min']),
                (x.target_has_spec == 0),
            )
            .distinct([ls.ls_id])   # avoid duplicates - trust the ls_id
        )
        if self.mask_sense is True:
            query = query.where(x.xmatch_flags.bin_and(self.bitmask) != 0)
        else:
            query = query.where(x.xmatch_flags.bin_and(self.bitmask) == 0)

        return query


#-------BHM SPIDERS eFEDS Clusters SDSS RedMapper------ #

class BhmSpidersClusEfedsSdssRedmapperCarton(BhmSpidersClusEfedsCarton):
    '''
    Carton-specific pseudo-SQL based on targetting bits in superset catalogue:
            AND WHERE x.xmatch_flags & (2^2 + 2^6 + 2^10 + 2^12 + 2^17 + 2^19) !=0
    '''

    name = 'bhm_spiders_clusters-efeds-sdss-redmapper'
    bitmask = 2**2 + 2**6 + 2**10 + 2**12 + 2**17 + 2**19


#-------BHM SPIDERS eFEDS Clusters HSC RedMapper------ #

class BhmSpidersClusEfedsHscRedmapperCarton(BhmSpidersClusEfedsCarton):
    '''
    Carton-specific pseudo-SQL based on targetting bits in superset catalogue:
            AND WHERE x.xmatch_flags & (2^5 + 2^11) !=0
    '''

    name = 'bhm_spiders_clusters-efeds-hsc-redmapper'
    bitmask = 2**5 + 2**11


#-------BHM SPIDERS eFEDS Clusters LS RedMapper------ #

class BhmSpidersClusEfedsLsRedmapperCarton(BhmSpidersClusEfedsCarton):
    '''
    Carton-specific pseudo-SQL based on targetting bits in superset catalogue:
            AND WHERE x.xmatch_flags & (2^1 + 2^3 + 2^8 + 2^15) !=0
    '''

    name = 'bhm_spiders_clusters-efeds-ls-redmapper'
    bitmask = 2**1 + 2**3 + 2**8 + 2**15


#-------BHM SPIDERS eFEDS Clusters eROSITA------ #

class BhmSpidersClusEfedsErositaCarton(BhmSpidersClusEfedsCarton):
    '''
    Carton-specific pseudo-SQL based on targetting bits in superset catalogue:
    AND WHERE x.xmatch_flags & (2^2 + 2^6 + 2^10 + 2^12 + 2^17    +
                                2^19 + 2^5 + 2^11  +
                                2^1 + 2^3 + 2^8 + 2^15)
                                = 0
    '''

    name = 'bhm_spiders_clusters-efeds-erosita'
    bitmask =  (2**2 + 2**6 + 2**10 + 2**12 + 2**17 + 2**19 +
                2**5 + 2**11 +
                2**1 + 2**3 + 2**8 + 2**15)
    mask_sense = False





'''

 target_selection --profile tunnel_operations_sdss --verbose run --include bhm_spiders_clusters_efeds_sdss_redmapper,bhm_spiders_clusters_efeds_hsc_redmapper,bhm_spiders_clusters_efeds_ls_redmapper,bhm_spiders_clusters_efeds_erosita --keep --overwrite '0.1.0' --write-table


Exporting from the temp table

\copy (SELECT * FROM sandbox.temp_bhm_spiders_clusters_efeds_sdss_redmapper)  TO '/home/tdwelly/scratch/targetdb/bhm_spiders_clusters_efeds_sdss_redmapper.csv' with csv header
\copy (SELECT * FROM sandbox.temp_bhm_spiders_clusters_efeds_hsc_redmapper)  TO '/home/tdwelly/scratch/targetdb/bhm_spiders_clusters_efeds_hsc_redmapper.csv' with csv header
\copy (SELECT * FROM sandbox.temp_bhm_spiders_clusters_efeds_ls_redmapper)  TO '/home/tdwelly/scratch/targetdb/bhm_spiders_clusters_efeds_ls_redmapper.csv' with csv header
\copy (SELECT * FROM sandbox.temp_bhm_spiders_clusters_efeds_erosita)  TO '/home/tdwelly/scratch/targetdb/bhm_spiders_clusters_efeds_erosita.csv' with csv header


for F in bhm_spiders_clusters_*.csv; do   stilts tpipe in=${F} out="${F%.*}.fits" ifmt=csv ofmt=fits-basic; done

'''
