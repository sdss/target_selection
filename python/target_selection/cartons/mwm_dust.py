#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-16
# @Filename: mwm_dust.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (GLIMPSE, AllWise, Catalog,
                                             CatalogToAllWise,
                                             CatalogToGLIMPSE, CatalogToTIC_v8,
                                             Gaia_DR2, TIC_v8, TwoMassPSC)

from . import BaseCarton


class MWM_Dust_Carton(BaseCarton):
    """MWM Dust Carton.

    Definition:

        A sample of bright, nearby, midplane giants that, when combined
        with the more distant and/or heavily reddened sample in other cartons,
        produces 100 stars per (100pc)^3 volume in the near disk.

        - Use Gaia parallax and (l,b) to get distance
        - Use H and 4.5 micron to get A_Ks and E_JKs, a la APOGEE
          (using RJCE dereddening prescription from Majewski et al. 2011).
        - Use J, K, and E_JKs to get (J-Ks)_0
        - Use distance, K, and A_Ks to get absolute mag M_K
        - 0 < Gaia parallax error/parallax < 0.2
        - M_K < 2.6
        - H < 11.2
        - distance < 5 kpc
        - |z| < 0.2 kpc
        - (J-Ks)_0 > 0.5
        - Can be spatially subselected to complement GG sampling in order to
          obtain 100 stars per (100pc)^3 volume.

    Non-SQL implementation:

        https://faun.rc.fas.harvard.edu/eschlafly/sdss5/dustsel.py

    """

    name = 'mwm_dust'
    mapper = 'MWM'
    category = 'science'
    program = 'Dust'

    def build_query(self, version_id, query_region=None):

        fn = peewee.fn

        gallong = Gaia_DR2.l
        gallat = Gaia_DR2.b

        ipar = 1. / Gaia_DR2.parallax / 1000  # kpc
        z = ipar * fn.sin(fn.radians(gallat))
        x = ipar * fn.cos(fn.radians(gallat)) * fn.cos(fn.radians(gallong))
        y = ipar * fn.cos(fn.radians(gallat)) * fn.sin(fn.radians(gallong))

        aks_glimpse = 0.918 * (GLIMPSE.mag_h - GLIMPSE.mag4_5 - 0.08)
        aks_allwise = 0.918 * (AllWise.h_m_2mass - AllWise.w2mpro - 0.08)
        aks = fn.coallesce(aks_glimpse, aks_allwise)

        Ej_ks = 1.5 * aks
        j_ks_glimpse = GLIMPSE.mag_j - GLIMPSE.mag_ksAllWise

        plxfracunc = Gaia_DR2. parallax_error / Gaia_DR2. parallax
        dm = 5 * fn.log(1000. / Gaia_DR2.parallax / 10)
        absmag = TwoMassPSC.k_m - aks - dm

        query = (Catalog
                 .select(Catalog.catalogid,
                         )
                 .join(CatalogToTIC_v8)
                 .join(TIC_v8)
                 .join(TwoMassPSC)
                 .join_from(TIC_v8, Gaia_DR2)
                 .join_from(Catalog, CatalogToAllWise)
                 .join(AllWise)
                 .where(CatalogToAllWise.version_id == version_id,
                        CatalogToAllWise.best >> True)
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True))

        if query_region:
            query = (query
                     .join_from(CatalogToAllWise, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query

    def post_process(self, model):

        pass
