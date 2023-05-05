#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-05
# @Filename: mwm_cb.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (HECATE_1_1, MILLIQUAS_7_7,
                                             UVOT_SSC_1, XMM_OM_SUSS_5_0,
                                             BailerJonesEDR3, Catalog,
                                             CatalogToGaia_DR3,
                                             CatalogToMILLIQUAS_7_7,
                                             CatalogToUVOT_SSC_1,
                                             CatalogToXMM_OM_SUSS_5_0,
                                             Gaia_DR3, Galex_GR7_Gaia_DR3,
                                             GUVCat)

from target_selection.cartons import BaseCarton


class MWM_CB_GALEX_Mag(BaseCarton):
    """Magnitude-limited GALEX FUV / Gaia selection.

    Simplified Description of selection criteria:
        Uses Nicola Gentile Fusillo's Gaia/GALEX cross-match (which should be already ingested).
        Select objects that lie below the main-sequence in an FUV/optical absolute magnitude-color
        diagram. Remove AGN interlopers. Remove the LMC and SMC, where a lot of hot young stars
        mimic the properties of this target selection. For more details, see here.

    Pseudo-code:
        Use the following short-hands:
            pm = sqrt(pmra*pmra+pmdec*pmdec)
            pm_error = sqrt((pmra_error*pmra/pm)*(pmra_error*pmra/pm)+
                            (pmdec_error*pmdec/pm)*(pmdec_error*pmdec/pm))
            logpmdpm = log10(pm/pm_error)
            poe = parallax / parallax_error
        Then combine the following cuts with &&
            fuv_mag + 5*log10(parallax/1000) + 5 > 1.5 + 1.28*(fuv_mag - phot_g_mean_mag)
            fuv_mag + 5*log10(parallax/1000) + 5 > 2.5*(fuv_mag - phot_g_mean_mag) - 0.5
            fuv_magerr < 0.3
            !(logpmdpm < 0.301 && poe > -1.75*logpmdpm - 2.8 && poe <  1.75*logpmdpm + 2.8 ||
              (logpmdpm-0.301)*(logpmdpm-0.301) / (0.33*0.33) + (poe * poe) / (3.2* 3.2)<=1)
            !((sqrt( pow(l-280.3,2) + 2*pow(b+33.0,2) ) < 8) ||
              sqrt(0.5*square(l-301) + 2*square(b+43.3) ) < 5)

    """

    name = 'mwm_cb_galex_mag'
    mapper = 'MWM'
    category = 'science'
    instrument = 'BOSS'
    cadence = None
    program = 'mwm_cb'
    priority = 1400
    can_offset = True

    def build_query(self, version_id, query_region=None):

        G3 = Gaia_DR3
        pm = fn.sqrt(G3.pmra * G3.pmra + G3.pmdec * G3.pmdec)
        pm_error = fn.sqrt((G3.pmra_error * G3.pmra / pm) * (G3.pmra_error * G3.pmra / pm) +
                           (G3.pmdec_error * G3.pmdec / pm) * (G3.pmdec_error * G3.pmdec / pm))
        logpmdpm = fn.log(pm / pm_error)
        poe = G3.parallax / G3.parallax_error

        fuv_mag = GUVCat.fuv_mag
        Gm = G3.phot_g_mean_mag
        parallax = G3.parallax
        logpar = fn.log(fn.abs(parallax / 1000))
        ll = G3.l
        bb = G3.b

        fuv_cut = (((fuv_mag + 5 * logpar + 5) > (1.5 + 1.28 * (fuv_mag - Gm))) &
                   ((fuv_mag + 5 * logpar + 5) > (2.5 * (fuv_mag - Gm) - 0.5)))

        astro_cut1 = ~(((logpmdpm < 0.301) &
                        (poe > (-1.75 * logpmdpm - 2.8)) &
                        (poe < (1.75 * logpmdpm + 2.8))) |
                       (((logpmdpm - 0.301) * (logpmdpm - 0.301) / (0.33 * 0.33) +
                         (poe * poe) / (3.2 * 3.2)) <= 1))

        astro_cut2 = ~((fn.sqrt(fn.pow(ll - 280.3, 2) + 2 * fn.pow(bb + 33.0, 2)) < 8) |
                       (fn.sqrt(0.5 * fn.pow(ll - 301, 2) + 2 * fn.pow(bb + 43.3, 2)) < 5))

        cadence = peewee.Case(None, ((Gm <= 16, 'bright_2x1'), (Gm > 16, 'dark_2x1')))

        query = (Gaia_DR3
                 .select(CatalogToGaia_DR3.catalogid,
                         G3.source_id,
                         ll,
                         bb,
                         parallax,
                         Gm,
                         fuv_mag,
                         cadence.alias('cadence'))
                 .join(Galex_GR7_Gaia_DR3,
                       on=(Gaia_DR3.source_id == Galex_GR7_Gaia_DR3.gaia_edr3_source_id))
                 .join(GUVCat,
                       on=(GUVCat.objid == Galex_GR7_Gaia_DR3.galex_objid))
                 .join_from(Gaia_DR3, CatalogToGaia_DR3)
                 .where(CatalogToGaia_DR3.version_id == version_id,
                        CatalogToGaia_DR3.best >> True,
                        fuv_cut,
                        fuv_mag > -999,
                        GUVCat.fuv_magerr < 0.3,
                        astro_cut1,
                        astro_cut2)
                 .order_by(CatalogToGaia_DR3.catalogid, Galex_GR7_Gaia_DR3.galex_separation)
                 .distinct(CatalogToGaia_DR3.catalogid))

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_CB_GALEX_Vol(BaseCarton):
    """Volume-limited GALEX FUV / Gaia selection.

    Simplified Description of selection criteria:
        Uses Nicola Gentile Fusillo's Gaia/GALEX cross-match (which should be already ingested).
        Select objects that lie below the main-sequence in an NUV/optical absolute magnitude-color
        diagram. Remove known AGN & galaxies from published catalogues. Remove AGN interlopers by
        parallax / proper motion. Remove the LMC and SMC, where a lot of hot young stars mimic the
        properties of this target selection. Makes use of milliquas (million quasar catalog
        https://quasars.org/milliquas.htm) and HECATE (Heraklion Extragalactic Catalog -
        https://hecate.ia.forth.gr) to remove contaminants.

    Pseudo-code:
        Select entries with maximum separation between Gaia and GALEX position of 2.5 arcsec.
        Retrieve Gaia distances from Bailer-Jones et al 2022 (referred to as either r_med_geo or
            r_geo)
        Eliminate known AGN, i.e. all milliquas entries with 100% probability being an AGN within a
            radius of 1 arcsec around Gaia position:
                java -jar stilts.jar tpipe in=milliquas.fits
                    cmd='select QPCT==100' out=milliquas_good1.fits
                java -jar -Xmx16G stilts.jar tmatch2 in1=INPUT in2=milliquas_good1.fits
                    matcher=sky params=1.0 values1='ra dec' values2='RA DEC'
                    join=1not2 find=all out=xxx)
        Eliminate nearby galaxies  (all entries in HECATE_v1.1.fits within 60*R1, R1 is
            radius of galaxies given in arcmin)
        Select entries with valid NUV magnitude only:
            select nuv_mag>-999 && nuv_magerr<0.2 && nuv_magerr>0
        Select objects with an UV excess via two cuts, one in the NUV-HRD, the other
            in the optical-UV color-color diagram:
                - (nuv_mag-5*log10(r_med_geo)+5>2*(nuv_mag-phot_g_mean_mag)+0.7) ||
                  (nuv_mag-5*log10(1000/parallax)+5>2*(nuv_mag-phot_g_mean_mag)+0.7)
                - nuv_mag - phot_g_mean_mag < 1.6 * (bp-rp) + 1.2
        Astrometric cut (same short-hands as in mwm_cb_galex_mag):
            To remove AGN contaminants:
                !(logpmdpm < 0.301 && poe > -1.75*logpmdpm - 2.8 && poe <  1.75*logpmdpm + 2.8 ||
                (logpmdpm-0.301)*(logpmdpm-0.301) / (0.33*0.33) + (poe * poe) / (3.2* 3.2)<=1)
        Distance cut: r_med_geo<=1500 ||  1000/parallax <= 1500
        Eliminate duplicates both ways (Gaia and GALEX), keep always nearest
            Gaia - GALEX match only.
        Eliminate spurious sources in the Magellanic Clouds:
            !((sqrt( pow(l-280.3,2) + 2*pow(b+33.0,2) ) < 8) ||
            sqrt(0.5*square(l-301) + 2*square(b+43.3) ) < 5)

    """

    name = 'mwm_cb_galex_vol'
    mapper = 'MWM'
    category = 'science'
    instrument = 'BOSS'
    cadence = None
    program = 'mwm_cb'
    priority = 1400
    can_offset = True

    def build_query(self, version_id, query_region=None):

        G3 = Gaia_DR3
        pm = fn.sqrt(G3.pmra * G3.pmra + G3.pmdec * G3.pmdec)
        pm_error = fn.sqrt((G3.pmra_error * G3.pmra / pm) * (G3.pmra_error * G3.pmra / pm) +
                           (G3.pmdec_error * G3.pmdec / pm) * (G3.pmdec_error * G3.pmdec / pm))
        logpmdpm = fn.log(pm / pm_error)
        poe = G3.parallax / G3.parallax_error
        logpar = fn.log(fn.abs(1000 / G3.parallax))

        nuv_mag = GUVCat.nuv_mag
        nuv_magerr = GUVCat.nuv_magerr
        Gm = G3.phot_g_mean_mag
        bp = G3.phot_bp_mean_mag
        rp = G3.phot_rp_mean_mag
        parallax = G3.parallax
        ll = G3.l
        bb = G3.b

        r_med_geo = BailerJonesEDR3.r_med_geo

        nuv_cut = ((((nuv_mag - 5 * fn.log(r_med_geo) + 5) > (2 * (nuv_mag - Gm) + 0.7)) |
                    ((nuv_mag - 5 * logpar + 5) > (2 * (nuv_mag - Gm) + 0.7))) &
                   ((nuv_mag - Gm) < (1.6 * (bp - rp) + 1.2)))

        astro_cut = ~(((logpmdpm < 0.301) &
                      (poe > (-1.75 * logpmdpm - 2.8)) &
                      (poe < (1.75 * logpmdpm + 2.8))) |
                      (((logpmdpm - 0.301) * (logpmdpm - 0.301) / (0.33 * 0.33) +
                        (poe * poe) / (3.2 * 3.2)) <= 1))

        magellanic_cut = ~((fn.sqrt(fn.pow(ll - 280.3, 2) + 2 * fn.pow(bb + 33.0, 2)) < 8) |
                           (fn.sqrt(0.5 * fn.pow(ll - 301, 2) + 2 * fn.pow(bb + 43.3, 2)) < 5))

        cadence = peewee.Case(None, ((Gm <= 16, 'bright_2x1'), (Gm > 16, 'dark_2x1')))

        q1 = (Gaia_DR3
              .select(CatalogToGaia_DR3.catalogid,
                      G3.source_id,
                      G3.ra,
                      G3.dec,
                      cadence.alias('cadence'))
              .join(Galex_GR7_Gaia_DR3,
                    on=(Gaia_DR3.source_id == Galex_GR7_Gaia_DR3.gaia_edr3_source_id))
              .join(GUVCat,
                    on=(GUVCat.objid == Galex_GR7_Gaia_DR3.galex_objid))
              .join_from(Gaia_DR3, BailerJonesEDR3,
                         on=(Gaia_DR3.source_id == BailerJonesEDR3.source_id))
              .join_from(Gaia_DR3, CatalogToGaia_DR3)
              .join_from(CatalogToGaia_DR3, CatalogToMILLIQUAS_7_7,
                         join_type=peewee.JOIN.LEFT_OUTER,
                         on=(CatalogToMILLIQUAS_7_7.catalogid == CatalogToGaia_DR3.catalogid))
              .join(MILLIQUAS_7_7, join_type=peewee.JOIN.LEFT_OUTER)
              .where(CatalogToGaia_DR3.version_id == version_id,
                     CatalogToGaia_DR3.best >> True,
                     Galex_GR7_Gaia_DR3.galex_separation < 2.5,
                     (CatalogToMILLIQUAS_7_7.catalogid.is_null() |
                      ((CatalogToMILLIQUAS_7_7.best >> True) & (MILLIQUAS_7_7.qpct != 100))),
                     nuv_mag > -999,
                     nuv_magerr < 0.2,
                     nuv_magerr > 0,
                     nuv_cut,
                     astro_cut,
                     (r_med_geo <= 1500) | (fn.abs(1000 / parallax) <= 1500),
                     magellanic_cut)
              .order_by(CatalogToGaia_DR3.catalogid, Galex_GR7_Gaia_DR3.galex_separation)
              .distinct(CatalogToGaia_DR3.catalogid))

        q1.create_table('mwm_cb_galex_vol_temp', temporary=True)

        temp = peewee.Table('mwm_cb_galex_vol_temp').bind(self.database)
        t2 = (temp
              .select(temp.c.catalogid.alias('c2'))
              .join(HECATE_1_1,
                    on=(fn.q3c_join(HECATE_1_1.ra, HECATE_1_1.dec,
                                    temp.c.ra, temp.c.dec,
                                    HECATE_1_1.r1 / 60.) &
                        (HECATE_1_1.r1 != peewee.SQL("'NaN'::numeric"))))
              .distinct(temp.c.catalogid)
              .alias('hecate_q3c_join'))

        query = (temp
                 .select(temp.c.catalogid,
                         temp.c.ra,
                         temp.c.dec,
                         temp.c.priority)
                 .join(t2,
                       join_type=peewee.JOIN.LEFT_OUTER,
                       on=(temp.c.catalogid == t2.c.c2))
                 .where(t2.c.c2.is_null()))

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_CB_XMMOM(BaseCarton):
    """Volume-limited XMM-OM / Gaia selection.

    Simplified Description of selection criteria:
        Will use matching XMM-OM / Gaia DR3 sources and slightly revised color cuts compared to
        V0.5. Select objects that lie below the main-sequence in an NUV/optical absolute
        magnitude-color diagram. Remove known AGN & galaxies from published catalogues. Remove AGN
        interlopers by parallax / proper motion. Remove the LMC and SMC, where a lot of hot young
        stars mimic the properties of this target selection. Makes use of milliquas (million quasar
        catalog https://quasars.org/milliquas.htm) and HECATE (Heraklion Extragalactic Catalog -
        https://hecate.ia.forth.gr).

    Pseudo-code:
        * Preparation:
            Trim milliquas catalog
                java -jar stilts.jar tpipe in=milliquas.fits cmd='select QPCT==100'
                    out=milliquas_good1.fits
            Reduce file size of XMMOM input by selecting only used columns and by applying
            quality cuts:
                java -jar stilts.jar tpipe in=XMM-OM-SUSS5.0.fits out=xmm_reduced1.fits
                    cmd="select UVM2_AB_MAG>=-100" cmd="keepcols 'SRCNUM OBSID RA DEC UVW1_AB_MAG
                    UVW2_AB_MAG UVM2_AB_MAG UVW1_QUALITY_FLAG_ST UVW2_QUALITY_FLAG_ST
                    UVM2_QUALITY_FLAG_ST UVW1_SIGNIF UVW2_SIGNIF UVM2_SIGNIF'"
        * Gaia match within 1.5" and Bailer-Jones distances
            java -jar stilts.jar cdsskymatch cdstable=I/350/gaiaedr3 find=all
                icmd='select (!(equals(substring(UVM2_QUALITY_FLAG_ST,0,1),\"T\")||
                        (equals(substring(UVM2_QUALITY_FLAG_ST,1,2),\"T\")&&UVM2_SIGNIF<10)||
                        (equals(substring(UVM2_QUALITY_FLAG_ST,2,3),\"T\")&&UVM2_SIGNIF<10)||
                        equals(substring(UVM2_QUALITY_FLAG_ST,6,7),\"T\")||
                        equals(substring(UVM2_QUALITY_FLAG_ST,7,8),\"T\")||
                        equals(substring(UVM2_QUALITY_FLAG_ST,8,9),\"T\")))'
                    in=xmm_reduced1.fits ra= dec=DEC radius=1.5 blocksize=50000
                    out=xmm_gaia.fits
            java -jar stilts.jar cdsskymatch cdstable=I/352/gedr3dis find=best
                in=xmm_gaia.fits ra=ra_cds dec=dec_cds radius=0.001
                blocksize=50000 out=xmm_gaia_dist.fits
        * Eliminate known AGN, i.e. all milliquas entries with 100% probability being an AGN
          and eliminate UV sources in nearby galaxies
            java -jar -Xmx16G stilts.jar tmatch2 in1=xmm_gaia_dist.fits in2=milliquas_good1.fits
                matcher=sky params=1.0 values1='ra_cds dec_cds' values2='RA DEC'
                join=1not2 find=all out=xmm_gaia_dist_noAGN.fits
            java -jar -Xmx16G stilts.jar tmatch2 in1=xmm_gaia_dist_noAGN.fits in2=HECATE_v1.1.fits
                matcher=skyerr params=10 values1='ra_cds dec_cds 0.0001' values2='RA DEC 60*R1'
                join=1not2 find=all out=reduced_XMM_all_data.fits
        * Select entries with maximum separation between Gaia and XMM-OM position of 1.5 arcsec
        * Select objects with an UV excess via two cuts, one in the in the NUV-HRD,
          the other in an UV-to-optical color-color diagram
          - UVM2_AB_MAG>-999 && UVM2_AB_MAG-5*log10(rgeo)+5>2*(UVM2_AB_MAG-phot_g_mean_mag)+0.7 ||
            (UVM2_AB_MAG-5*log10(parallax/1000)+5>2*(UVM2_AB_MAG-phot_g_mean_mag)+0.7)
          - nuv_mag - phot_g_mean_mag < 1.6 * (bp-rp) + 1.2
        * Astrometric cut: to remove AGN contaminants:
            !(logpmdpm < 0.301 && poe > -1.75*logpmdpm - 2.8 && poe <  1.75*logpmdpm + 2.8 ||
            (logpmdpm-0.301)*(logpmdpm-0.301) / (0.33*0.33) + (poe * poe) / (3.2* 3.2)<=1)
        * Distance cut: r_med_geo<=1500 ||  1000/parallax <= 1500
        * Eliminate duplicates both ways (Gaia and XMMOM), keep always nearest
          Gaia - XMMOM match only.

    """

    name = 'mwm_cb_xmmom'
    mapper = 'MWM'
    category = 'science'
    instrument = 'BOSS'
    cadence = None
    program = 'mwm_cb'
    priority = 1400
    can_offset = True

    def build_query(self, version_id, query_region=None):

        G3 = Gaia_DR3
        pm = fn.sqrt(G3.pmra * G3.pmra + G3.pmdec * G3.pmdec)
        pm_error = fn.sqrt((G3.pmra_error * G3.pmra / pm) * (G3.pmra_error * G3.pmra / pm) +
                           (G3.pmdec_error * G3.pmdec / pm) * (G3.pmdec_error * G3.pmdec / pm))
        logpmdpm = fn.log(pm / pm_error)
        poe = G3.parallax / G3.parallax_error

        Gm = G3.phot_g_mean_mag
        parallax = G3.parallax
        logpar = fn.log(fn.abs(parallax / 1000))

        r_med_geo = BailerJonesEDR3.r_med_geo

        XMM = XMM_OM_SUSS_5_0
        uvm2 = XMM.uvm2_ab_mag

        qcut = ~((fn.substring(XMM.uvm2_quality_flag_st, 0, 1) == 'T') |
                 ((fn.substring(XMM.uvm2_quality_flag_st, 1, 2) == 'T') & (XMM.uvm2_signif < 10)) |
                 ((fn.substring(XMM.uvm2_quality_flag_st, 2, 3) == 'T') & (XMM.uvm2_signif < 10)) |
                 (fn.substring(XMM.uvm2_quality_flag_st, 6, 7) == 'T') |
                 (fn.substring(XMM.uvm2_quality_flag_st, 7, 8) == 'T') |
                 (fn.substring(XMM.uvm2_quality_flag_st, 8, 9) == 'T'))

        color_cuts = (((uvm2 > -999) &
                       ((uvm2 - 5 * fn.log10(r_med_geo) + 5) > (2 * (uvm2 - Gm) + 0.7)) |
                       ((uvm2 - 5 * logpar + 5) > (2 * (uvm2 - Gm) + 0.7))))

        astro_cut = ~(((logpmdpm < 0.301) &
                      (poe > (-1.75 * logpmdpm - 2.8)) &
                      (poe < (1.75 * logpmdpm + 2.8))) |
                      (((logpmdpm - 0.301) * (logpmdpm - 0.301) / (0.33 * 0.33) +
                        (poe * poe) / (3.2 * 3.2)) <= 1))

        cadence = peewee.Case(None, ((Gm <= 16, 'bright_2x1'), (Gm > 16, 'dark_2x1')))

        q1 = (Gaia_DR3
              .select(CatalogToGaia_DR3.catalogid,
                      G3.source_id,
                      G3.ra,
                      G3.dec,
                      cadence.alias('cadence'))
              .join_from(Gaia_DR3, BailerJonesEDR3,
                         on=(Gaia_DR3.source_id == BailerJonesEDR3.source_id))
              .join_from(Gaia_DR3, CatalogToGaia_DR3)
              .join_from(CatalogToGaia_DR3, CatalogToMILLIQUAS_7_7,
                         join_type=peewee.JOIN.LEFT_OUTER,
                         on=(CatalogToMILLIQUAS_7_7.catalogid == CatalogToGaia_DR3.catalogid))
              .join(MILLIQUAS_7_7, join_type=peewee.JOIN.LEFT_OUTER)
              .join_from(CatalogToGaia_DR3, CatalogToXMM_OM_SUSS_5_0,
                         on=(CatalogToGaia_DR3.catalogid == CatalogToXMM_OM_SUSS_5_0.catalogid))
              .join(XMM_OM_SUSS_5_0)
              .where(CatalogToGaia_DR3.version_id == version_id,
                     CatalogToGaia_DR3.best >> True,
                     CatalogToXMM_OM_SUSS_5_0.best >> True,
                     (CatalogToMILLIQUAS_7_7.catalogid.is_null() |
                         ((CatalogToMILLIQUAS_7_7.best >> True) & (MILLIQUAS_7_7.qpct != 100))),
                     (r_med_geo <= 1500) | (fn.abs(1000 / parallax) <= 1500),
                     uvm2 >= -100,
                     qcut,
                     color_cuts,
                     astro_cut)
              .order_by(CatalogToGaia_DR3.catalogid, CatalogToXMM_OM_SUSS_5_0.distance)
              .distinct(CatalogToGaia_DR3.catalogid))

        q1.create_table('mwm_cb_xmmom_temp', temporary=True)

        temp = peewee.Table('mwm_cb_xmmom_temp').bind(self.database)
        t2 = (temp
              .select(temp.c.catalogid.alias('c2'))
              .join(HECATE_1_1,
                    on=(fn.q3c_join(HECATE_1_1.ra, HECATE_1_1.dec,
                                    temp.c.ra, temp.c.dec,
                                    HECATE_1_1.r1 / 60.) &
                        (HECATE_1_1.r1 != peewee.SQL("'NaN'::numeric"))))
              .distinct(temp.c.catalogid)
              .alias('hecate_q3c_join'))

        query = (temp
                 .select(temp.c.catalogid,
                         temp.c.ra,
                         temp.c.dec,
                         temp.c.priority)
                 .join(t2,
                       join_type=peewee.JOIN.LEFT_OUTER,
                       on=(temp.c.catalogid == t2.c.c2))
                 .where(t2.c.c2.is_null()))

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_CB_SWIFTUVOT_Carton(BaseCarton):
    """MWM Compact Binaries UV excess 4.

    Simplified Description of selection criteria:
        Will use matching Swift-UVOT / Gaia DR3 sources and slightly revised color cuts compared to
        V0.5. Select objects that lie below the main-sequence in an NUV/optical absolute
        magnitude-color diagram. Remove known AGN & galaxies from published catalogues. Remove AGN
        interlopers by parallax / proper motion. Remove the LMC and SMC, where a lot of hot young
        stars mimic the properties of this target selection. Makes use of milliquas (million quasar
        catalog https://quasars.org/milliquas.htm) and HECATE (Heraklion Extragalactic Catalog -
        https://hecate.ia.forth.gr)  ← ingestion 'done' in the catalog ingestion list on the
        timeline + tasks for v1 page.  Full details, see here.

    Pseudo-code:
        - Target selection analog to mwm_cb_xmmom. Names of input columns in swiftuvot catalog
          slightly different:
        - Quality cuts in input catalog:
            select (!(equals(substring(UVM2_QUALITY_FLAG_BIN,8,9),"1") ||
            (equals(substring(UVM2_QUALITY_FLAG_BIN,7,8),"1") && UVM2_SIGNIF<10) ||
            (equals(substring(UVM2_QUALITY_FLAG_BIN,6,7),"1") && UVM2_SIGNIF<10) ||
             equals(substring(UVM2_QUALITY_FLAG_BIN,3,4),"1") && UVM2_SIGNIF<10) ||
             equals(substring(UVM2_QUALITY_FLAG_BIN,2,3),"1"))"'
        - Gaia match within 1.5" and Bailer-Jones distances (r_med_geo)
        - Eliminate known AGN, i.e. all milliquas entries with 100% probability
          being an AGN and eliminate UV sources in nearby galaxies
            java -jar -Xmx16G stilts.jar tmatch2 in1=xmm_gaia_dist.fits in2=milliquas_good1.fits
              matcher=sky params=1.0 values1='ra_cds dec_cds' values2='RA DEC' join=1not2
              find=all out=swiftuvot_gaia_dist_noAGN.fits
            java -jar -Xmx16G stilts.jar tmatch2 in1=swiftuvot_gaia_dist_noAGN.fits
              in2=HECATE_v1.1.fits matcher=skyerr params=10 values1='ra_cds dec_cds 0.0001'
              values2='RA DEC 60*R1' join=1not2 find=all out=reduced_swift_all_data.fits
        - Select objects with an UV excess via two cuts, one in the in the NUV-HRD, the other
          in an UV-to-optical color-color diagram:
            UVM2_ABMAG>-999 && UVM2_ABMAG-5*log10(r_med_geo)+5>2*(UVM2_ABMAG-phot_g_mean_mag)+0.7
                || (UVM2_ABMAG-5*log10(parallax/1000)+5>2*(UVM2_ABMAG-phot_g_mean_mag)+0.7)
            UVM2_ABMAG - phot_g_mean_mag < 1.6 * (bp-rp) + 1.2
        - Astrometric cut: to remove AGN contaminants:
            !(logpmdpm < 0.301 && poe > -1.75*logpmdpm - 2.8 && poe <  1.75*logpmdpm + 2.8 ||
            (logpmdpm-0.301)*(logpmdpm-0.301) / (0.33*0.33) + (poe * poe) / (3.2* 3.2)<=1)
        - Distance cut: r_med_geo<=1500 ||  1000/parallax <= 1500
        - Eliminate duplicates both ways (Gaia and Swift), keep always nearest
          Gaia - Swiftuvot match only.

    """

    name = 'mwm_cb_swiftuvot'
    mapper = 'MWM'
    category = 'science'
    program = 'mwm_cb'
    instrument = 'BOSS'
    cadence = None
    priority = 1400

    def build_query(self, version_id, query_region=None):

        G3 = Gaia_DR3
        pm = fn.sqrt(G3.pmra * G3.pmra + G3.pmdec * G3.pmdec)
        pm_error = fn.sqrt((G3.pmra_error * G3.pmra / pm) * (G3.pmra_error * G3.pmra / pm) +
                           (G3.pmdec_error * G3.pmdec / pm) * (G3.pmdec_error * G3.pmdec / pm))
        logpmdpm = fn.log(pm / pm_error)
        poe = G3.parallax / G3.parallax_error

        Gm = G3.phot_g_mean_mag
        parallax = G3.parallax
        logpar = fn.log(fn.abs(parallax / 1000))

        r_med_geo = BailerJonesEDR3.r_med_geo

        uvw2_ab_mag = UVOT_SSC_1.uvw2_ab
        uvm2 = UVOT_SSC_1.uvm2_ab
        uvw1_ab_mag = UVOT_SSC_1.uvw1_ab

        uvm2_quality_flag = UVOT_SSC_1.fuvm2
        uvm2_signif = UVOT_SSC_1.suvm2

        qcut = (
            (uvw2_ab_mag > -100) | (uvw1_ab_mag > -100) | (uvm2 > -100),
            ~((uvm2_quality_flag.bin_and(2**0) > 0) |
              (uvm2_quality_flag.bin_and(2**1) > 0) & (uvm2_signif < 10) |
              (uvm2_quality_flag.bin_and(2**2) > 0) & (uvm2_signif < 10) |
              (uvm2_quality_flag.bin_and(2**5) > 0) & (uvm2_signif < 10) |
              (uvm2_quality_flag.bin_and(2**6) > 0))
        )

        color_cuts = (((uvm2 > -999) &
                       ((uvm2 - 5 * fn.log10(r_med_geo) + 5) > (2 * (uvm2 - Gm) + 0.7)) |
                       ((uvm2 - 5 * logpar + 5) > (2 * (uvm2 - Gm) + 0.7))))

        astro_cut = ~(((logpmdpm < 0.301) &
                      (poe > (-1.75 * logpmdpm - 2.8)) &
                      (poe < (1.75 * logpmdpm + 2.8))) |
                      (((logpmdpm - 0.301) * (logpmdpm - 0.301) / (0.33 * 0.33) +
                        (poe * poe) / (3.2 * 3.2)) <= 1))

        cadence = peewee.Case(None, ((Gm <= 16, 'bright_2x1'), (Gm > 16, 'dark_2x1')))

        q1 = (Gaia_DR3
              .select(CatalogToGaia_DR3.catalogid,
                      G3.source_id,
                      G3.ra,
                      G3.dec,
                      cadence.alias('cadence'))
              .join_from(Gaia_DR3, BailerJonesEDR3,
                         on=(Gaia_DR3.source_id == BailerJonesEDR3.source_id))
              .join_from(Gaia_DR3, CatalogToGaia_DR3)
              .join_from(CatalogToGaia_DR3, CatalogToMILLIQUAS_7_7,
                         join_type=peewee.JOIN.LEFT_OUTER,
                         on=(CatalogToMILLIQUAS_7_7.catalogid == CatalogToGaia_DR3.catalogid))
              .join(MILLIQUAS_7_7, join_type=peewee.JOIN.LEFT_OUTER)
              .join_from(CatalogToGaia_DR3, CatalogToUVOT_SSC_1,
                         on=(CatalogToGaia_DR3.catalogid == CatalogToUVOT_SSC_1.catalogid))
              .join(UVOT_SSC_1)
              .where(CatalogToGaia_DR3.version_id == version_id,
                     CatalogToGaia_DR3.best >> True,
                     CatalogToUVOT_SSC_1.best >> True,
                     (CatalogToMILLIQUAS_7_7.catalogid.is_null() |
                         ((CatalogToMILLIQUAS_7_7.best >> True) & (MILLIQUAS_7_7.qpct != 100))),
                     (r_med_geo <= 1500) | (fn.abs(1000 / parallax) <= 1500),
                     uvm2 >= -100,
                     *qcut,
                     color_cuts,
                     astro_cut)
              .order_by(CatalogToGaia_DR3.catalogid, CatalogToUVOT_SSC_1.distance)
              .distinct(CatalogToGaia_DR3.catalogid))

        q1.create_table('mwm_cb_swiftuvot_temp', temporary=True)

        temp = peewee.Table('mwm_cb_swiftuvot_temp').bind(self.database)
        t2 = (temp
              .select(temp.c.catalogid.alias('c2'))
              .join(HECATE_1_1,
                    on=(fn.q3c_join(HECATE_1_1.ra, HECATE_1_1.dec,
                                    temp.c.ra, temp.c.dec,
                                    HECATE_1_1.r1 / 60.) &
                        (HECATE_1_1.r1 != peewee.SQL("'NaN'::numeric"))))
              .distinct(temp.c.catalogid)
              .alias('hecate_q3c_join'))

        query = (temp
                 .select(temp.c.catalogid,
                         temp.c.ra,
                         temp.c.dec,
                         temp.c.priority)
                 .join(t2,
                       join_type=peewee.JOIN.LEFT_OUTER,
                       on=(temp.c.catalogid == t2.c.c2))
                 .where(t2.c.c2.is_null()))

        if query_region:
            query = (query
                     .join_from(CatalogToGaia_DR3, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
