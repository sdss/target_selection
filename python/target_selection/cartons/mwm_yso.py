#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Pramod Gupta (psgupta@uw.edu)
# @Date: 2020-06-10
# @Filename: mwm_yso.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee

from sdssdb.peewee.sdss5db.catalogdb import (MIPSGAL, AllWise, Catalog,
                                             CatalogToTIC_v8, Gaia_DR2, TIC_v8,
                                             TwoMassPSC, YSO_Clustering)

from target_selection.cartons import BaseCarton


# See catalog.py for the name of peewee model names corresponding
# to postgres table names:
# https://github.com/sdss/sdssdb/blob/master/python/sdssdb/peewee/sdss5db/catalogdb.py


class MWM_YSO_S1_Carton(BaseCarton):
    """YSOs - S1 (IR excess).

    Shorthand name: mwm_yso_s1
    Simplified Description of selection criteria:
    selection of YSOs based on IR excess,
    with WISE colors W1-W2>0.25, W2-W3>0.5, W3-W4>1.5,
    closer than parallax>0.3, and brighter than H<13
    (should have ~21.5K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: Gaia, 2mass, allwise
    Additional cross-matching needed:
    Note: Using the Gaia xmatch somehow misses half the sources.
    Selection was done on the allwise catalog that
    had 2mass photometry,
    and then the resulting selection was crossmatched against against
    Gaia with 1" search radius.
    Return columns: Gaia id, 2mass id, allwise id, G, BP, RP,
    J, H, K, W1, W2, W3, W4,parallax
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation: h_m<13 and w1mpro-w2mpro>0.25 and
    w2mpro-w3mpro>0.5 and w3mpro-w4mpro>1.5 and parallax>0.3

    """

    name = 'mwm_yso_s1'
    category = 'science'
    cadence = None
    program = 'YSO'
    mapper = 'MWM'
    priority = 2700

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid)
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .switch(TIC_v8)
                 .join(AllWise, on=(TIC_v8.allwise == AllWise.designation))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        TwoMassPSC.h_m < 13,
                        (AllWise.w1mpro - AllWise.w2mpro) > 0.25,
                        (AllWise.w2mpro - AllWise.w3mpro) > 0.50,
                        (AllWise.w3mpro - AllWise.w4mpro) > 1.50,
                        Gaia_DR2.parallax > 0.3))

        # Gaia_DR2 pweewee model class corresponds to
        # table catalogdb.gaia_dr2_source.
        #
        # All values of TIC_v8.plx (for non-null entries) are not the same as
        # values of Gaia_DR2.parallax.
        # Hence, in the above query, we cannot use TIC_v8.plx instead
        # of Gaia_DR2.parallax.

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_YSO_S2_Carton(BaseCarton):
    """YSOs - S2 (optically invisible).

    Shorthand name: mwm_yso_s2
    Simplified Description of selection criteria:
    selection of YSOs, brighter than H<13, fainter than G>15 or
    without gaia detection,
    colors J-H>0,5, W1-W2>0.5, W2-W3>1, W3-W4>1.5, and
     relates (W3-W4)>(W1-W2)*0.5+1.1
    (should have ~11.6K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: 2mass+allwise, gaia
    (allow sources that lack gaia xmatch)
    Additional cross-matching needed:
    Note: Using the Gaia xmatch somehow misses half the sources.
    Selection was done on the allwise catalog
    that had 2mass photometry,
    and then the resulting selection was crossmatched
    against against Gaia with 1" search radius.
    Return columns: Gaia id, 2mass id, allwise id, G, BP, RP,
    J, H, K, W1, W2, W3, W4
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation: h_m<13 and
    (phot_g_mean_mag>18.5 or phot_g_mean_mag is null)
    and j_m-h_m>1
    and h_m-ks_m>0.5
    and w1mpro-w2mpro>0.5
    and w2mpro-w3mpro>1
    and w3mpro-w4mpro>1.5
    and w3mpro-w4mpro>(w1mpro-w2mpro)*0.8+1.1

    """

    name = 'mwm_yso_s2'
    category = 'science'
    cadence = None
    program = 'YSO'
    mapper = 'MWM'
    priority = 2700

    def build_query(self, version_id, query_region=None):

        query = (AllWise
                 .select(CatalogToTIC_v8.catalogid)
                 .join(TIC_v8, on=(TIC_v8.allwise == AllWise.designation))
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .switch(TIC_v8)
                 .join(Gaia_DR2, peewee.JOIN.LEFT_OUTER,
                       on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(CatalogToTIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        TwoMassPSC.h_m < 13,
                        (Gaia_DR2.phot_g_mean_mag > 18.5) |
                        (Gaia_DR2.phot_g_mean_mag >> None),
                        (AllWise.j_m_2mass - AllWise.h_m_2mass) > 1.0,
                        (AllWise.h_m_2mass - AllWise.k_m_2mass) > 0.5,
                        (AllWise.w1mpro - AllWise.w2mpro) > 0.50,
                        (AllWise.w2mpro - AllWise.w3mpro) > 1.00,
                        (AllWise.w3mpro - AllWise.w4mpro) > 1.50,
                        (AllWise.w3mpro - AllWise.w4mpro) >
                        (AllWise.w1mpro - AllWise.w2mpro) * 0.8 + 1.1))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
        return query


class MWM_YSO_S2_5_Carton(BaseCarton):
    """YSOs - S2.5 (optically invisible, WISE saturated).

    Shorthand name: mwm_yso_s2_5
    Simplified Description of selection criteria:
    selection of YSOs, brighter than H<15,
    saturated (blank) W4 with W2-W3>4,
    or saturated W3 and W2, with J-H>1.1.
    Some contaminants from scanning are
    filtered on the plane of the sky:
    all the targets should be within 5 deg of the plane+
    few sources that can be located
    further south of the plane if l>180
    (should have ~1.2K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: 2mass, allwise
    Additional cross-matching needed:
    Return columns: 2mass id, allwise id, J, H, K, W1, W2, W3, W4
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation: h_m<13 and
    (w2mpro-w3mpro>4 and w4mpro is null) or
    (w3mpro is null and w4mpro is null and j_m-h_m>1.1)
    and (b>-5 or l>180) and b<-5

    """

    name = 'mwm_yso_s2_5'
    category = 'science'
    cadence = None
    program = 'YSO'
    mapper = 'MWM'
    priority = 2700

    # Above implementation has below clause
    # and (b>-5 or l>180) and b<-5
    # Replace (b>-5 or l>180) and b<-5 as below based on the text.
    # In words:
    # all the targets should be within 5 deg of the plane+
    # few sources that can be
    # located further south of the plane if l>180
    # Hence:
    # ((b>-5) and (b<5)) or ((b<-5) and (l > 180))
    #  l, b in Gaia_DR2 are gallong and gallat in TIC_v8.
    # We are using the values from Gaia since
    # TIC propagates the coordinates back to epoch 2000.0
    # (b>-5 or l>180) and b<-5
    # S2_5 query below has the same part before where() as S2 query.

    def build_query(self, version_id, query_region=None):

        query = (AllWise
                 .select(CatalogToTIC_v8.catalogid)
                 .join(TIC_v8, on=(TIC_v8.allwise == AllWise.designation))
                 .join(TwoMassPSC,
                       on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .switch(TIC_v8)
                 .join(Gaia_DR2, peewee.JOIN.LEFT_OUTER,
                       on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .switch(TIC_v8)
                 .join(CatalogToTIC_v8,
                       on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        TwoMassPSC.h_m < 13,
                        (((AllWise.w2mpro - AllWise.w3mpro) > 4) &
                         (AllWise.w4mpro >> None)) |
                        ((AllWise.w3mpro >> None) &
                         (AllWise.w4mpro >> None) &
                         ((AllWise.j_m_2mass - AllWise.h_m_2mass) > 1.1)),
                        ((Gaia_DR2.b > -5) & (Gaia_DR2.b < 5)) |
                        ((Gaia_DR2.b < -5) & (Gaia_DR2.l > 180)) |
                        ((Gaia_DR2.b >> None) & (Gaia_DR2.l >> None))))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
        return query


class MWM_YSO_S3_Carton(BaseCarton):
    """YSOs - S3 (pre-main sequence optical variables).

    Shorthand name: mwm_yso_s3
    Simplified Description of selection criteria:
    selection of YSOs brighter than H<13, closer than parallax>0.3.
    Filter on the position of the HR diagram to
    select cool pre-main sequence stars,
    with BP-RP>13, (BP-RP)*2.5+2.5>M_G, (BP-RP)*2.5-1<M_G,
    requiring variability in g,bp,rp>0.02
    (with var_x defined as
    sqrt(phot_x_n_obs)/phot_x_mean_flux_over_error),
    have relations in variability of
    var_g<var_bp<var_g^0.75, 0.75*var_g<var_rp<var_g^0.95,
    and log10(var_bp)*5+11<M_BP, in which M_x is the absolute mag
    (should have ~52.7K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: 2mass, gaia
    Additional cross-matching needed:
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation:
    phot_g_mean_mag < 18.5 and h_m <13 and parallax >0.3 and
    bp_rp*2.5+2.5 > phot_g_mean_mag-5*(log10(1000/parallax)-1) and
    bp_rp*2.5-1 < phot_g_mean_mag-5*(log10(1000/parallax)-1) and
    sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error>
    sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error and
    sqrt(phot_rp_n_obs)/phot_rp_mean_flux_over_error>
    sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error*0.75 and
    sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error<
    power(sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error,0.75) and
    sqrt(phot_rp_n_obs)/phot_rp_mean_flux_over_error<
    power(sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error,0.95) and
    log10(sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error)*5+11<
    phot_bp_mean_mag-5*(log10(1000/parallax)-1) and
    bp_rp>1.3 and sqrt(phot_g_n_obs)/phot_g_mean_flux_over_error>0.02
    and
    sqrt(phot_bp_n_obs)/phot_bp_mean_flux_over_error>0.02 and
    sqrt(phot_rp_n_obs)/phot_rp_mean_flux_over_error>0.02

    """

    name = 'mwm_yso_s3'
    category = 'science'
    cadence = None
    program = 'YSO'
    mapper = 'MWM'
    priority = 2700

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid)
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .switch(TIC_v8)
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        Gaia_DR2.phot_g_mean_mag < 18.5,
                        TwoMassPSC.h_m < 13,
                        Gaia_DR2.parallax > 0.3,
                        Gaia_DR2.bp_rp * 2.5 + 2.5 >
                        Gaia_DR2.phot_g_mean_mag -
                        5 * (peewee.fn.log(1000 / Gaia_DR2.parallax) - 1),
                        Gaia_DR2.bp_rp * 2.5 - 1 <
                        Gaia_DR2.phot_g_mean_mag -
                        5 * (peewee.fn.log(1000 / Gaia_DR2.parallax) - 1),
                        peewee.fn.sqrt(Gaia_DR2.phot_bp_n_obs) /
                        Gaia_DR2.phot_bp_mean_flux_over_error >
                        peewee.fn.sqrt(Gaia_DR2.phot_g_n_obs) /
                        Gaia_DR2.phot_g_mean_flux_over_error,
                        peewee.fn.sqrt(Gaia_DR2.phot_rp_n_obs) /
                        Gaia_DR2.phot_rp_mean_flux_over_error >
                        peewee.fn.sqrt(Gaia_DR2.phot_g_n_obs) /
                        Gaia_DR2.phot_g_mean_flux_over_error * 0.75,
                        peewee.fn.sqrt(Gaia_DR2.phot_bp_n_obs) /
                        Gaia_DR2.phot_bp_mean_flux_over_error <
                        peewee.fn.power(
                            peewee.fn.sqrt(Gaia_DR2.phot_g_n_obs) /
                            Gaia_DR2.phot_g_mean_flux_over_error, 0.75),
                        peewee.fn.sqrt(Gaia_DR2.phot_rp_n_obs) /
                        Gaia_DR2.phot_rp_mean_flux_over_error <
                        peewee.fn.power(
                            peewee.fn.sqrt(Gaia_DR2.phot_g_n_obs) /
                            Gaia_DR2.phot_g_mean_flux_over_error, 0.95),
                        peewee.fn.log(
                            peewee.fn.sqrt(Gaia_DR2.phot_bp_n_obs) /
                            Gaia_DR2.phot_bp_mean_flux_over_error) * 5 + 11 <
                        Gaia_DR2.phot_bp_mean_mag -
                        5 * (peewee.fn.log(1000 / Gaia_DR2.parallax) - 1),
                        Gaia_DR2.bp_rp > 1.3,
                        peewee.fn.sqrt(Gaia_DR2.phot_g_n_obs) /
                        Gaia_DR2.phot_g_mean_flux_over_error > 0.02,
                        peewee.fn.sqrt(Gaia_DR2.phot_bp_n_obs) /
                        Gaia_DR2.phot_bp_mean_flux_over_error > 0.02,
                        peewee.fn.sqrt(Gaia_DR2.phot_rp_n_obs) /
                        Gaia_DR2.phot_rp_mean_flux_over_error > 0.02))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))
        return query


class MWM_YSO_OB_Carton(BaseCarton):
    """YSOs - Upper (pre-)Main Sequence.

    Shorthand name: mwm_yso_ob
    Simplified Description of selection criteria:
    Selecting the OB stars at the tip of the main sequence,
    brighter than H<13, G<18 mag, closer than parallax>0.3,
    color -0.2<BP-RP<1.1, and M_G<(BP-RP)*1.6-2.2
    (should have ~8.7K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: 2mass, gaia
    Additional cross-matching needed:
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation: h_m<13 and bp_rp between -0.2 and 1.1 and
    phot_g_mean_mag<18 and
    phot_g_mean_mag-5*(log10(1000/parallax)-1) <
    1.6*bp_rp-2.2 and parallax>0.3

    """

    name = 'mwm_yso_ob'
    category = 'science'
    cadence = None
    program = 'YSO'
    mapper = 'MWM'
    priority = 2700

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid)
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(TwoMassPSC, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .switch(TIC_v8)
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        TwoMassPSC.h_m < 13,
                        (Gaia_DR2.bp_rp > -0.2) & (Gaia_DR2.bp_rp < 1.1),
                        Gaia_DR2.phot_g_mean_mag < 18,
                        Gaia_DR2.phot_g_mean_mag -
                        5 * (peewee.fn.log(1000 / Gaia_DR2.parallax) - 1) <
                        1.6 * Gaia_DR2.bp_rp - 2.2,
                        Gaia_DR2.parallax > 0.3))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_YSO_CMZ_Carton(BaseCarton):
    """YSOs - Central Molecular Zone.

    Shorthand name: mwm_yso_cmz
    Simplified Description of selection criteria:
    selection of sources in the central molecular zone
    based on spitzer fluxes from mipsgal.
    Sources are within 2 degrees in l and
     1 degree in b from the galactic center,
     brighter than H<13, have color 8.0-24>2.5, and
     have parallax<0.2 or lack a Gaia xmatch.
    (should have ~3.2K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: mipsgal
    Additional cross-matching needed: the table has xmatch included
    Return columns:
    mipsgal id, 2mass id, j, h, k, 3.6, 4.8, 8.0, 24 mag
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation: Hmag<13 and (l> 358 or l< 2) and
    b between -1 and 1 and _8_0_-_24_>2.5 and
    (parallax<0.2 or parallax is null)

    For CMZ, the raw sql query would be:
    select ct.catalogid from mipsgal m
    join twomass_psc t on twomass_name = designation
    join tic_v8 tic on tic.twomass_psc = t.designation
    left outer join gaia_dr2_source g on g.source_id = tic.gaia_int
    join catalog_to_tic_v8 ct on ct.target_id = tic.id
    where m.hmag < 13 and
    (m.glon > 358 or m.glon < 2) and
    (m.glat > -1 and m.glat < 1) and
    (m.mag_8_0 - m.mag_24) > 2.5 and
    (g.parallax < 0.2 or g.parallax is null)
    and ct.version_id = 13 and ct.best is true;

    Note you only need one left outer join between TIC and Gaia
    (all MIPSGAL targets have a counterpart in 2MASS,
    and all 2MASS have an entry in TIC,
    but not all the TIC entries have a Gaia counterpart).

    """

    name = 'mwm_yso_cmz'
    category = 'science'
    cadence = None
    program = 'YSO'
    mapper = 'MWM'
    priority = 2700

    # l is glon (galactic longitude)
    # b is glat (galactic latitude)
    # mipsgal is a subset of 2MASS
    # mipsgal can be joined to twomass_psc via
    # mipsgal.twomass_name = TwoMassPSC.designation.
    # Then join via TIC and catalog_to_tic.
    #
    # mipsgal is a subset of 2MASS
    # 2MASS is a subset of TIC_v8
    # Gaia_DR2 is a subset of TIC_v8
    #
    # 2MASS is not a subset of Gaia_DR2
    # Gaia_DR2 is not a subset of 2MASS
    #
    # table catalogdb.mipsgal
    # Foreign-key constraints:
    #    "twomass_name_fk" FOREIGN KEY (twomass_name)
    # REFERENCES twomass_psc(designation)
    #
    # Due to below, we do not need a between to Catalog and CatalogToTIC_v8
    # Catalog.catalogid == CatalogToTIC_v8.catalogid
    # We can remove the join with Catalog in all the cartons
    # since catalogid is completely unique (even across different version_id)
    # so the join with Catalog doesn't give us anything extra and it's a costly join.

    def build_query(self, version_id, query_region=None):

        query = (MIPSGAL.select(CatalogToTIC_v8.catalogid)
                 .join(TwoMassPSC, on=(MIPSGAL.twomass_name == TwoMassPSC.designation))
                 .join(TIC_v8, on=(TIC_v8.twomass_psc == TwoMassPSC.designation))
                 .join(Gaia_DR2, peewee.JOIN.LEFT_OUTER,
                       on=(Gaia_DR2.source_id == TIC_v8.gaia_int))
                 .switch(TIC_v8)
                 .join(CatalogToTIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        MIPSGAL.hmag < 13,
                        (MIPSGAL.glon > 358) | (MIPSGAL.glon < 2),
                        (MIPSGAL.glat > -1) & (MIPSGAL.glat < 1),
                        (MIPSGAL.mag_8_0 - MIPSGAL.mag_24) > 2.5,
                        (Gaia_DR2.parallax < 0.2) |
                        (Gaia_DR2.parallax >> None)))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query


class MWM_YSO_Cluster_Carton(BaseCarton):
    """YSOs - Cluster Catalog
    Shorthand name: mwm_yso_cluster
    Simplified Description of selection criteria:
    Selecting the clustered sources from
    the catalog of clustered structures,
    with age<7.5 dex and brighter than H<13 mag.
    (should have ~45.5K sources)
    Wiki page:
    https://wiki.sdss.org/display/MWM/YSO+selection+function
    Additional source catalogs needed: Kounkel+20 clustered catalog
    Additional cross-matching needed:
    Return columns: Gaia id, 2mass id, G, BP, RP, J, H, K, parallax
    cadence options for these targets
    (list all options,
    even though no single target will receive more than one):
    Pseudo SQL (optional):
    Implementation: age<7.5 and h<13

    """

    name = 'mwm_yso_cluster'
    category = 'science'
    cadence = None
    program = 'YSO'
    mapper = 'MWM'
    priority = 2700

    # yso_clustering is a subset of gaia and
    # can be joined to gaia_dr2_source via source_id.
    #
    # table catalogdb.yso_clustering
    # Foreign-key constraints:
    #    "yso_clustering_source_id_fkey" FOREIGN KEY (source_id)
    # REFERENCES gaia_dr2_source(source_id)

    def build_query(self, version_id, query_region=None):

        query = (CatalogToTIC_v8
                 .select(CatalogToTIC_v8.catalogid)
                 .join(TIC_v8, on=(CatalogToTIC_v8.target_id == TIC_v8.id))
                 .join(Gaia_DR2, on=(TIC_v8.gaia_int == Gaia_DR2.source_id))
                 .join(YSO_Clustering,
                       on=(Gaia_DR2.source_id == YSO_Clustering.source_id))
                 .where(CatalogToTIC_v8.version_id == version_id,
                        CatalogToTIC_v8.best >> True,
                        YSO_Clustering.h < 13,
                        YSO_Clustering.age < 7.5))

        if query_region:
            query = (query
                     .join_from(CatalogToTIC_v8, Catalog)
                     .where(peewee.fn.q3c_radial_query(Catalog.ra,
                                                       Catalog.dec,
                                                       query_region[0],
                                                       query_region[1],
                                                       query_region[2])))

        return query
