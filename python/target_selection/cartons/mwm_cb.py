#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-06-05
# @Filename: mwm_cb.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import peewee
from peewee import fn

from sdssdb.peewee.sdss5db.catalogdb import (
    HECATE_1_1,
    MILLIQUAS_7_7,
    UVOT_SSC_1,
    XMM_OM_SUSS_5_0,
    BailerJonesEDR3,
    Catalog,
    CatalogToGaia_DR3,
    CatalogToMILLIQUAS_7_7,
    CatalogToUVOT_SSC_1,
    CatalogToXMM_OM_SUSS_5_0,
    Gaia_DR3,
    Galex_GR7_Gaia_DR3,
    GUVCat,
)

from target_selection import log
from target_selection.cartons import BaseCarton
from target_selection.cartons.tools import create_table_as


class MWM_CB_GALEX_Mag_boss_Carton(BaseCarton):
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

    name = "mwm_cb_galex_mag_boss"
    mapper = "MWM"
    category = "science"
    instrument = "BOSS"
    cadence = None  # cadence is set in peewee.Case() below
    program = "mwm_cb"
    priority = 1820
    can_offset = True

    # cadence is set in peewee.Case() below
    # For cadence:
    # G < =16, bright_flexible_2x1
    # G > 16, dark_flexible_2x1
    def build_query(self, version_id, query_region=None):
        G3 = Gaia_DR3
        pm = fn.sqrt(G3.pmra * G3.pmra + G3.pmdec * G3.pmdec)
        pm_error = fn.sqrt(
            (G3.pmra_error * G3.pmra / pm) * (G3.pmra_error * G3.pmra / pm)
            + (G3.pmdec_error * G3.pmdec / pm) * (G3.pmdec_error * G3.pmdec / pm)
        )
        logpmdpm = fn.log(pm / pm_error)
        poe = G3.parallax / G3.parallax_error

        fuv_mag = GUVCat.fuv_mag
        Gm = G3.phot_g_mean_mag
        parallax = G3.parallax
        logpar = fn.log(fn.abs(parallax / 1000))
        ll = G3.l
        bb = G3.b

        fuv_cut = ((fuv_mag + 5 * logpar + 5) > (1.5 + 1.28 * (fuv_mag - Gm))) & (
            (fuv_mag + 5 * logpar + 5) > (2.5 * (fuv_mag - Gm) - 0.5)
        )

        astro_cut1 = ~(
            (
                (logpmdpm < 0.301)
                & (poe > (-1.75 * logpmdpm - 2.8))
                & (poe < (1.75 * logpmdpm + 2.8))
            )
            | (
                (
                    (logpmdpm - 0.301) * (logpmdpm - 0.301) / (0.33 * 0.33)
                    + (poe * poe) / (3.2 * 3.2)
                )
                <= 1
            )
        )

        astro_cut2 = ~(
            (fn.sqrt(fn.pow(ll - 280.3, 2) + 2 * fn.pow(bb + 33.0, 2)) < 8)
            | (fn.sqrt(0.5 * fn.pow(ll - 301, 2) + 2 * fn.pow(bb + 43.3, 2)) < 5)
        )

        cadence = peewee.Case(
            None, ((Gm <= 16, "bright_flexible_2x1"), (Gm > 16, "dark_flexible_2x1"))
        )

        query = (
            Gaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                G3.source_id,
                ll,
                bb,
                parallax,
                Gm,
                fuv_mag,
                cadence.alias("cadence"),
            )
            .join(
                Galex_GR7_Gaia_DR3,
                on=(Gaia_DR3.source_id == Galex_GR7_Gaia_DR3.gaia_edr3_source_id),
            )
            .join(GUVCat, on=(GUVCat.objid == Galex_GR7_Gaia_DR3.galex_objid))
            .join_from(Gaia_DR3, CatalogToGaia_DR3)
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                fuv_cut,
                fuv_mag > -999,
                GUVCat.fuv_magerr < 0.3,
                astro_cut1,
                astro_cut2,
            )
            .order_by(CatalogToGaia_DR3.catalogid, Galex_GR7_Gaia_DR3.galex_separation)
            .distinct(CatalogToGaia_DR3.catalogid)
        )

        if query_region:
            query = query.join_from(CatalogToGaia_DR3, Catalog).where(
                peewee.fn.q3c_radial_query(
                    Catalog.ra,
                    Catalog.dec,
                    query_region[0],
                    query_region[1],
                    query_region[2],
                )
            )

        return query


class MWM_CB_GALEX_Vol_boss_Carton(BaseCarton):
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
        1. Eliminate AGN - 80 632 559
        Condition: 1" coordinate match, select all objects that do not match an entry in
        MILLIQUAS_7_7 or that match and have QPCT!=100.

        2. Eliminate Galaxies - 77 074 394
        Condition: Eliminate targets that are within 60*R1 of a HECATE Galaxy

        3. Keep matches with GAIA-GALEX coordinate separation below 2.5" - 70 608 597
        Condition: GALEX_Separation<2.5

        4. Valid nuv_mags and nuv error - 30 204 026
        Condition: nuv_mag>-999 && nuv_magerr<0.2 && nuv_magerr>0

        5. HRD/CCD cut - 2 012 017
        Condition: ((nuv_mag-5*log10(r_med_geo)+5>2*(nuv_mag-phot_g_mean_mag)+0.7)
        or(nuv_mag-5*log10(abs(1000/parallax))+5>2*(nuv_mag-phot_g_mean_mag)+0.7)) and
        (nuv_mag-phot_g_mean_mag<=1.6*bp_rp+1.2)"

        6. Astrometry cut - 1 031 545
        Condition:
        (NOT(((log10(sqrt(pmra*pmra+pmdec*pmdec)/sqrt((pmra_error*pmra/sqrt(pmra*pmra+pmdec*pmdec))
        *(pmra_error*pmra/sqrt(pmra*pmra+pmdec*pmdec))+(pmdec_error*pmdec/sqrt(pmra*pmra+pmdec*
        pmdec))*(pmdec_error*pmdec/sqrt(pmra*pmra+pmdec*pmdec))))<0.301 AND
        (parallax/parallax_error)>-1.75*log10(sqrt(pmra*pmra+pmdec*pmdec)/
        sqrt((pmra_error*pmra/sqrt(pmra*pmra+pmdec*pmdec))*(pmra_error*pmra/
        sqrt(pmra*pmra+pmdec*pmdec))+(pmdec_error*pmdec/sqrt(pmra*pmra+pmdec*pmdec))*
        (pmdec_error*pmdec/sqrt(pmra*pmra+pmdec*pmdec))))-2.8 AND
        (parallax/parallax_error)<1.75*log10(sqrt(pmra*pmra+pmdec*pmdec)/
        sqrt((pmra_error*pmra/sqrt(pmra*pmra+pmdec*pmdec))*(pmra_error*pmra/
        sqrt(pmra*pmra+pmdec*pmdec))+(pmdec_error*pmdec/sqrt(pmra*pmra+pmdec*pmdec))*
        (pmdec_error*pmdec/sqrt(pmra*pmra+pmdec*pmdec))))+2.8) OR
        (log10(sqrt(pmra*pmra+pmdec*pmdec)/sqrt((pmra_error*pmra/sqrt(pmra*pmra+pmdec*pmdec))*
        (pmra_error*pmra/sqrt(pmra*pmra+pmdec*pmdec))+(pmdec_error*pmdec/
        sqrt(pmra*pmra+pmdec*pmdec))*(pmdec_error*pmdec/sqrt(pmra*pmra+pmdec*pmdec))))-0.301)*
        (log10(sqrt(pmra*pmra+pmdec*pmdec)/sqrt((pmra_error*pmra/sqrt(pmra*pmra+pmdec*pmdec))*
        (pmra_error*pmra/sqrt(pmra*pmra+pmdec*pmdec))+(pmdec_error*pmdec/
        sqrt(pmra*pmra+pmdec*pmdec))*(pmdec_error*pmdec/sqrt(pmra*pmra+pmdec*pmdec))))-0.301)/
        (0.33*0.33)+(parallax/parallax_error)*(parallax/parallax_error)/(3.2*3.2)<=1)))

        7. Distance cut - 564 537
        Condition: (r_med_geo<=1500|| abs(1000/parallax)<=1500)

        8a. Eliminate GAIA duplicates - 408 446
        Condition: Sort result from step 7 from lowest to largest matching distance and
        keep only closest match for each GAIA source_id

        8b. Eliminate GALEX duplicates - 560 021
        Condition: Sort result from step 7 from lowest to largest matching distance and
        keep only closest match for each GALEX ID

        9. Combine unique catalogues, keeping only single detections both ways - 405 156
        Condition: Match result catalogues from steps 8a and 8b via GAIA source_id -> Numbers in
        agreement with what is shown in uvex_selection_final_numbers.txt

        10. Remove objects in Magellanic Clouds - 402 832
        Condition: !((sqrt( pow(l-280.3,2) + 2*pow(b+33.0,2) ) < 8) || sqrt(0.5*square(l-301) +
                                                                            2*square(b+43.3) ) < 5)

    """

    name = "mwm_cb_galex_vol_boss"
    mapper = "MWM"
    category = "science"
    instrument = "BOSS"
    cadence = None  # cadence is set in peewee.Case() below
    program = "mwm_cb"
    priority = 1820
    can_offset = True

    # cadence is set in peewee.Case() below
    # For cadence:
    # G <= 16, bright_flexible_2x1
    # G > 16, dark_flexible_2x1
    def build_query(self, version_id, query_region=None):
        self.database.execute_sql('SET work_mem = "10GB";')

        GR7 = Galex_GR7_Gaia_DR3
        G3 = Gaia_DR3

        nuv_cut = (
            GR7.select(GR7.gaia_edr3_source_id, GR7.galex_objid, GR7.galex_separation, GR7.nuv_mag)
            .where(
                GR7.nuv_mag > -999,
                GR7.nuv_magerr < 0.2,
                GR7.nuv_magerr > 0,
                GR7.galex_separation <= 2.5,
            )
            .cte("nuv_cut")
        )

        Gm = G3.phot_g_mean_mag
        cadence = peewee.Case(
            None, ((Gm <= 16, "bright_flexible_2x1"), (Gm > 16, "dark_flexible_2x1"))
        )

        gaia_query = (
            G3.select(
                G3.source_id,
                G3.pmra,
                G3.pmdec,
                G3.ra,
                G3.dec,
                G3.l,
                G3.b,
                G3.parallax,
                G3.parallax_error,
                G3.pmra_error,
                G3.pmdec_error,
                G3.phot_g_mean_mag,
                G3.bp_rp,
                cadence.alias("cadence"),
                nuv_cut.c.nuv_mag,
                nuv_cut.c.galex_objid,
                nuv_cut.c.galex_separation,
                BailerJonesEDR3.r_med_geo,
            )
            .join(nuv_cut, on=(nuv_cut.c.gaia_edr3_source_id == G3.source_id))
            .join(BailerJonesEDR3, on=(BailerJonesEDR3.source_id == G3.source_id))
            .with_cte(nuv_cut)
        )

        log.debug("Collecting Gaia-GALEX information ...")

        cb_temp1, _ = create_table_as(
            gaia_query, "cb_temp1", temporary=True, database=self.database, overwrite=True
        )

        log.debug(f"Gaia-Galex NUV cut returned {cb_temp1.select().count():,} rows.")

        sqrt = fn.sqrt
        pow = fn.pow
        log10 = fn.log10
        pabs = fn.abs

        pmra = cb_temp1.c.pmra
        pmdec = cb_temp1.c.pmdec
        pmra_error = cb_temp1.c.pmra_error
        pmdec_error = cb_temp1.c.pmdec_error
        parallax = cb_temp1.c.parallax
        parallax_error = cb_temp1.c.parallax_error

        pm = sqrt(pmra * pmra + pmdec * pmdec)
        pm_error = sqrt(
            (pmra_error * pmra / pm) * (pmra_error * pmra / pm)
            + (pmdec_error * pmdec / pm) * (pmdec_error * pmdec / pm)
        )
        perror2 = (parallax / parallax_error) * (parallax / parallax_error)

        # Astro cut
        astro_cut = ~(
            (
                (log10(pm / pm_error) < 0.301)
                & ((parallax / parallax_error) > -1.75 * log10(pm / pm_error) - 2.8)
                & ((parallax / parallax_error) < 1.75 * log10(pm / pm_error) + 2.8)
            )
            | (
                (log10(pm / pm_error) - 0.301) * (log10(pm / pm_error) - 0.301) / (0.33**2)
                + perror2 / (3.2**2)
                <= 1
            )
        )

        nuv_mag = cb_temp1.c.nuv_mag
        r_med_geo = cb_temp1.c.r_med_geo
        phot_g_mean_mag = cb_temp1.c.phot_g_mean_mag
        bp_rp = cb_temp1.c.bp_rp

        d_pc = pabs(1000 / parallax)
        hrd_bj = (nuv_mag - 5 * log10(r_med_geo) + 5) > (2 * (nuv_mag - phot_g_mean_mag) + 0.7)
        hrd_gaia = (nuv_mag - 5 * log10(d_pc) + 5) > (2 * (nuv_mag - phot_g_mean_mag) + 0.7)

        # HRD/CCD cut
        hrd_cut = (hrd_bj | hrd_gaia) & (nuv_mag - phot_g_mean_mag <= (1.6 * bp_rp + 1.2))

        # Distance cut
        distance_cut = (d_pc <= 1500) | (cb_temp1.c.r_med_geo <= 1500)

        # MC cut
        ll = cb_temp1.c.l
        bb = cb_temp1.c.b
        mc_cut = ~(
            (sqrt(pow(ll - 280.3, 2) + 2 * pow(bb + 33.0, 2)) < 8)
            | (sqrt(0.5 * pow(ll - 301, 2) + 2 * pow(bb + 43.3, 2)) < 5)
        )

        cb_temp2_query = cb_temp1.select(
            cb_temp1.c.source_id,
            cb_temp1.c.ra,
            cb_temp1.c.dec,
            cb_temp1.c.galex_objid,
            cb_temp1.c.galex_separation,
            cb_temp1.c.cadence,
        ).where(astro_cut, hrd_cut, distance_cut, mc_cut)

        log.debug("Applying astrometric and colour cuts ...")

        cb_temp2, _ = create_table_as(
            cb_temp2_query, "cb_temp2", temporary=True, database=self.database, overwrite=True
        )

        log.debug(f"Astrometric and colour cuts returned {cb_temp2.select().count():,} rows.")

        # Remove AGNs, although this doesn't seem to actually reject any.
        agn_query = (
            cb_temp2.select(
                CatalogToGaia_DR3.catalogid,
                cb_temp2.c.source_id,
                cb_temp2.c.ra,
                cb_temp2.c.dec,
                cb_temp2.c.galex_objid,
                cb_temp2.c.galex_separation,
                cb_temp2.c.cadence,
            )
            .join(CatalogToGaia_DR3, on=(CatalogToGaia_DR3.target_id == cb_temp2.c.source_id))
            .join(
                CatalogToMILLIQUAS_7_7,
                on=(CatalogToGaia_DR3.catalogid == CatalogToMILLIQUAS_7_7.catalogid),
                join_type=peewee.JOIN.LEFT_OUTER,
            )
            .join(
                MILLIQUAS_7_7,
                on=(MILLIQUAS_7_7.pk == CatalogToMILLIQUAS_7_7.target_id),
                join_type=peewee.JOIN.LEFT_OUTER,
            )
            .where(
                CatalogToMILLIQUAS_7_7.best.is_null()
                | ((CatalogToMILLIQUAS_7_7.best >> True) & (MILLIQUAS_7_7.qpct != 100))
            )
            .where(CatalogToGaia_DR3.best >> True, CatalogToGaia_DR3.version_id == version_id)
        )

        cb_temp3, _ = create_table_as(
            agn_query,
            "cb_temp3",
            temporary=True,
            database=self.database,
            overwrite=True,
            indices=["q3c_ang2ipix(ra,dec)"],
        )

        log.debug(f"Remaining after rejecting AGNs: {cb_temp3.select().count():,} rows")

        # Select targets without nearby galaxies. First create a list of galaxies
        # within hecate_1_1.r1 arcmin of each of the targets.
        gal_query = (
            cb_temp3.select(cb_temp3.c.source_id)
            .join(HECATE_1_1, join_type=peewee.JOIN.CROSS)
            .where(
                fn.q3c_join(
                    HECATE_1_1.ra,
                    HECATE_1_1.dec,
                    cb_temp3.c.ra,
                    cb_temp3.c.dec,
                    HECATE_1_1.r1 / 60.0,
                ),
                HECATE_1_1.r1 != "NaN",
            )
            .cte("hecate_query", materialized=True)
        )

        # Now keep only targets that are not in the previous list.
        gal_free_query = (
            cb_temp3.select(
                cb_temp3.c.catalogid,
                cb_temp3.c.source_id,
                cb_temp3.c.ra,
                cb_temp3.c.dec,
                cb_temp3.c.galex_objid,
                cb_temp3.c.galex_separation,
                cb_temp3.c.cadence,
            )
            .join(
                gal_query,
                on=(gal_query.c.source_id == cb_temp3.c.source_id),
                join_type=peewee.JOIN.LEFT_OUTER,
            )
            .where(gal_query.c.source_id.is_null())
            .with_cte(gal_query)
        )

        cb_temp4, _ = create_table_as(
            gal_free_query, "cb_temp4", temporary=True, database=self.database, overwrite=True
        )

        log.debug(
            "Remaining after filtering for nearby galaxies: "
            f"{cb_temp4.select().count():,} rows."
        )

        # Get unique Gaia and GALEX targets with minimum separaton, independently, and keep
        # source_ids that are in both groups.

        # First we create a subquery partitioned by source_id and rank each row by separation.
        gaia_ranked = cb_temp4.select(
            cb_temp4.c.source_id,
            fn.rank()
            .over(
                partition_by=[cb_temp4.c.source_id], order_by=[cb_temp4.c.galex_separation.asc()]
            )
            .alias("sep_rank"),
        ).alias("gaia_ranked")

        # Now we select only the entry with sep_rank=1 (minimum galex_separation).
        gaia_unique = (
            cb_temp4.select(cb_temp4.c.source_id)
            .join(gaia_ranked, on=(cb_temp4.c.source_id == gaia_ranked.c.source_id))
            .where(gaia_ranked.c.sep_rank == 1)
            .distinct(cb_temp4.c.source_id)
            .alias("gaia_unique")
        )

        # Do the same for galex_objid but keep the Gaia source since that's the one we'll join on.
        galex_ranked = cb_temp4.select(
            cb_temp4.c.source_id,
            fn.rank()
            .over(
                partition_by=[cb_temp4.c.galex_objid], order_by=[cb_temp4.c.galex_separation.asc()]
            )
            .alias("sep_rank"),
        ).alias("galex_ranked")

        galex_unique = (
            cb_temp4.select(cb_temp4.c.source_id)
            .join(galex_ranked, on=(cb_temp4.c.source_id == galex_ranked.c.source_id))
            .where(galex_ranked.c.sep_rank == 1)
            .distinct(cb_temp4.c.galex_objid)
            .alias("galex_unique")
        )

        # Now keep source_ids that exist in both subsamples. Note we need to distinct on
        # source_id again because galex_unique has duplicate source_ids (different galex_objid
        # can have the same source_id).
        unique = (
            cb_temp4.select(
                cb_temp4.c.catalogid,
                cb_temp4.c.source_id,
                cb_temp4.c.ra,
                cb_temp4.c.dec,
                cb_temp4.c.galex_objid,
                cb_temp4.c.galex_separation,
                cb_temp4.c.cadence,
            )
            .join(gaia_unique, on=(gaia_unique.c.source_id == cb_temp4.c.source_id))
            .join(galex_unique, on=(galex_unique.c.source_id == cb_temp4.c.source_id))
            .distinct(cb_temp4.c.source_id)
        )

        log.debug("Removing duplicates ...")

        return unique


class MWM_CB_XMMOM_boss_Carton(BaseCarton):
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
          - UVM2_AB_mag - phot_g_mean_mag < 1.6 * (bp-rp) + 1.2
        * Astrometric cut: to remove AGN contaminants:
            !(logpmdpm < 0.301 && poe > -1.75*logpmdpm - 2.8 && poe <  1.75*logpmdpm + 2.8 ||
            (logpmdpm-0.301)*(logpmdpm-0.301) / (0.33*0.33) + (poe * poe) / (3.2* 3.2)<=1)
        * Distance cut: r_med_geo<=1500 ||  1000/parallax <= 1500
        * Eliminate duplicates both ways (Gaia and XMMOM), keep always nearest
          Gaia - XMMOM match only.

    """

    name = "mwm_cb_xmmom_boss"
    mapper = "MWM"
    category = "science"
    instrument = "BOSS"
    cadence = None  # cadence is set in peewee.Case() below
    program = "mwm_cb"
    priority = 1820
    can_offset = True

    # cadence is set in peewee.Case() below
    # For cadence:
    # G <= 16, bright_flexible_2x1
    # G > 16, dark_flexible_2x1
    def build_query(self, version_id, query_region=None):
        G3 = Gaia_DR3
        pm = fn.sqrt(G3.pmra * G3.pmra + G3.pmdec * G3.pmdec)
        pm_error = fn.sqrt(
            (G3.pmra_error * G3.pmra / pm) * (G3.pmra_error * G3.pmra / pm)
            + (G3.pmdec_error * G3.pmdec / pm) * (G3.pmdec_error * G3.pmdec / pm)
        )
        logpmdpm = fn.log(pm / pm_error)
        poe = G3.parallax / G3.parallax_error

        Gm = G3.phot_g_mean_mag
        bp_rp = G3.bp_rp
        parallax = G3.parallax
        logpar = fn.log(fn.abs(1000 / parallax))

        r_med_geo = BailerJonesEDR3.r_med_geo

        XMM = XMM_OM_SUSS_5_0
        uvm2 = XMM.uvm2_ab_mag

        qcut = ~(
            (fn.substring(XMM.uvm2_quality_flag_st, 1, 1) == "T")
            | ((fn.substring(XMM.uvm2_quality_flag_st, 2, 1) == "T") & (XMM.uvm2_signif < 10))
            | ((fn.substring(XMM.uvm2_quality_flag_st, 3, 1) == "T") & (XMM.uvm2_signif < 10))
            | (fn.substring(XMM.uvm2_quality_flag_st, 7, 1) == "T")
            | (fn.substring(XMM.uvm2_quality_flag_st, 8, 1) == "T")
            | (fn.substring(XMM.uvm2_quality_flag_st, 9, 1) == "T")
        )

        color_cuts = [
            (uvm2 >= -100),
            (
                ((uvm2 - 5 * fn.log10(r_med_geo) + 5) > (2 * (uvm2 - Gm) + 0.7))
                | ((uvm2 - 5 * logpar + 5) > (2 * (uvm2 - Gm) + 0.7))
            ),
            (uvm2 - Gm < 1.6 * bp_rp + 1.2),
        ]

        astro_cut = ~(
            (logpmdpm < 0.301) & (poe > (-1.75 * logpmdpm - 2.8)) & (poe < (1.75 * logpmdpm + 2.8))
            | (
                (
                    (logpmdpm - 0.301) * (logpmdpm - 0.301) / (0.33 * 0.33)
                    + (poe * poe) / (3.2 * 3.2)
                )
                <= 1
            )
        )

        cadence = peewee.Case(
            None, ((Gm <= 16, "bright_flexible_2x1"), (Gm > 16, "dark_flexible_2x1"))
        )

        # MC cut
        ll = G3.l
        bb = G3.b
        mc_cut = ~(
            (fn.sqrt(fn.pow(ll - 280.3, 2) + 2 * fn.pow(bb + 33.0, 2)) < 8)
            | (fn.sqrt(0.5 * fn.pow(ll - 301, 2) + 2 * fn.pow(bb + 43.3, 2)) < 5)
        )

        q1 = (
            Gaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                G3.source_id,
                G3.ra,
                G3.dec,
                CatalogToXMM_OM_SUSS_5_0.distance.alias("separation"),
                XMM_OM_SUSS_5_0.srcnum.alias("srcnum"),
                cadence.alias("cadence"),
            )
            .join_from(
                Gaia_DR3, BailerJonesEDR3, on=(Gaia_DR3.source_id == BailerJonesEDR3.source_id)
            )
            .join_from(Gaia_DR3, CatalogToGaia_DR3)
            .join_from(
                CatalogToGaia_DR3,
                CatalogToMILLIQUAS_7_7,
                join_type=peewee.JOIN.LEFT_OUTER,
                on=(CatalogToMILLIQUAS_7_7.catalogid == CatalogToGaia_DR3.catalogid),
            )
            .join(MILLIQUAS_7_7, join_type=peewee.JOIN.LEFT_OUTER)
            .join_from(
                CatalogToGaia_DR3,
                CatalogToXMM_OM_SUSS_5_0,
                on=(CatalogToGaia_DR3.catalogid == CatalogToXMM_OM_SUSS_5_0.catalogid),
            )
            .join(XMM_OM_SUSS_5_0)
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                CatalogToXMM_OM_SUSS_5_0.distance <= 1.5 / 3600.0,
                (
                    CatalogToMILLIQUAS_7_7.catalogid.is_null()
                    | ((CatalogToMILLIQUAS_7_7.best >> True) & (MILLIQUAS_7_7.qpct != 100))
                ),
                (r_med_geo <= 1500) | (fn.abs(1000 / parallax) <= 1500),
                qcut,
                *color_cuts,
                astro_cut,
                mc_cut,
            )
        )

        xmm_temp, _ = create_table_as(
            q1,
            "mwm_cb_xmmom_temp",
            temporary=True,
            database=self.database,
            indices=["q3c_ang2ipix(ra,dec)"],
        )

        # Select targets without nearby galaxies. First create a list of galaxies
        # within hecate_1_1.r1 arcmin of each of the targets.
        gal_query = (
            xmm_temp.select(xmm_temp.c.source_id)
            .join(HECATE_1_1, join_type=peewee.JOIN.CROSS)
            .where(
                fn.q3c_join(
                    HECATE_1_1.ra,
                    HECATE_1_1.dec,
                    xmm_temp.c.ra,
                    xmm_temp.c.dec,
                    HECATE_1_1.r1 / 60.0,
                ),
                HECATE_1_1.r1 != "NaN",
            )
            .cte("hecate_query", materialized=True)
        )

        # Now keep only targets that are not in the previous list.
        gal_free_query = (
            xmm_temp.select(
                xmm_temp.c.catalogid,
                xmm_temp.c.source_id,
                xmm_temp.c.srcnum,
                xmm_temp.c.separation,
                xmm_temp.c.ra,
                xmm_temp.c.dec,
                xmm_temp.c.cadence,
            )
            .join(
                gal_query,
                on=(gal_query.c.source_id == xmm_temp.c.source_id),
                join_type=peewee.JOIN.LEFT_OUTER,
            )
            .where(gal_query.c.source_id.is_null())
            .with_cte(gal_query)
        )

        xmm_temp2, _ = create_table_as(
            gal_free_query, "mwm_cb_xmmom_temp2", temporary=True, database=self.database
        )

        # Get unique Gaia and XMMOM targets with minimum separaton, independently, and keep
        # source_ids that are in both groups.

        # First we create a subquery partitioned by source_id and rank each row by separation.
        gaia_ranked = xmm_temp2.select(
            xmm_temp2.c.source_id,
            fn.rank()
            .over(partition_by=[xmm_temp2.c.source_id], order_by=[xmm_temp2.c.separation.asc()])
            .alias("sep_rank"),
        ).alias("gaia_ranked")

        # Now we select only the entry with sep_rank=1 (minimum separation).
        gaia_unique = (
            xmm_temp2.select(xmm_temp2.c.source_id)
            .join(gaia_ranked, on=(xmm_temp2.c.source_id == gaia_ranked.c.source_id))
            .where(gaia_ranked.c.sep_rank == 1)
            .distinct(xmm_temp2.c.source_id)
            .alias("gaia_unique")
        )

        # Do the same for XMMOM id but keep the Gaia source since that's the one we'll join on.
        xmm_ranked = xmm_temp2.select(
            xmm_temp2.c.source_id,
            fn.rank()
            .over(partition_by=[xmm_temp2.c.srcnum], order_by=[xmm_temp2.c.separation.asc()])
            .alias("sep_rank"),
        ).alias("xmm_ranked")

        xmm_unique = (
            xmm_temp2.select(xmm_temp2.c.source_id)
            .join(xmm_ranked, on=(xmm_temp2.c.source_id == xmm_ranked.c.source_id))
            .where(xmm_ranked.c.sep_rank == 1)
            .distinct(xmm_temp2.c.srcnum)
            .alias("xmm_unique")
        )

        # Now keep source_ids that exist in both subsamples. Note we need to distinct on
        # source_id again because xmm_unique has duplicate source_ids (different srcnum
        # can have the same source_id).
        unique = (
            xmm_temp2.select(
                xmm_temp2.c.catalogid,
                xmm_temp2.c.source_id,
                xmm_temp2.c.ra,
                xmm_temp2.c.dec,
                xmm_temp2.c.srcnum,
                xmm_temp2.c.separation,
                xmm_temp2.c.cadence,
            )
            .join(gaia_unique, on=(gaia_unique.c.source_id == xmm_temp2.c.source_id))
            .join(xmm_unique, on=(xmm_unique.c.source_id == xmm_temp2.c.source_id))
            .distinct(xmm_temp2.c.source_id)
        )

        return unique


class MWM_CB_SWIFTUVOT_boss_Carton(BaseCarton):
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

    name = "mwm_cb_swiftuvot_boss"
    mapper = "MWM"
    category = "science"
    program = "mwm_cb"
    instrument = "BOSS"
    cadence = None  # cadence is set in peewee.Case() below
    priority = 1820

    # cadence is set in peewee.Case() below
    # For cadence:
    # G <= 16, bright_flexible_2x1
    # G > 16, dark_flexible_2x1
    def build_query(self, version_id, query_region=None):
        G3 = Gaia_DR3
        pm = fn.sqrt(G3.pmra * G3.pmra + G3.pmdec * G3.pmdec)
        pm_error = fn.sqrt(
            (G3.pmra_error * G3.pmra / pm) * (G3.pmra_error * G3.pmra / pm)
            + (G3.pmdec_error * G3.pmdec / pm) * (G3.pmdec_error * G3.pmdec / pm)
        )
        logpmdpm = fn.log(pm / pm_error)
        poe = G3.parallax / G3.parallax_error

        Gm = G3.phot_g_mean_mag
        bp_rp = G3.bp_rp

        parallax = G3.parallax
        logpar = fn.log(fn.abs(1000 / parallax))

        r_med_geo = BailerJonesEDR3.r_med_geo

        uvw2_ab_mag = UVOT_SSC_1.uvw2_ab
        uvm2 = UVOT_SSC_1.uvm2_ab
        uvw1_ab_mag = UVOT_SSC_1.uvw1_ab

        uvm2_quality_flag = UVOT_SSC_1.fuvm2
        uvm2_signif = UVOT_SSC_1.suvm2

        qcut = (
            (uvw2_ab_mag > -100) | (uvw1_ab_mag > -100) | (uvm2 > -100),
            ~(
                (uvm2_quality_flag.bin_and(2**0) > 0)
                | (uvm2_quality_flag.bin_and(2**1) > 0) & (uvm2_signif < 10)
                | (uvm2_quality_flag.bin_and(2**2) > 0) & (uvm2_signif < 10)
                | (uvm2_quality_flag.bin_and(2**5) > 0) & (uvm2_signif < 10)
                | (uvm2_quality_flag.bin_and(2**6) > 0)
            ),
        )

        color_cuts = [
            uvm2 > -100,
            ((uvm2 - 5 * fn.log10(r_med_geo) + 5) > (2 * (uvm2 - Gm) + 0.7))
            | ((uvm2 - 5 * logpar + 5) > (2 * (uvm2 - Gm) + 0.7)),
            (uvm2 - Gm < 1.6 * bp_rp + 1.2),
        ]

        astro_cut = ~(
            (logpmdpm < 0.301) & (poe > (-1.75 * logpmdpm - 2.8)) & (poe < (1.75 * logpmdpm + 2.8))
            | (
                (
                    (logpmdpm - 0.301) * (logpmdpm - 0.301) / (0.33 * 0.33)
                    + (poe * poe) / (3.2 * 3.2)
                )
                <= 1
            )
        )

        cadence = peewee.Case(
            None, ((Gm <= 16, "bright_flexible_2x1"), (Gm > 16, "dark_flexible_2x1"))
        )

        # MC cut
        ll = G3.l
        bb = G3.b
        mc_cut = ~(
            (fn.sqrt(fn.pow(ll - 280.3, 2) + 2 * fn.pow(bb + 33.0, 2)) < 8)
            | (fn.sqrt(0.5 * fn.pow(ll - 301, 2) + 2 * fn.pow(bb + 43.3, 2)) < 5)
        )

        q1 = (
            Gaia_DR3.select(
                CatalogToGaia_DR3.catalogid,
                G3.source_id,
                G3.ra,
                G3.dec,
                UVOT_SSC_1.id.alias("uvot_id"),
                CatalogToUVOT_SSC_1.distance.alias("separation"),
                cadence.alias("cadence"),
            )
            .join_from(
                Gaia_DR3, BailerJonesEDR3, on=(Gaia_DR3.source_id == BailerJonesEDR3.source_id)
            )
            .join_from(Gaia_DR3, CatalogToGaia_DR3)
            .join_from(
                CatalogToGaia_DR3,
                CatalogToMILLIQUAS_7_7,
                join_type=peewee.JOIN.LEFT_OUTER,
                on=(CatalogToMILLIQUAS_7_7.catalogid == CatalogToGaia_DR3.catalogid),
            )
            .join(MILLIQUAS_7_7, join_type=peewee.JOIN.LEFT_OUTER)
            .join_from(
                CatalogToGaia_DR3,
                CatalogToUVOT_SSC_1,
                on=(CatalogToGaia_DR3.catalogid == CatalogToUVOT_SSC_1.catalogid),
            )
            .join(UVOT_SSC_1)
            .where(
                CatalogToGaia_DR3.version_id == version_id,
                CatalogToGaia_DR3.best >> True,
                (
                    CatalogToMILLIQUAS_7_7.catalogid.is_null()
                    | ((CatalogToMILLIQUAS_7_7.best >> True) & (MILLIQUAS_7_7.qpct != 100))
                ),
                CatalogToUVOT_SSC_1.distance <= 1.5 / 3600.0,
                (r_med_geo <= 1500) | (fn.abs(1000 / parallax) <= 1500),
                *qcut,
                *color_cuts,
                astro_cut,
                mc_cut,
            )
            .order_by(CatalogToGaia_DR3.catalogid, CatalogToUVOT_SSC_1.distance)
            .distinct(CatalogToGaia_DR3.catalogid)
        )

        uvot_temp, _ = create_table_as(
            q1,
            "mwm_cb_uvot_temp",
            temporary=True,
            database=self.database,
            indices=["q3c_ang2ipix(ra,dec)"],
        )

        # Select targets without nearby galaxies. First create a list of galaxies
        # within hecate_1_1.r1 arcmin of each of the targets.
        gal_query = (
            uvot_temp.select(uvot_temp.c.source_id)
            .join(HECATE_1_1, join_type=peewee.JOIN.CROSS)
            .where(
                fn.q3c_join(
                    HECATE_1_1.ra,
                    HECATE_1_1.dec,
                    uvot_temp.c.ra,
                    uvot_temp.c.dec,
                    HECATE_1_1.r1 / 60.0,
                ),
                HECATE_1_1.r1 != "NaN",
            )
            .cte("hecate_query", materialized=True)
        )

        # Now keep only targets that are not in the previous list.
        gal_free_query = (
            uvot_temp.select(
                uvot_temp.c.catalogid,
                uvot_temp.c.source_id,
                uvot_temp.c.uvot_id,
                uvot_temp.c.separation,
                uvot_temp.c.ra,
                uvot_temp.c.dec,
                uvot_temp.c.cadence,
            )
            .join(
                gal_query,
                on=(gal_query.c.source_id == uvot_temp.c.source_id),
                join_type=peewee.JOIN.LEFT_OUTER,
            )
            .where(gal_query.c.source_id.is_null())
            .with_cte(gal_query)
        )

        uvot_temp2, _ = create_table_as(
            gal_free_query, "mwm_cb_uvotom_temp2", temporary=True, database=self.database
        )

        # Get unique Gaia and UVOT targets with minimum separaton, independently, and keep
        # source_ids that are in both groups.

        # First we create a subquery partitioned by source_id and rank each row by separation.
        gaia_ranked = uvot_temp2.select(
            uvot_temp2.c.source_id,
            fn.rank()
            .over(partition_by=[uvot_temp2.c.source_id], order_by=[uvot_temp2.c.separation.asc()])
            .alias("sep_rank"),
        ).alias("gaia_ranked")

        # Now we select only the entry with sep_rank=1 (minimum separation).
        gaia_unique = (
            uvot_temp2.select(uvot_temp2.c.source_id)
            .join(gaia_ranked, on=(uvot_temp2.c.source_id == gaia_ranked.c.source_id))
            .where(gaia_ranked.c.sep_rank == 1)
            .distinct(uvot_temp2.c.source_id)
            .alias("gaia_unique")
        )

        # Do the same for UVOT id but keep the Gaia source since that's the one we'll join on.
        uvot_ranked = uvot_temp2.select(
            uvot_temp2.c.source_id,
            fn.rank()
            .over(partition_by=[uvot_temp2.c.uvot_id], order_by=[uvot_temp2.c.separation.asc()])
            .alias("sep_rank"),
        ).alias("uvot_ranked")

        uvot_unique = (
            uvot_temp2.select(uvot_temp2.c.source_id)
            .join(uvot_ranked, on=(uvot_temp2.c.source_id == uvot_ranked.c.source_id))
            .where(uvot_ranked.c.sep_rank == 1)
            .distinct(uvot_temp2.c.uvot_id)
            .alias("uvot_unique")
        )

        # Now keep source_ids that exist in both subsamples. Note we need to distinct on
        # source_id again because uvot_unique has duplicate source_ids (different uvot_objid
        # can have the same source_id).
        unique = (
            uvot_temp2.select(
                uvot_temp2.c.catalogid,
                uvot_temp2.c.source_id,
                uvot_temp2.c.ra,
                uvot_temp2.c.dec,
                uvot_temp2.c.uvot_id,
                uvot_temp2.c.separation,
                uvot_temp2.c.cadence,
            )
            .join(gaia_unique, on=(gaia_unique.c.source_id == uvot_temp2.c.source_id))
            .join(uvot_unique, on=(uvot_unique.c.source_id == uvot_temp2.c.source_id))
            .distinct(uvot_temp2.c.source_id)
        )

        return unique
