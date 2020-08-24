#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Tom Dwelly
# @Date: 2020-03-31
# @Filename: bhm_rm.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
# derived from guide.py

# flake8: noqa
# isort: skip_file

import peewee
import sdssdb
from peewee import JOIN

#from sdssdb.peewee.sdss5db.catalogdb import (Catalog, BHM_RM_v0, CatalogToBHM_RM_v0)
from sdssdb.peewee.sdss5db.catalogdb import (Catalog,
                                             BHM_RM_v0_2,
                                             CatalogToBHM_RM_v0,
                                             CatalogToSDSS_DR16_SpecObj,
                                             SDSS_DR16_SpecObj)

from target_selection.cartons.base import BaseCarton


# Carton descriptions here:
# https://wiki.sdss.org/display/OPS/Defining+target+selection+and+cadence+algorithms#Definingtargetselectionandcadencealgorithms-ReverberationMapping(RM)

## This module provides the following BHM cartons:
#  bhm_rm_core
#  bhm_rm_known-spec
#  bhm_rm_var
#  bhm_rm_ancillary

pmsig_min = -3.0
plxsig_min = -3.0


class BhmRmBaseCarton(BaseCarton):
    '''
    This class provides common setting and the masking routines used by all RM cartons
    '''

    name = 'bhm_rm_base'
    base_name = 'bhm_rm_base'
    category = 'science'
    mapper = 'BHM'
    program = 'bhm_rm'
    cadence = 'bhm_rm_174x8'
    tile = False
    priority = None
    alias_c = None
    alias_t = None

    def get_fieldlist(self):
        '''Read the RM field centres from the yaml'''
        fieldlist = []
        base_parameters = self.config['parameters'].get(self.base_name, None)
        if base_parameters:
            fieldlist = base_parameters['fieldlist'];
        return fieldlist


    def append_spatial_query(self, query, fieldlist):
        '''extend the peewee query using a list of field centres'''
        if fieldlist is None :
            return query
        elif len(fieldlist) == 0 :
            return query

        q = False
        for f in fieldlist:
            q = ( q | peewee.fn.q3c_radial_query(self.alias_c.ra,
                                                 self.alias_c.dec,
                                                 f['racen'],
                                                 f['deccen'],
                                                 f['radius']))
        return query.where(q)


    def build_query(self, version_id, query_region=None):
        c = Catalog.alias()
        c2t = CatalogToBHM_RM_v0.alias()
        #t = BHM_RM_v0.alias()
        t = BHM_RM_v0_2.alias()
        self.alias_c = c
        self.alias_t = t
        c2s = CatalogToSDSS_DR16_SpecObj.alias()
        s = SDSS_DR16_SpecObj.alias()

        # set the Carton priority+values here - read from yaml
        #priority = peewee.Value(int(self.parameters.get('priority', 10000))).alias('priority')

        # fold in tiers of magnitude-based priority
        priority_mag_step = 0.5
        priority_mag_bright = 17.0
        priority_mag_faint = 22.0
        priority_mag_bright_known_spec = 20.5
        priority_floor = self.parameters.get('priority', 10000)
        priority = peewee.Case(
            None,
            (
                ((t.mi <= priority_mag_bright),
                 priority_floor+0),
                (((self.name == 'bhm_rm_known_spec') &
                  ~(t.field_name.contains('SDSS-RM')) &
                  (t.mi <= priority_mag_bright_known_spec)),
                 priority_floor+0),
                ((t.mi <= priority_mag_faint),
                 (priority_floor +
                  5*(1+peewee.fn.floor((t.mi-priority_mag_bright)/priority_mag_step).cast('int')))),
                ((t.mi > priority_mag_faint),
                 priority_floor+95),
            ),
            None
        ).cast('int')


        value = peewee.Value(self.parameters.get('value', 1.0)).cast('float').alias('value')
        pmra = peewee.Value(0.0).cast('float').alias('pmra')
        pmdec = peewee.Value(0.0).cast('float').alias('pmdec')
        parallax = peewee.Value(0.0).cast('float').alias('parallax')

        ### this leads to some scattered nulls in photometry - prefer to use best estimate photom instead
        # magnitude_g = peewee.Case(None,
        #                           (
        #                             ((t.optical_survey.contains('SDSS')) & (t.psfmag_sdss[1] > 0.0), t.psfmag_sdss[1]),
        #                             ((t.optical_survey.contains('PS1')) & (t.psfmag_ps1[0] > 0.0), t.psfmag_ps1[0]),
        #                             ((t.optical_survey.contains('DES')) & (t.psfmag_des[0] > 0.0), t.psfmag_des[0]),
        #                             ((t.optical_survey.contains('NSC')) & (t.mag_nsc[0] > 0.0), t.mag_nsc[0]),
        #                           ),
        #                           None) ## should never get here
        # magnitude_r = peewee.Case(None,
        #                           (
        #                             ((t.optical_survey.contains('SDSS')) & (t.psfmag_sdss[2] > 0.0), t.psfmag_sdss[2]),
        #                             ((t.optical_survey.contains('PS1')) & (t.psfmag_ps1[1] > 0.0), t.psfmag_ps1[1]),
        #                             ((t.optical_survey.contains('DES')) & (t.psfmag_des[1] > 0.0), t.psfmag_des[1]),
        #                             ((t.optical_survey.contains('NSC')) & (t.mag_nsc[1] > 0.0), t.mag_nsc[1]),
        #                           ),
        #                           None) ## should never get here
        # magnitude_i = peewee.Case(None,
        #                           (
        #                             ((t.optical_survey.contains('SDSS')) & (t.psfmag_sdss[3] > 0.0), t.psfmag_sdss[3]),
        #                             ((t.optical_survey.contains('PS1')) & (t.psfmag_ps1[2] > 0.0), t.psfmag_ps1[2]),
        #                             ((t.optical_survey.contains('DES')) & (t.psfmag_des[2] > 0.0), t.psfmag_des[2]),
        #                             ((t.optical_survey.contains('NSC')) & (t.mag_nsc[2] > 0.0), t.mag_nsc[2]),
        #                           ),
        #                           t.mi) ## should never get here
        # magnitude_z = peewee.Case(None,
        #                           (
        #                             ((t.optical_survey.contains('SDSS')) & (t.psfmag_sdss[4] > 0.0), t.psfmag_sdss[4]),
        #                             ((t.optical_survey.contains('PS1')) & (t.psfmag_ps1[3] > 0.0), t.psfmag_ps1[3]),
        #                             ((t.optical_survey.contains('DES')) & (t.psfmag_des[3] > 0.0), t.psfmag_des[3]),
        #                             ((t.optical_survey.contains('NSC')) & (t.mag_nsc[3] > 0.0), t.mag_nsc[3]),
        #                           ),
        #                           None) ## should never get here

        # Photometric precedence: DES>PS1>SDSS(>Gaia)>NSC.

        magnitude_g = peewee.Case(None,
                                  (
                                    ((t.des == 1) & (t.psfmag_des[0] > 0.0), t.psfmag_des[0]),
                                    ((t.ps1 == 1) & (t.psfmag_ps1[0] > 0.0), t.psfmag_ps1[0]),
                                    ((t.sdss == 1) & (t.psfmag_sdss[1] > 0.0), t.psfmag_sdss[1]),
                                    ((t.nsc == 1) & (t.mag_nsc[0] > 0.0), t.mag_nsc[0]),
                                  ),
                                  None) ## should never get here
        magnitude_r = peewee.Case(None,
                                  (
                                    ((t.des == 1) & (t.psfmag_des[1] > 0.0), t.psfmag_des[1]),
                                    ((t.ps1 == 1) & (t.psfmag_ps1[1] > 0.0), t.psfmag_ps1[1]),
                                    ((t.sdss == 1) & (t.psfmag_sdss[2] > 0.0), t.psfmag_sdss[2]),
                                    ((t.nsc == 1) & (t.mag_nsc[1] > 0.0), t.mag_nsc[1]),
                                  ),
                                  None) ## should never get here
        magnitude_i = peewee.Case(None,
                                  (
                                    ((t.des == 1) & (t.psfmag_des[2] > 0.0), t.psfmag_des[2]),
                                    ((t.ps1 == 1) & (t.psfmag_ps1[2] > 0.0), t.psfmag_ps1[2]),
                                    ((t.sdss == 1) & (t.psfmag_sdss[3] > 0.0), t.psfmag_sdss[3]),
                                    ((t.nsc == 1) & (t.mag_nsc[2] > 0.0), t.mag_nsc[2]),
                                  ),
                                  t.mi) ## should never get here
        magnitude_z = peewee.Case(None,
                                  (
                                    ((t.des == 1) & (t.psfmag_des[3] > 0.0), t.psfmag_des[3]),
                                    ((t.ps1 == 1) & (t.psfmag_ps1[3] > 0.0), t.psfmag_ps1[3]),
                                    ((t.sdss == 1) & (t.psfmag_sdss[4] > 0.0), t.psfmag_sdss[4]),
                                    ((t.nsc == 1) & (t.mag_nsc[3] > 0.0), t.mag_nsc[3]),
                                  ),
                                  None) ## should never get here


        ## now gaia mags
        #magnitude_bp = peewee.Case(None,
        #                           (
        #                            ((t.gaia == 1), t.mag_gaia[1]),
        #                           ),
        #                           None)
        #magnitude_rp = peewee.Case(None,
        #                           (
        #                            ((t.gaia == 1), t.mag_gaia[2]),
        #                           ),
        #                           None)



        query = (
            c
            .select(c.catalogid,
                    priority.alias('priority'),
                    value,
                    pmra,
                    pmdec,
                    parallax,
                    magnitude_g.alias('g'),
                    magnitude_r.alias('r'),
                    magnitude_i.alias('i'),
                    magnitude_z.alias('z'),
            )
            .join(c2t)
            # The following is needed because we are using c2t for Catalog_to_BHM_RM_v0
            # rather than a native c2t for Catalog_to_BHM_RM_v0_2
            .join(t, on=(c2t.target_id == t.pk))
            .where(c.version_id == version_id,
                   c2t.version_id == version_id,
                   (
                       (c2s.version_id == version_id) |
                       (c2s.version_id.is_null())
                   ),
                   c2t.best == True)
            .where
            (
                (t.mi >= self.parameters['mag_i_min']),
                (t.mi <  self.parameters['mag_i_max']),
            )
            .switch(c)
            .join(c2s, JOIN.LEFT_OUTER)
            .join(s, JOIN.LEFT_OUTER)
            .where(
                (s.specobjid.is_null()) |    # no match in sdss_dr16_specobj
                ~(
                    (s.class_.contains('STAR')) & # Reject objects where the best spectrum for
                    (s.scienceprimary > 0)      # this object is classified as STAR
                )
            )
            .distinct([t.pk])   # avoid duplicates - trust the RM parent sample
        )
        query = self.append_spatial_query(query, self.get_fieldlist())

        return query




class BhmRmCoreCarton(BhmRmBaseCarton):
    '''
    bhm_rm_core: select all photometric QSO targets with the likelihood method (Skewt), flux-limited to 21.5 in i-band PSF mag

    SELECT * FROM bhm_rm
            WHERE skewt_qso = 1
            AND WHERE mi BETWEEN 15.0 AND 21.5
            AND WHERE pmsig < 5.0
            AND WHERE plxsig < 5.0


    also require t.skewt_qso_prior == 1 in CVZ-S
    '''

    name = 'bhm_rm_core'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (t.skewt_qso == 1) &
            (
                ~(t.field_name.contains('S-CVZ')) |
                ((t.field_name.contains('S-CVZ')) & (t.skewt_qso_prior == 1))
            ) &
            (t.pmsig <  self.parameters['pmsig_max']) &
            (t.plxsig <  self.parameters['plxsig_max']) &
            (
                #catch bad photometry - require at least gri detections in at least one system
                ((t.sdss == 1) &
                 ( t.psfmag_sdss[1] > 0.0) &
                 ( t.psfmag_sdss[2] > 0.0) &
                 ( t.psfmag_sdss[3] > 0.0)) |
                ((t.ps1 == 1) &
                 ( t.psfmag_ps1[0] > 0.0) &
                 ( t.psfmag_ps1[1] > 0.0) &
                 ( t.psfmag_ps1[2] > 0.0)) |
                ((t.des == 1) &
                 ( t.psfmag_des[0] > 0.0) &
                 ( t.psfmag_des[1] > 0.0) &
                 ( t.psfmag_des[2] > 0.0)) |
                ((t.nsc == 1) &
                 ( t.mag_nsc[0] > 0.0) &
                 ( t.mag_nsc[1] > 0.0) &
                 ( t.mag_nsc[2] > 0.0))
            ) &
#            (t.pmsig > pmsig_min) &  # this catches cases with NULL=-9
#            (t.plxsig > plxsig_min) &
            ~(t.field_name.contains('SDSS-RM')) # ignore this carton in the SDSS-RM field
        )

        return query

class BhmRmKnownSpecCarton(BhmRmBaseCarton):
    '''
    bhm_rm_known-spec:  select all spectroscopically confirmed QSOs where redshift is extragalactic

    SELECT * FROM bhm_rm
        WHERE specz > 0.005
        AND WHERE mi BETWEEN 15.0 AND 21.5   # <- exact limits TBD

    For SDSS-RM field select only QSOs that were part of the original SDSS-III/IV RM program
    Do that via a bit in the spec_bitmask field. Also restrict to mi < 21.0
    '''

    name = 'bhm_rm_known-spec'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        spec_bitmask_sdss_rm_qso = 2**3
        query = query.where(
            (t.spec_q == 1 ) &
            (
                ~(t.field_name.contains('SDSS-RM')) |
                (
                    # include extra constraints on SDSS-RM targets
                    (t.spec_bitmask.bin_and(spec_bitmask_sdss_rm_qso) != 0) &
                    (t.mi < self.parameters['mag_i_max_sdss_rm'] )
                )
            ) &
            (
                ~(t.field_name.contains('COSMOS')) |
                (
                    # include extra constraints on COSMOS targets
                    (t.mi < self.parameters['mag_i_max_cosmos'] )
                )
            ) &
            (
                ~(t.field_name.contains('XMM-LSS')) |
                (
                    # include extra constraints on XMM-LSS targets
                    (t.mi < self.parameters['mag_i_max_xmm_lss'] )
                )
            ) &
            (t.specz >= self.parameters['specz_min']) &
            (t.specz <= self.parameters['specz_max'])
        )

        return query




class BhmRmVarCarton(BhmRmBaseCarton):
    '''bhm_rm_var: selected based on g-band variability > 0.05 mag and bright enough to be detected by Gaia (G<~21)

    SELECT * FROM bhm_rm
        WHERE ( (des_var_sn[0] > 5.0 AND des_var_rms[0] > 0.05)  OR
                (ps1_var_sn[0]>5.0 AND ps1_var_rms[0]>0.05))
        AND WHERE mi BETWEEN 15.0 AND 21.5    # <- exact limits TBD
        AND WHERE pmsig < x.x
        AND WHERE plxsig < x.x
        AND WHERE gaia = 1

    #debug select t.pk, t.ra, t.dec, t.mi, t.psfmag_sdss[4] as
    psfmag_i,t.pmsig,t.ps1_var_sn[1],t.ps1_var_rms[1],t.des_var_sn[1],t.des_var_rms[1]
    from bhm_rm_v0 as t where (t.gaia = 1 AND t.mi < 21.5 AND t.pmsig
    < 5.0 AND t.plxsig < 5.0 AND t.ps1_var_sn[1] > 5.0 AND
    t.ps1_var_rms[1] > 0.05 ) limit 10;

    '''

    name = 'bhm_rm_var'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t
        query = query.where(
            (
                (
                    (t.des_var_sn[0] > self.parameters['des_var_sn_min']) &
                    (t.des_var_rms[0] > self.parameters['des_var_rms_min'])
                ) |
                (
                    (t.ps1_var_sn[0] > self.parameters['ps1_var_sn_min']) &
                    (t.ps1_var_rms[0] > self.parameters['ps1_var_rms_min'])
                )
            ) &
            (t.source_id_gaia > 0 ) &
            (t.pmsig < self.parameters['pmsig_max']) &
            (t.plxsig < self.parameters['plxsig_max']) &
            (t.pmsig > pmsig_min) &  # this catches cases with NULL=-9
            (t.plxsig > plxsig_min) &
            ~(t.field_name.contains('SDSS-RM')) # ignore this carton in the SDSS-RM field
        )
            #& (t.mg <  self.parameters['g_mag_max'])  # TBD

        return query



class BhmRmAncillaryCarton(BhmRmBaseCarton):
    '''
    bhm_rm_ancillary: from the Gaia_unWISE AGN catalog or the XDQSO catalog,
                      but requiring no proper motion/parallax detection from Gaia DR2

    SELECT * FROM bhm_rm
        WHERE photo_bitmask & (2^0+2^1) != 0
        AND WHERE mi BETWEEN 15.0 AND 21.5    ← TBD
        AND WHERE pmsig < 5.0
        AND WHERE plxsig < 5.0
    '''

    name = 'bhm_rm_ancillary'

    def build_query(self, version_id, query_region=None):
        query = super().build_query(version_id, query_region)
        t = self.alias_t

        query = query.where(
            (t.photo_bitmask.bin_and(self.parameters['photo_bitmask']) != 0 ) &
            (t.pmsig <  self.parameters['pmsig_max']) &
            (t.plxsig <  self.parameters['plxsig_max']) &
            (t.pmsig > pmsig_min) &  # this catches cases with NULL=-9
            (t.plxsig > plxsig_min) &
            ~(t.field_name.contains('SDSS-RM')) # ignore this carton in the SDSS-RM field
        )

        return query










#######################################################################
#######################################################################
#######################################################################
#######################################################################
# Notes and debug
#
# '''
# t = BHM_RM_v0.alias()
# for f in t._meta.fields:
#     print (f)
# '''
#
# sql, param = qq.sql()
# print (sql.replace("%s","{}").format(*param))
#    # some functional SQL for one field:
#    ##SELECT c.*,t.specz,t.mi FROM catalog AS C INNER JOIN catalog_to_bhm_rm_v0 AS c2t ON c.catalogid = c2t.catalogid INNER JOIN bhm_rm_v0 AS t ON c2t.target_id = t.pk WHERE (q3c_radial_query(c.ra, c.dec, 36.45, -4.6, 1.49) AND c.version_id = 13 AND t.specz > 0.005 AND mi < 21.1 AND mi > 16.0);
#    SELECT DISTINCT ON (t.pk) c.*,t.mg,t.mi,t.specz,t.spec_q,t.spec_strmask,specz_ref,q3c_dist(c.ra,c.dec,t.ra,t.dec)*3600.::FLOAT as sep FROM catalog AS C INNER JOIN catalog_to_bhm_rm_v0 AS c2t ON c.catalogid = c2t.catalogid INNER JOIN bhm_rm_v0 AS t ON c2t.target_id = t.pk WHERE (q3c_radial_query(c.ra, c.dec, 36.45, -4.6, 1.49) AND c.version_id = 13 AND t.specz > 0.005 AND mi < 21.1 AND mi > 16.0);
#    \copy (SELECT DISTINCT ON (t.pk) c.*,t.mg,t.mi,t.specz,t.spec_q,t.spec_strmask,specz_ref,q3c_dist(c.ra,c.dec,t.ra,t.dec)*3600.::FLOAT as sep FROM catalog AS C INNER JOIN catalog_to_bhm_rm_v0 AS c2t ON c.catalogid = c2t.catalogid INNER JOIN bhm_rm_v0 AS t ON c2t.target_id = t.pk WHERE (q3c_radial_query(c.ra, c.dec, 35.5, -4.6, 1.49) AND c.version_id = 13 AND t.specz > 0.005 AND mi < 22.0 AND mi > 16.0)) to ./test_faint.csv csv header;
#



'''
# for testing do domething like this

import peewee
import sdssdb
from sdssdb.peewee.sdss5db.catalogdb import database
database.set_profile('tunnel_operations')
from target_selection.cartons.bhm_rm import *
c = BhmRmCoreCarton(targeting_plan='0.1.0-beta.1')
c = BhmRmKnownSpecCarton(targeting_plan='0.1.0-beta.1')
c = BhmRmVarCarton(targeting_plan='0.1.0-beta.1')
c = BhmRmAncillaryCarton(targeting_plan='0.1.0-beta.1')
q = c.build_query(version_id=13)
q.count()
for r in q.limit(5).namedtuples():
    print(r)

'''


'''
#target_selection  --profile tunnel_operations --verbose run --include bhm_rm_known_spec --keep --overwrite '0.1.0-beta.1' --no-load

target_selection  --profile tunnel_operations --verbose run --include bhm_rm_core,bhm_rm_var,bhm_rm_ancillary,bhm_rm_known_spec --keep --overwrite '0.1.0-beta.1' --no-load

# to get a psql prompt:
#> psql -d sdss5db -U sdss_user -h localhost -p 7502

# Exporting from the temp table
# in psql terminal:

\copy (SELECT * FROM sandbox.temp_bhm_rm_core)  TO '/home/tdwelly/scratch/targetdb/bhm_rm_core.csv' with csv header
\copy (SELECT * FROM sandbox.temp_bhm_rm_known_spec)  TO '/home/tdwelly/scratch/targetdb/bhm_rm_known_spec.csv' with csv header
\copy (SELECT * FROM sandbox.temp_bhm_rm_var)  TO '/home/tdwelly/scratch/targetdb/bhm_rm_var.csv' with csv header
\copy (SELECT * FROM sandbox.temp_bhm_rm_ancillary)  TO '/home/tdwelly/scratch/targetdb/bhm_rm_ancillary.csv' with csv header

for F in bhm_rm_*.csv; do   stilts tpipe in=${F} out="${F%.*}.fits" ifmt=csv ofmt=fits-basic; done

head -n 1 bhm_rm_core.csv  > bhm_rm.csv
tail -q -n +2 bhm_rm_*.csv  >> bhm_rm.csv
stilts tpipe in=bhm_rm.csv out=bhm_rm.fits ifmt=csv ofmt=fits-basic
ftcopy "bhm_rm.fits[1][col *,bhm_rm_core(L)=priority==1002?1:0,bhm_rm_known_spec(L)=priority==1001?1:0,bhm_rm_var(L)=priority==1003?1:0,bhm_rm_ancillary(L)=priority==1004?1:0]" bhm_rm.fits clobber=yes mode=q
ftsort bhm_rm.fits bhm_rm_unique.fits catalogid method=heap unique=yes clobber=yes mode=q

gnuplot plot_rm_target_stats.plot
tar -cvzf bhm_core_cartons.tar.gz bhm_rm*.fits bhm_aqmes_*fits bhm_spiders_*fits bhm_csc_*fits

'''

#######################################################################
#######################################################################
#######################################################################
