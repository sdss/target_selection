#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Zach Way (zway1@gsu.edu)
# @Date: 2024-04-23
# @Filename: append_to_sdss_id.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import os
import peewee
from peewee import fn, JOIN

from sdssdb.peewee.sdss5db.catalogdb import database, Catalog
from create_catalogidx_to_catalogidy import MetaXMatch, create_unique_from_region
from create_catalogidx_to_catalogidy import TempMatch, UniqueMatch

database.connect(dbname="sdss5db", user="sdss_user")


class temp_catalogid_v21(peewee.Model):
    """Model for a temporary table of all catalogids within a version."""
    pk = peewee.PrimaryKeyField()
    catalogid21 = peewee.BigIntegerField(index=True, null=False)

    class Meta:
        database = database
        schema = "sandbox"
        table_name = "temp_catalogid_v21"


class temp_catalogid_v25(peewee.Model):
    """Model for a temporary table of all catalogids within a version."""
    pk = peewee.PrimaryKeyField()
    catalogid25 = peewee.BigIntegerField(index=True, null=False)

    class Meta:
        database = database
        schema = "sandbox"
        table_name = "temp_catalogid_v25"


class temp_catalogid_v31(peewee.Model):
    """Model for a temporary table of all catalogids within a version."""
    pk = peewee.PrimaryKeyField()
    catalogid31 = peewee.BigIntegerField(index=True, null=False)

    class Meta:
        database = database
        schema = "sandbox"
        table_name = "temp_catalogid_v31"


class sdss_id_stacked(peewee.Model):
    """ Model for catalogdb sdss_id_stacked"""

    sdss_id = peewee.BigAutoField(primary_key=True)
    catalogid21 = peewee.BigIntegerField()
    catalogid25 = peewee.BigIntegerField()
    catalogid31 = peewee.BigIntegerField()
    ra_sdss_id = peewee.DoubleField()
    dec_sdss_id = peewee.DoubleField()

    class Meta:
        database = database
        schema = "catalogdb"  # "sandbox"
        table_name = "sdss_id_stacked"


class sdss_id_flat(peewee.Model):
    """ Model for catalogdb sdss_id_flat"""

    sdss_id = peewee.BigIntegerField()
    catalogid = peewee.BigIntegerField()
    version_id = peewee.SmallIntegerField()
    ra_sdss_id = peewee.DoubleField()
    dec_sdss_id = peewee.DoubleField()
    n_associated = peewee.SmallIntegerField()
    ra_catalogid = peewee.DoubleField()
    dec_catalogid = peewee.DoubleField()
    pk = peewee.BigAutoField(primary_key=True)

    class Meta:
        database = database
        schema = "catalogdb"  # "sandbox"
        table_name = "sdss_id_flat"


class sdss_id_stacked_to_add(peewee.Model):
    """ Model for addendum to sdss_id_stacked"""

    catalogid21 = peewee.BigIntegerField(index=True)
    catalogid25 = peewee.BigIntegerField(index=True)
    catalogid31 = peewee.BigIntegerField(index=True)
    ra_sdss_id = peewee.DoubleField()
    dec_sdss_id = peewee.DoubleField()

    class Meta:
        database = database
        schema = "sandbox"


class sdss_id_flat_to_add(peewee.Model):
    """ Model for rows to be added to sdss_id_flat """

    sdss_id = peewee.BigIntegerField()
    catalogid = peewee.BigIntegerField()
    version_id = peewee.SmallIntegerField()
    ra_sdss_id = peewee.DoubleField()
    dec_sdss_id = peewee.DoubleField()
    n_associated = peewee.SmallIntegerField(null=True)
    ra_catalogid = peewee.DoubleField(null=True)
    dec_catalogid = peewee.DoubleField(null=True)

    class Meta:
        database = database
        schema = "sandbox"
        table_name = "sdss_id_flat_to_add"


class append_to_tables:
    """foobar"""

    def __init__(self, database, individual_table=None, catalogid_list=None):
        self.database = database
        self.individual_table = individual_table
        self.ind_table_clean = self.individual_table
        self.catalogid_list = catalogid_list
        self.dir_path = os.path.dirname(__file__)

        if self.catalogid_list is not None:
            config = {"version_ids_to_match": [21, 25, 31],
                      "individual_xmatch_config":
                      self.dir_path + "/config/individual_crossmatches.yml",
                      "log_file": "catalogidx_to_catalogidy_from_list.log",
                      "show_first": 20,
                      "split_insert_nunmber": 100000,
                      "database_options": {"enable_hashjoin": "false"},
                      # "split_query": [['panstarrs1',522000000000000,5000000000000]],
                      "catalogid_list": catalogid_list}

        else:
            if "." in self.individual_table:
                self.ind_table_clean = self.individual_table.split(".")[-1]

            config = {"version_ids_to_match": [21, 25, 31],
                      "individual_xmatch_config":
                      self.dir_path + "/config/individual_crossmatches.yml",
                      "log_file": f"catalogidx_to_catalogidy_{self.ind_table_clean}.log",
                      "show_first": 20,
                      "split_insert_nunmber": 100000,
                      "database_options": {"enable_hashjoin": "false"},
                      "split_query": [['panstarrs1', 522000000000000, 5000000000000]],
                      "individual_table": individual_table}
        self.config = config

    def run_MetaXMatch(self, database):
        """ This takes in an individual catalog_to_? table (labeled 'individual_table') or a
         list of catalogids and creates the catalogidx_to_catalogidy_?
         and catalogidx_to_catalogidy_?_unique tables.
        """
        metax = MetaXMatch(config_filename=None, database=database, from_yaml=False,
                           from_dict=True, config_dict=self.config)

        metax.run()
        create_unique_from_region(metax.output_name)
        return metax.output_name

    def create_temp_catalogid_lists(self, database, output_name):
        """ This method takes the newly created catalogidx_to_catalogidy table and
        creates three aggregate tables in the sandbox containing all of the possible
        catalogids to be matched. This is also where the process restricts the matching of
        catalogids to 1) those which are labeled a 'best' match and 2) catalogids which have an
        inter-crossmatch match on that catalogid's lead catalog """

        TempMatch._meta.table_name = output_name
        UniqueMatch._meta.table_name = output_name + "_unique"

        v21_cid_query_x = (TempMatch
                           .select(TempMatch.catalogidx.alias('catalogid21'))
                           .join(Catalog, on=((Catalog.catalogid == TempMatch.catalogidx) &
                                              (Catalog.lead == TempMatch.lead)))
                           .switch(TempMatch)
                           .join(sdss_id_stacked, join_type=JOIN.LEFT_OUTER,
                                 on=(TempMatch.catalogidx == sdss_id_stacked.catalogid21))
                           .where((TempMatch.version_idx == 21) &
                                  (sdss_id_stacked.catalogid21.is_null())))

        v25_cid_query_x = (TempMatch
                           .select(TempMatch.catalogidx.alias('catalogid25'))
                           .join(Catalog, on=((Catalog.catalogid == TempMatch.catalogidx) &
                                              (Catalog.lead == TempMatch.lead)))
                           .switch(TempMatch)
                           .join(sdss_id_stacked, join_type=JOIN.LEFT_OUTER,
                                 on=(TempMatch.catalogidx == sdss_id_stacked.catalogid25))
                           .where((TempMatch.version_idx == 25) &
                                  (sdss_id_stacked.catalogid25.is_null())))

        v31_cid_query_x = (TempMatch
                           .select(TempMatch.catalogidx.alias('catalogid31'))
                           .join(Catalog, on=((Catalog.catalogid == TempMatch.catalogidx) &
                                              (Catalog.lead == TempMatch.lead)))
                           .switch(TempMatch)
                           .join(sdss_id_stacked, join_type=JOIN.LEFT_OUTER,
                                 on=(TempMatch.catalogidx == sdss_id_stacked.catalogid31))
                           .where((TempMatch.version_idx == 31) &
                                  (sdss_id_stacked.catalogid31.is_null())))

        v21_cid_query_y = (TempMatch
                           .select(TempMatch.catalogidy.alias('catalogid21'))
                           .join(Catalog, on=((Catalog.catalogid == TempMatch.catalogidy) &
                                              (Catalog.lead == TempMatch.lead)))
                           .switch(TempMatch)
                           .join(sdss_id_stacked, join_type=JOIN.LEFT_OUTER,
                                 on=(TempMatch.catalogidy == sdss_id_stacked.catalogid21))
                           .where((TempMatch.version_idy == 21) &
                                  (sdss_id_stacked.catalogid21.is_null())))

        v25_cid_query_y = (TempMatch
                           .select(TempMatch.catalogidy.alias('catalogid25'))
                           .join(Catalog, on=((Catalog.catalogid == TempMatch.catalogidy) &
                                              (Catalog.lead == TempMatch.lead)))
                           .switch(TempMatch)
                           .join(sdss_id_stacked, join_type=JOIN.LEFT_OUTER,
                                 on=(TempMatch.catalogidy == sdss_id_stacked.catalogid25))
                           .where((TempMatch.version_idy == 25) &
                                  (sdss_id_stacked.catalogid25.is_null())))

        v31_cid_query_y = (TempMatch
                           .select(TempMatch.catalogidy.alias('catalogid31'))
                           .join(Catalog, on=((Catalog.catalogid == TempMatch.catalogidy) &
                                              (Catalog.lead == TempMatch.lead)))
                           .switch(TempMatch)
                           .join(sdss_id_stacked, join_type=JOIN.LEFT_OUTER,
                                 on=(TempMatch.catalogidy == sdss_id_stacked.catalogid31))
                           .where((TempMatch.version_idy == 31) &
                                  (sdss_id_stacked.catalogid31.is_null())))

        for table in [temp_catalogid_v21, temp_catalogid_v25, temp_catalogid_v31]:
            if table.table_exists():
                database.drop_tables([table])

        self.database.create_tables([temp_catalogid_v21, temp_catalogid_v25, temp_catalogid_v31])

        temp_catalogid_v21.insert_from(v21_cid_query_x, [temp_catalogid_v21.catalogid21]).execute()
        temp_catalogid_v21.insert_from(v21_cid_query_y, [temp_catalogid_v21.catalogid21]).execute()
        temp_catalogid_v25.insert_from(v25_cid_query_x, [temp_catalogid_v25.catalogid25]).execute()
        temp_catalogid_v25.insert_from(v25_cid_query_y, [temp_catalogid_v25.catalogid25]).execute()
        temp_catalogid_v31.insert_from(v31_cid_query_x, [temp_catalogid_v31.catalogid31]).execute()
        temp_catalogid_v31.insert_from(v31_cid_query_y, [temp_catalogid_v31.catalogid31]).execute()

        if self.catalogid_list is not None:
            catid_input_query = (Catalog.select(Catalog.catalogid, Catalog.version_id)
                                 .join(sdss_id_flat, join_type=JOIN.LEFT_OUTER,
                                       on=(sdss_id_flat.catalogid == Catalog.catalogid))
                                 .where(Catalog.catalogid << self.catalogid_list)
                                 .where(sdss_id_flat.catalogid.is_null()))

            for row in catid_input_query:
                if row.version_id == 21:
                    temp_catalogid_v21.insert(catalogid21=row.catalogid).execute()
                elif row.version_id == 25:
                    temp_catalogid_v25.insert(catalogid25=row.catalogid).execute()
                elif row.version_id == 31:
                    temp_catalogid_v31.insert(catalogid31=row.catalogid).execute()
        else:
            try:
                ind_table = self.database.models[self.individual_table]
            except:
                raise ValueError("Could not find the table in database model: " +
                                 self.individual_table)

            ind_table_input_query = (Catalog.select(Catalog.catalogid, Catalog.version_id)
                                     .join(ind_table,
                                           on=(Catalog.catalogid == ind_table.catalogid))
                                     .join(sdss_id_flat, join_type=JOIN.LEFT_OUTER,
                                           on=(sdss_id_flat.catalogid == Catalog.catalogid))
                                     .where(sdss_id_flat.catalogid.is_null()))

            for row in ind_table_input_query:
                if row.version_id == 21:
                    temp_catalogid_v21.insert(catalogid21=row.catalogid).execute()
                elif row.version_id == 25:
                    temp_catalogid_v25.insert(catalogid25=row.catalogid).execute()
                elif row.version_id == 31:
                    temp_catalogid_v31.insert(catalogid31=row.catalogid).execute()

    def create_sdss_id_stacked_to_add(self, database, output_name):
        """ This method matches the previous catalogids with the 'sdss_id' process
        and places them in a temporary table in the sandbox (to be added later). At
        this point, the ra/dec_sdss_id are assigned, which are just the coordinates of
        the most recent catalogid in the match. """

        large_cte_query = f"""DROP TABLE sandbox.sdss_id_stacked_to_add;
        CREATE TABLE sandbox.sdss_id_stacked_to_add as (
            with cat_ids21 as (select distinct tc.catalogid21, cat.lead from
                sandbox.temp_catalogid_v21 tc join catalog cat on cat.catalogid=tc.catalogid21),
            cat_ids25 as (select distinct tc.catalogid25, cat.lead from
                sandbox.temp_catalogid_v25 tc join catalog cat on cat.catalogid=tc.catalogid25),
            cat_ids31 as (select distinct tc.catalogid31, cat.lead from
                sandbox.temp_catalogid_v31 tc join catalog cat on cat.catalogid=tc.catalogid31),
            sq21_25 as (select cc.catalogidx, cc.catalogidy from sandbox.{output_name} cc
                    join catalog cat on cat.catalogid=cc.catalogidx and cat.lead=cc.lead
                    where cc.version_idx=21 and cc.version_idy=25),
            sq25_31 as (select cc.catalogidx, cc.catalogidy from sandbox.{output_name} cc
                    join catalog cat on cat.catalogid=cc.catalogidx and cat.lead=cc.lead
                    where cc.version_idx=25 and cc.version_idy=31),
            left21_to_25 as (select cat_ids21.catalogid21, sq21_25.catalogidy as catalogid25 from
                cat_ids21 left join sq21_25 on cat_ids21.catalogid21=sq21_25.catalogidx),
            add_outer25 as (select left21_to_25.catalogid21,
                cat_ids25.catalogid25 as catalogid25 from
                left21_to_25 full outer join cat_ids25 on
                left21_to_25.catalogid25=cat_ids25.catalogid25),
            left_25_to_31 as (select add_outer25.catalogid21,
                add_outer25.catalogid25, sq25_31.catalogidy as catalogid31 from
                add_outer25 left join sq25_31 on add_outer25.catalogid25=sq25_31.catalogidx),
            add_outer31 as (select left_25_to_31.catalogid21, left_25_to_31.catalogid25,
                cat_ids31.catalogid31 from
                left_25_to_31 full outer join cat_ids31 on
                left_25_to_31.catalogid31=cat_ids31.catalogid31)
            select * from add_outer31);"""
        self.database.execute_sql(large_cte_query)

        add_ra_dec_columns = """ ALTER TABLE sandbox.sdss_id_stacked_to_add
                                    ADD COLUMN ra_sdss_id double precision;
                                 ALTER TABLE sandbox.sdss_id_stacked_to_add
                                    ADD COLUMN dec_sdss_id double precision; """
        self.database.execute_sql(add_ra_dec_columns)
        sdss_id_stacked_to_add._meta.table_name = "sdss_id_stacked_to_add"

        ra_dec_update31 = (sdss_id_stacked_to_add
                           .update(ra_sdss_id=Catalog.ra, dec_sdss_id=Catalog.dec)
                           .from_(Catalog)
                           .where(sdss_id_stacked_to_add.catalogid31 == Catalog.catalogid)
                           .where(sdss_id_stacked_to_add.catalogid31.is_null(False)))

        ra_dec_update25 = (sdss_id_stacked_to_add
                           .update(ra_sdss_id=Catalog.ra, dec_sdss_id=Catalog.dec)
                           .from_(Catalog)
                           .where(sdss_id_stacked_to_add.catalogid25 == Catalog.catalogid)
                           .where(sdss_id_stacked_to_add.catalogid31.is_null())
                           .where(sdss_id_stacked_to_add.catalogid25.is_null(False)))

        ra_dec_update21 = (sdss_id_stacked_to_add
                           .update(ra_sdss_id=Catalog.ra, dec_sdss_id=Catalog.dec)
                           .from_(Catalog)
                           .where(sdss_id_stacked_to_add.catalogid21 == Catalog.catalogid)
                           .where(sdss_id_stacked_to_add.catalogid31.is_null())
                           .where(sdss_id_stacked_to_add.catalogid25.is_null())
                           .where(sdss_id_stacked_to_add.catalogid21.is_null(False)))

        ra_dec_update31.execute()
        ra_dec_update25.execute()
        ra_dec_update21.execute()

    def add_to_sdss_id_stacked(self, database):
        """ This method adds the previously made sdss_id_stached_to_add to sdss_id_stacked """

        sid_stacked_f = [sdss_id_stacked.catalogid21,
                         sdss_id_stacked.catalogid25,
                         sdss_id_stacked.catalogid31,
                         sdss_id_stacked.ra_sdss_id,
                         sdss_id_stacked.dec_sdss_id]

        to_add = sdss_id_stacked_to_add.select(sdss_id_stacked_to_add.catalogid21,
                                               sdss_id_stacked_to_add.catalogid25,
                                               sdss_id_stacked_to_add.catalogid31,
                                               sdss_id_stacked_to_add.ra_sdss_id,
                                               sdss_id_stacked_to_add.dec_sdss_id).tuples()
        with self.database.atomic():
            sdss_id_stacked.insert_many(to_add, fields=sid_stacked_f).execute()

    def create_sdss_id_flat_to_add(self, database):
        """ This method takes the new sdss_id rows in sdss_id_stacked and
        adds them to a temporary table in the sandbox in the sdss_id_flat format.
        Here, the ra/dec_catalogid are populated with the specific values associated
        with the catalogid. """

        if sdss_id_flat_to_add.table_exists():
            self.database.drop_tables([sdss_id_flat_to_add])
            self.database.create_tables([sdss_id_flat_to_add])
        else:
            self.database.create_tables([sdss_id_flat_to_add])

        sid_flat_fields = [sdss_id_flat_to_add.sdss_id,
                           sdss_id_flat_to_add.catalogid,
                           sdss_id_flat_to_add.version_id,
                           sdss_id_flat_to_add.ra_sdss_id,
                           sdss_id_flat_to_add.dec_sdss_id]

        max_sdss_id = sdss_id_flat.select(fn.MAX(sdss_id_flat.sdss_id)).scalar()

        sid_flat_add_v21 = (sdss_id_stacked.select(sdss_id_stacked.sdss_id,
                                                   sdss_id_stacked.catalogid21,
                                                   peewee.Value(21),
                                                   sdss_id_stacked.ra_sdss_id,
                                                   sdss_id_stacked.dec_sdss_id)
                            .where(sdss_id_stacked.catalogid21.is_null(False))
                            .where(sdss_id_stacked.sdss_id > max_sdss_id).tuples())

        with self.database.atomic():
            sdss_id_flat_to_add.insert_many(sid_flat_add_v21, fields=sid_flat_fields).execute()

        sid_flat_add_v25 = (sdss_id_stacked.select(sdss_id_stacked.sdss_id,
                                                   sdss_id_stacked.catalogid25,
                                                   peewee.Value(25),
                                                   sdss_id_stacked.ra_sdss_id,
                                                   sdss_id_stacked.dec_sdss_id)
                            .where(sdss_id_stacked.catalogid25.is_null(False))
                            .where(sdss_id_stacked.sdss_id > max_sdss_id).tuples())

        with self.database.atomic():
            sdss_id_flat_to_add.insert_many(sid_flat_add_v25, fields=sid_flat_fields).execute()

        sid_flat_add_v31 = (sdss_id_stacked.select(sdss_id_stacked.sdss_id,
                                                   sdss_id_stacked.catalogid31,
                                                   peewee.Value(31),
                                                   sdss_id_stacked.ra_sdss_id,
                                                   sdss_id_stacked.dec_sdss_id)
                            .where(sdss_id_stacked.catalogid31.is_null(False))
                            .where(sdss_id_stacked.sdss_id > max_sdss_id).tuples())
        with self.database.atomic():
            sdss_id_flat_to_add.insert_many(sid_flat_add_v31, fields=sid_flat_fields).execute()

        cte = (sdss_id_flat_to_add
               .select(sdss_id_flat_to_add.catalogid, fn.COUNT("*").alias('ct'))
               .group_by(sdss_id_flat_to_add.catalogid)
               .cte('catid_n_associated', columns=("catalogid", "ct")))

        (sdss_id_flat_to_add.update(n_associated=cte.c.ct)
                            .from_(cte)
                            .where(sdss_id_flat_to_add.catalogid == cte.c.catalogid)
                            .with_cte(cte)
                            .execute())

        (sdss_id_flat_to_add.update(ra_catalogid=Catalog.ra,
                                    dec_catalogid=Catalog.dec)
                            .from_(Catalog)
                            .where(sdss_id_flat_to_add.catalogid == Catalog.catalogid)
                            .execute())

    def add_to_sdss_id_flat(self, database):
        """ This method adds sdss_id_flat_to_add to sdss_id_flat """

        sid_flat_f = [sdss_id_flat.sdss_id,
                      sdss_id_flat.catalogid,
                      sdss_id_flat.version_id,
                      sdss_id_flat.ra_sdss_id,
                      sdss_id_flat.dec_sdss_id,
                      sdss_id_flat.n_associated,
                      sdss_id_flat.ra_catalogid,
                      sdss_id_flat.dec_catalogid]

        to_add = sdss_id_flat_to_add.select(sdss_id_flat_to_add.sdss_id,
                                            sdss_id_flat_to_add.catalogid,
                                            sdss_id_flat_to_add.version_id,
                                            sdss_id_flat_to_add.ra_sdss_id,
                                            sdss_id_flat_to_add.dec_sdss_id,
                                            sdss_id_flat_to_add.n_associated,
                                            sdss_id_flat_to_add.ra_catalogid,
                                            sdss_id_flat_to_add.dec_catalogid).tuples()

        with self.database.atomic():
            sdss_id_flat.insert_many(to_add, fields=sid_flat_f).execute()


"""
EXAMPLE OF USAGE
foo = append_to_tables(database, catalogid_list=[list, of, integer, catalogids])
output_name = foo.run_MetaXMatch(database)
foo.create_temp_catalogid_lists(database, output_name)
foo.create_sdss_id_stacked_to_add(database, output_name)
foo.add_to_sdss_id_stacked(database)
foo.create_sdss_id_flat_to_add(database)
foo.add_to_sdss_id_flat(database)
"""
