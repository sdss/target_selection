import os
import time
from itertools import combinations

import numpy as np
import peewee
import yaml
from peewee import fn, JOIN
from playhouse.postgres_ext import ArrayField

from sdssdb.peewee.sdss5db.catalogdb import database, Catalog
from sdssdb.peewee.sdss5db.targetdb import Target

import target_selection

from create_catalogidx_to_catalogidy import MetaXMatch, create_unique_from_region, TempMatch, UniqueMatch

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
        schema = "catalogdb"
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
        schema = "catalogdb"
        table_name = "sdss_id_flat"

class sdss_id_stacked_to_add(peewee.Model):
    """ Model for addendum to sdss_id_stacked"""

    #sdss_id = peewee.BigAutoField(index=True, primary_key=True)
    catalogid21 = peewee.BigIntegerField(index=True)
    catalogid25 = peewee.BigIntegerField(index=True)
    catalogid31 = peewee.BigIntegerField(index=True)
    ra_sdss_id = peewee.DoubleField()
    dec_sdss_id = peewee.DoubleField()

    class Meta:
        database = database
        schema = "sandbox"

class sdss_id_stacked_new(peewee.Model):
    """ Model for new sdss_id_stacked """

    sdss_id = peewee.BigAutoField(primary_key=True)
    catalogid21 = peewee.BigIntegerField()
    catalogid25 = peewee.BigIntegerField()
    catalogid31 = peewee.BigIntegerField()
    ra_sdss_id = peewee.DoubleField()
    dec_sdss_id = peewee.DoubleField()

    class Meta:
        database = database
        schema = "sandbox"
        table_name = "sdss_id_stacked_new"

class sdss_id_flat_new(peewee.Model):
    """ Model for new sdss_id_flat"""

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
        schema = "sandbox"
        table_name = "sdss_id_flat_new"

#class sdss_id_flat_to_add(peewee.Model):
#    """ Model for addendum to sdss_id_flat"""
#
#    sdss_id = peewee.BigIntegerField(index=True)
#    catalogid = peewee.BigIntegerField(index=True)
#    version_id = peewee.SmallIntegerField(index=True)
#    ra_sdss_id = peewee.DoubleField()
#    dec_sdss_id = peewee.DoubleField()
#    n_associated = peewee.SmallIntegerField(index=True)
#    ra_catalogid = peewee.DoubleField()
#    dec_catalogid = peewee.DoubleField()
#    pk = peewee.BigIntegerField(primary_key=True)
#
#    class Meta:
#        database = database
#        schema = "sandbox"


def duplicate_sdss_id_stacked():
    """ Method to duplicate sdss_id_stacked to the
    sandbox.
    """
    
    #fields = [f.name for f in sdss_id_stacked._meta.fields]
    
    with database:
        database.execute_sql('CREATE TABLE IF NOT EXISTS sandbox.sdss_id_stacked_new AS SELECT * FROM catalogdb.sdss_id_stacked')
        database.execute_sql('ALTER TABLE sandbox.sdss_id_stacked_new ALTER COLUMN sdss_id set not null')
        database.execute_sql('CREATE SEQUENCE sandbox.sdss_id_stacked_new_sdss_id__seq')
        database.execute_sql('ALTER TABLE sandbox.sdss_id_stacked_new ALTER COLUMN sdss_id SET DEFAULT nextval(\'sandbox.sdss_id_stacked_new_sdss_id__seq\')')
        database.execute_sql('ALTER SEQUENCE sandbox.sdss_id_stacked_new_sdss_id__seq OWNED BY sandbox.sdss_id_stacked_new.sdss_id')
        database.execute_sql('SELECT setval(\'sandbox.sdss_id_stacked_new_sdss_id__seq\', (SELECT MAX(sdss_id) FROM sandbox.sdss_id_stacked_new))')

class append_to_tables:
    """foobar"""

    def __init__(self, individual_table, database):
        self.database = database
        self.individual_table = individual_table
        self.ind_table_clean = self.individual_table
        if "." in self.individual_table:
            self.ind_table_clean = self.individual_table.split(".")[-1]

        config = {"version_ids_to_match" : [21,25,31],
                         "individual_xmatch_config" : "individual_crossmatches.yml",
                         "log_file" : f"catalogidx_to_catalogidy_{self.ind_table_clean}.log",
                         "show_first" : 20,
                         "split_insert_nunmber" : 100000,
                         "database_options" : {"enable_hashjoin" : "false"},
                         "split_query" : [['panstarrs1',522000000000000,5000000000000]],
                         "individual_table" : individual_table}
    
    def run_MetaXMatch(self, database, config_dict):
        """ This takes in an individual catalog_to_? table (labeled 'individual_table') and
        creates the catalogidx_to_catalogidy_? and catalogidx_to_catalogidy_?_unique tables
        in the sandox.

        """
        metax = MetaXMatch(config_filename=None, database=database, from_yaml=False, config_dict=config_dict)
    
        metax.run()
        create_unique_from_region(metax.output_name)
        return metax.output_name
    

    def create_temp_catalogid_lists(self, database, output_name):
        """ foobar """
    
        TempMatch._meta.table_name = output_name
        UniqueMatch._meta.table_name = output_name+"_unique"
        
        v21_cid_query_x = (TempMatch
                             .select(TempMatch.catalogidx.alias('catalogid21'))
                             .join(Catalog, on=((Catalog.catalogid==TempMatch.catalogidx)&(Catalog.lead==TempMatch.lead)))
                             .switch(TempMatch)
                             .join(sdss_id_stacked, join_type = JOIN.LEFT_OUTER, 
                                   on=(TempMatch.catalogidx==sdss_id_stacked.catalogid21))
                             .where((TempMatch.version_idx==21)&(sdss_id_stacked.catalogid21.is_null())))
        
        v25_cid_query_x = (TempMatch
                             .select(TempMatch.catalogidx.alias('catalogid25'))
                             .join(Catalog, on=((Catalog.catalogid==TempMatch.catalogidx)&(Catalog.lead==TempMatch.lead)))
                             .switch(TempMatch)
                             .join(sdss_id_stacked, join_type = JOIN.LEFT_OUTER, 
                                   on=(TempMatch.catalogidx==sdss_id_stacked.catalogid25))
                             .where((TempMatch.version_idx==25)&(sdss_id_stacked.catalogid25.is_null())))
        
        v31_cid_query_x = (TempMatch
                             .select(TempMatch.catalogidx.alias('catalogid31'))
                             .join(Catalog, on=((Catalog.catalogid==TempMatch.catalogidx)&(Catalog.lead==TempMatch.lead)))
                             .switch(TempMatch)
                             .join(sdss_id_stacked, join_type = JOIN.LEFT_OUTER,
                                   on=(TempMatch.catalogidx==sdss_id_stacked.catalogid31))
                             .where((TempMatch.version_idx==31)&(sdss_id_stacked.catalogid31.is_null())))

        v21_cid_query_y = (TempMatch
                             .select(TempMatch.catalogidy.alias('catalogid21'))
                             .join(Catalog, on=((Catalog.catalogid==TempMatch.catalogidy)&(Catalog.lead==TempMatch.lead)))
                             .switch(TempMatch)
                             .join(sdss_id_stacked, join_type = JOIN.LEFT_OUTER,
                                   on=(TempMatch.catalogidy==sdss_id_stacked.catalogid21))
                             .where((TempMatch.version_idy==21)&(sdss_id_stacked.catalogid21.is_null())))

        v25_cid_query_y = (TempMatch
                             .select(TempMatch.catalogidy.alias('catalogid25'))
                             .join(Catalog, on=((Catalog.catalogid==TempMatch.catalogidy)&(Catalog.lead==TempMatch.lead)))
                             .switch(TempMatch)
                             .join(sdss_id_stacked, join_type = JOIN.LEFT_OUTER,
                                   on=(TempMatch.catalogidy==sdss_id_stacked.catalogid25))
                             .where((TempMatch.version_idy==25)&(sdss_id_stacked.catalogid25.is_null())))

        v31_cid_query_y = (TempMatch
                             .select(TempMatch.catalogidy.alias('catalogid31'))
                             .join(Catalog, on=((Catalog.catalogid==TempMatch.catalogidy)&(Catalog.lead==TempMatch.lead)))
                             .switch(TempMatch)
                             .join(sdss_id_stacked, join_type = JOIN.LEFT_OUTER,
                                   on=(TempMatch.catalogidy==sdss_id_stacked.catalogid31))
                             .where((TempMatch.version_idy==31)&(sdss_id_stacked.catalogid31.is_null())))

        for table in [temp_catalogid_v21, temp_catalogid_v25, temp_catalogid_v31]:
            if table.table_exists():
                self.database.drop_tables([table])
            
        temp_catalogid_v21.insert_from(v21_cid_query_x, [temp_catalogid_v21.catalogid21]).execute()
        temp_catalogid_v21.insert_from(v21_cid_query_y, [temp_catalogid_v21.catalogid21]).execute()    
        temp_catalogid_v25.insert_from(v25_cid_query_x, [temp_catalogid_v25.catalogid25]).execute()
        temp_catalogid_v25.insert_from(v25_cid_query_y, [temp_catalogid_v25.catalogid25]).execute()
        temp_catalogid_v31.insert_from(v31_cid_query_x, [temp_catalogid_v31.catalogid31]).execute()
        temp_catalogid_v31.insert_from(v31_cid_query_y, [temp_catalogid_v31.catalogid31]).execute()
        
    def create_sdss_id_stacked_to_add(self, database):
        """ foobar """

        cat_ids21 = (temp_catalogid_v21.select(temp_catalogid_v21.catalogid21, Catalog.lead).distinct()
                                         .join(Catalog, on=(Catalog.catalogid==temp_catalogid_v21.catalogid21)))
        cat_ids25 = (temp_catalogid_v25.select(temp_catalogid_v25.catalogid25, Catalog.lead).distinct()
                                         .join(Catalog, on=(Catalog.catalogid==temp_catalogid_v25.catalogid25)))
        cat_ids31 = (temp_catalogid_v31.select(temp_catalogid_v31.catalogid31, Catalog.lead).distinct()
                                         .join(Catalog, on=(Catalog.catalogid==temp_catalogid_v31.catalogid31)))
        
        sq21_25 = (TempMatch.select(TempMatch.catalogidx, TempMatch.catalogidy)
                                .join(Catalog, on=((Catalog.catalogid==TempMatch.catalogidx)&(Catalog.lead==TempMatch.lead)))
                                .where((TempMatch.version_idx==21)&(TempMatch.version_idy==25)))
        sq25_31 = (TempMatch.select(TempMatch.catalogidx, TempMatch.catalogidy)
                                .join(Catalog, on=((Catalog.catalogid==TempMatch.catalogidx)&(Catalog.lead==TempMatch.lead)))
                                .where((TempMatch.version_idx==25)&(TempMatch.version_idy==31)))

        left21_to_25 = (cat_ids21.select(cat_ids21.catalogid21, sq21_25.catalogidy.alias("catalogid25"))
                                .join(sq21_25, join_type = JOIN.LEFT_OUTER, on=(cat_ids21.catalogid21==sq21_25.catalogidx))
                                .with_cte(cat_ids21, sq21_25))
        
        add_outer25 = (left21_to_25.select(left21_to_25.catalogid21, cat_ids25.catalogid25.alias("catalogid25"))
                                .join(cat_ids25, join_type = JOIN.FULL_OUTER, on=(left21_to_25.catalogid25==cat_ids25.catalogid25))
                                .with_cte(cat_ids25, left21_to_25))
        
        left_25_to_31 = (add_outer25.select(add_outer25.catalogid21, add_outer25.catalogid25, sq25_31.catalogidy.alias("catalogid31"))
                                .join(sq25_31, join_typye = JOIN.LEFT_OUTER, on=(add_outer25.catalogid25==sq25_31.catalogidx))
                                .with_cte(add_outer25, sq25_31))
        
        add_outer31 = (left_25_to_31.select(left_25_to_31.catalogid21, left_25_to_31.catalogid25, cat_ids31.catalogid31)
                                .join(cat_ids31, join_type=JOIN.FULL_OUTER, on=(left_25_to_31.catalogid31==cat_ids31.catalogid31)))

        sdss_id_stacked_to_add._meta.table_name = "sdss_id_stacked_to_add" 
        insert_all = (sdss_id_stacked_to_add.insert_from(add_outer31, 
                                    [sdss_id_stacked_to_add.catalogid21, 
                                     sdss_id_stacked_to_add.catalogid25, 
                                     sdss_id_stacked_to_add.catalogi31])
                                .with_cte(add_outer31))
        insert_all.execute()
        
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
        """ foobar """
        
        if sdss_id_stacked_new.table_exists():
            self.database.drop_tables([sdss_id_stacked_new])
            duplicate_sdss_id_stacked()
        
        sid_stacked_to_add_f = [sdss_id_stacked_to_add.catalogid21, 
                                sdss_id_stacked_to_add.catalogid25, 
                                sdss_id_stacked_to_add.catalogid31,
                                sdss_id_stacked_to_add.ra_sdss_id,
                                sdss_id_stacked_to_add.dec_sdss_id]
        sid_stacked_new_f = [sdss_id_stacked_new.catalogid21,
                                sdss_id_stacked_new.catalogid25,
                                sdss_id_stacked_new.catalogid31,
                                sdss_id_stacked_new.ra_sdss_id,
                                sdss_id_stacked_new.dec_sdss_id]
        insert_stacked_to_add = (sdss_id_stacked_new.
                                     .insert_from(sdss_id_stacked_to_add.select(sid_stacked_to_add_f),
                                         fields=sid_stacked_new_f)
                                     .execute())
        
    def create_sdss_id_flat_new(self, database):
        """ foobar """
        
        if sdss_id_flat_new.table_exists():
            self.database.drop_tables([sdss_id_flat_new])
            self.database.create_tables([sdss_id_flat_new])

        query = (sdss_id_flat_new
                     .insert_from(sdss_id_stacked_new.select([sdss_id_stacked_new.sdss_id,
                                                              sdss_id_stacked_new.catalogid21.alias('catalogid'),
                                                              peewee.Value(21).alias('version_id'),
                                                              sdss_id_stacked_new.ra_sdss_id,
                                                              sdss_id_stacked_new.dec_sdss_id])
                                                     .where(sdss_id_stacked_new.catalogid21.is_null(False)))
                     .execute())
        
        query = (sdss_id_flat_new
                     .insert_from(sdss_id_stacked_new.select([sdss_id_stacked_new.sdss_id,
                                                              sdss_id_stacked_new.catalogid25.alias('catalogid'),
                                                              peewee.Value(25).alias('version_id'),
                                                              sdss_id_stacked_new.ra_sdss_id,
                                                              sdss_id_stacked_new.dec_sdss_id])
                                                     .where(sdss_id_stacked_new.catalogid25.is_null(False)))
                     .execute())
        
        query = (sdss_id_flat_new
                     .insert_from(sdss_id_stacked_new.select([sdss_id_stacked_new.sdss_id,
                                                              sdss_id_stacked_new.catalogid31.alias('catalogid'),
                                                              peewee.Value(31).alias('version_id'),
                                                              sdss_id_stacked_new.ra_sdss_id,
                                                              sdss_id_stacked_new.dec_sdss_id])
                                                     .where(sdss_id_stacked_new.catalogid31.is_null(False)))
                     .execute())

        with_cte = (sdss_id_flat_new
                        .select(sdss_id_flat_new.catalogid, fn.COUNT().alias('ct'))
                        .group_by(sdss_id_flat_new.catalogid))

        query = (sdss_id_flat_new
                     .update(sdss_id_flat_new.n_associated=with_cte.ct)
                     .from_(with_cte)
                     .where(sdss_id_flat_new.catalogid == with_cte.c.catalogid)
                     .execute())

