import os
import time

import peewee
import yaml
from peewee import Case, Model, fn

from sdssdb.peewee.sdss5db.catalogdb import database

import target_selection


schema = 'sandbox'


class dup1(Model):
    """Model for the first output table."""

    catalogid = peewee.BigIntegerField(null=False)
    target_id = peewee.BigIntegerField(null=False)
    version_id = peewee.SmallIntegerField(null=False)
    distance = peewee.DoubleField(null=True)
    best = peewee.BooleanField(null=False)
    dup_catalogid = peewee.BooleanField(null=False)
    dup_targetid = peewee.BooleanField(null=False)

    class Meta:
        table_name = 'v1_duplicates'
        database = database
        schema = schema
        primary_key = False


class dup1text(dup1):
    """ Model for tables that have text primary keys like tycho2"""
    catalogid = peewee.BigIntegerField(null=False)
    target_id = peewee.TextField(null=False)

    class Meta:
        table_name = 'v1_duplicates'


class BaseOutput(Model):
    """Base Model for the second and final output tables."""

    input_table = peewee.TextField(null=False)
    catalogid = peewee.BigIntegerField(null=False)
    target_id = peewee.BigIntegerField(null=False)
    version_id = peewee.SmallIntegerField(null=False)
    distance = peewee.DoubleField(null=True)
    best = peewee.BooleanField(null=False)
    dup_catalogid = peewee.BooleanField(null=False)
    dup_targetid = peewee.BooleanField(null=False)
    dup_type_catid = peewee.TextField(null=False)
    dup_type_tarid = peewee.TextField(null=False)

    class Meta:
        database = database
        schema = schema
        primary_key = False


class BaseOutputText(Model):
    """Base Model for the second and final output tables for tables with text pks."""

    input_table = peewee.TextField(null=False)
    catalogid = peewee.BigIntegerField(null=False)
    target_id = peewee.TextField(null=False)
    version_id = peewee.SmallIntegerField(null=False)
    distance = peewee.DoubleField(null=True)
    best = peewee.BooleanField(null=False)
    dup_catalogid = peewee.BooleanField(null=False)
    dup_targetid = peewee.BooleanField(null=False)
    dup_type_catid = peewee.TextField(null=False)
    dup_type_tarid = peewee.TextField(null=False)

    class Meta:
        database = database
        schema = schema
        primary_key = False


class dup2(BaseOutput):

    class Meta:
        table_name = 'v1_duplicates2'


class dup2text(dup2):
    """ Model for tables that have text primary keys like tycho2"""
    input_table = peewee.TextField(null=False)
    catalogid = peewee.BigIntegerField(null=False)
    target_id = peewee.TextField(null=False)

    class Meta:
        table_name = 'v1_duplicates2'


class dup_final(BaseOutput):
    """ Model for final output table"""
    class Meta:
        table_name = 'v1_duplicates_info'


class dup_final_text(BaseOutputText):
    """ Model for tables that have text primary keys like tycho2"""

    class Meta:
        table_name = 'v1_duplicates_info_text'


class Duplicates:
    """Calculates duplicated targetids and catalogids in relational tables.

    This class calculates the duplicated target_id and catalogid values in the v1 entries of the
    relational tables that were added as part of the v1 crossmatch (although it should be easy
    to modify to use for another crossmatch run).

    For each table crossmatched in v1 it looks all the entries in the relational table were
    either the catalogid or the target_id is repeated in the table. Then it stores all those
    duplicated entries from each relational table in an intermediate result table called
    v1_duplicates that contain the information from the relational table along with the
    two booleans ``dup_catalogid`` and ``dup_targetid`` to indicate if the catalogid and/or
    the target_id is duplicates.

    Then, the second query add the text columns ``dup_type_catid`` and ``dup_type_tarid``.
    For duplicates originated in phase_1 where we dont have distance information to compare
    the duplicates, there columns get value ``phase_1``. For duplicates originated in phase_2 of
    spatial matching the distance information is used to assign value ``best`` or ``not best``.
    If there are multiple target_ids associated with the same catalogid the closest target_id
    will get dup_type_catid=best while all the others will get dup_type_catid=not best.
    On the other hand if there are multiple catalogids associated with the same target_id,
    the closest catalogid will get dup_type_tarid=best. Finally entries where the catalogid or
    target_id is not duplicated will get value ``None`` in dup_type_catid or dup_type_tarid
    respectively.

    To correlate this information with the existing ``best`` column in the relational table,
    entries with best=True correspond to entries with dup_type_catid=best and a fraction of
    the stars with dup_type_catid=phase_1 (the ones that correspond to the lowest target_id
    linked).

    Thus, these results provide two advantages over the existing information. 1 they
    discriminate the entries with best=True were this actually means that is the closest match
    with the ones where the value got assigned in phase_1 and thus has no meaning, and 2 it
    provides the information in the other direction for the cases where one target_id is
    matched to multiple catalogids which actually happens more often because we ingest first
    the higher resolution tables that often resolve stars in the lower resolution tables
    ingested afterwards.

    Final results are stored in a table called sandbox.v1_duplicates_info. However, there are
    some relational tables that have text content in target_id instead of numbers, like tycho2,
    so those results are stored in a separate table called sandbox.v1_duplicates_info_text

    Parameters
    ----------
    config_filename : str
        Configuration file used for the crossmatch run to analyze.

    plan: str
        crossmatch plan to analyze corresponding to the plan value in
        table catalogdb.version

    database : peewee:PostgresqlDatabase
        A `peewee:PostgresqlDatabase` to the database the tables to
        cross-match.

    log_filename : str
        Name of log file to use. Optional, and if it is not set it uses
        duplicates_from_xmatch_plan_(plan).log
    """

    def __init__(self, config_filename, plan, database, log_filename=None):
        if log_filename is None:
            log_filename = f'./duplicates_from_xmatch_plan_{plan}.log'
        self.database = database
        self.plan = plan
        config = yaml.load(open(config_filename, 'r'), Loader=yaml.SafeLoader)[plan]
        self.log = target_selection.log
        log_path_and_name = os.path.realpath('./' + log_filename)
        self.log.start_file_logger(log_path_and_name, rotating=False, mode='a')
        self.log.sh.setLevel(0)

        version_model = database.models['catalogdb.version']
        query_version = version_model.select().where(version_model.plan == plan)
        self.version_id = query_version.dicts()[0]['id']
        self.xmatched_tables = config['order']

    def run_query1(self, model, rel_model):

        ti = time.time()

        # If the target_id is a text column then we need to use a special model
        if rel_model._meta.fields['target_id'].field_type == 'TEXT':
            dup1m = dup1text
        else:
            dup1m = dup1

        cte_catalogid = (rel_model
                         .select(rel_model.catalogid, fn.COUNT(rel_model.target_id))
                         .where(rel_model.version_id == self.version_id)
                         .group_by(rel_model.catalogid)
                         .having(fn.COUNT(rel_model.target_id) > 1)
                         .cte('cte_dupcat')
                         )

        cte_targetid = (rel_model
                        .select(rel_model.target_id, fn.COUNT(rel_model.catalogid))
                        .where(rel_model.version_id == self.version_id)
                        .group_by(rel_model.target_id)
                        .having(fn.COUNT(rel_model.catalogid) > 1)
                        .cte('cte_duptar')
                        )

        query1 = (rel_model
                  .select(rel_model.catalogid, rel_model.target_id,
                          rel_model.version_id, rel_model.distance, rel_model.best,
                          Case(None,
                               [(cte_catalogid.c.catalogid.is_null(), False)],
                               True)
                          .alias('dup_targetid'),
                          Case(None,
                               [(cte_targetid.c.target_id.is_null(), False)],
                               True)
                          .alias('dup_catalogid'))
                  .join(cte_catalogid, 'LEFT OUTER JOIN',
                        on=(rel_model.catalogid == cte_catalogid.c.catalogid))
                  .join(cte_targetid, 'LEFT OUTER JOIN',
                        on=(rel_model.target_id == cte_targetid.c.target_id))
                  .where(rel_model.version_id == self.version_id,
                         (cte_catalogid.c.catalogid.is_null(False)
                          | cte_targetid.c.target_id.is_null(False)))
                  .with_cte(cte_catalogid, cte_targetid)
                  )

        if dup1m.table_exists():
            database.drop_tables([dup1m])
            self.log.info(f'Dropped table {dup1m._meta.table_name}')
        dup1m.create_table()
        self.log.info(f'Creating table {dup1m._meta.table_name}')
        dup1_fields = [dup1m._meta.fields[field] for field in dup1m._meta.fields]

        n_ins1 = dup1m.insert_from(query1, dup1_fields).returning().execute()
        tf = time.time()
        self.log.info(f'The first query took {(tf-ti):.2f} seconds and selected {n_ins1:,} '
                      'duplicates')

    def run_query2(self, model, rel_model):

        ti = time.time()

        # If the target_id is a text column then we need to use a special model
        if rel_model._meta.fields['target_id'].field_type == 'TEXT':
            dup1m = dup1text
            dup2m = dup2text
        else:
            dup1m = dup1
            dup2m = dup2

        tablename = model._meta.table_name
        q2 = (dup1m.
              select(peewee.Value(tablename), dup1m.catalogid, dup1m.target_id, dup1m.version_id,
                     dup1m.distance, dup1m.best, dup1m.dup_catalogid, dup1m.dup_targetid,
                     Case(None, [((dup1m.dup_catalogid >> True) & (dup1m.distance.is_null()),
                                  'phase_1'),
                                 (((dup1m.dup_catalogid >> True)
                                   & (dup1m.distance.is_null(False)
                                   & (fn.first_value(dup1m.distance)
                                      .over(partition_by=[dup1m.catalogid],
                                            order_by=[dup1m.distance.asc()]) == dup1m.distance))),
                                  'best'),
                                 (((dup1m.dup_catalogid >> True)
                                   & (dup1m.distance.is_null(False)
                                   & (fn.first_value(dup1m.distance)
                                      .over(partition_by=[dup1m.catalogid],
                                            order_by=[dup1m.distance.asc()]) != dup1m.distance))),
                                  'not best')],
                          'None').alias('dup_type_catid'),
                     Case(None, [((dup1m.dup_targetid >> True) & (dup1m.distance.is_null()),
                                  'phase_1'),
                                 (((dup1m.dup_targetid >> True)
                                   & (dup1m.distance.is_null(False)
                                   & (fn.first_value(dup1m.distance)
                                      .over(partition_by=[dup1m.target_id],
                                            order_by=[dup1m.distance.asc()]) == dup1m.distance))),
                                  'best'),
                                 (((dup1m.dup_targetid >> True)
                                   & (dup1m.distance.is_null(False)
                                   & (fn.first_value(dup1m.distance)
                                      .over(partition_by=[dup1m.target_id],
                                            order_by=[dup1m.distance.asc()]) != dup1m.distance))),
                                  'not best')],
                          'None').alias('dup_type_tarid'))
              )

        if dup2m.table_exists():
            database.drop_tables([dup2m])
            self.log.info(f'Dropped table {dup2m._meta.table_name}')
        dup2m.create_table()
        self.log.info(f'Creating table {dup2m._meta.table_name}')
        dup2_fields = [dup2m._meta.fields[field] for field in dup2m._meta.fields]

        dup2m.insert_from(q2, dup2_fields).execute()

        tf = time.time()
        self.log.info(f'The second query took {(tf-ti):.2f} seconds')

    def process_table(self, tablename):

        """
        This method runs the 2 queries for a given table and stores the results
        in the final table.
        """
        ti = time.time()
        models = self.database.models
        model = models[f'catalogdb.{tablename}']
        rel_model = models[f'catalogdb.catalog_to_{tablename}']
        self.log.info(f'Processing table {tablename}')
        self.run_query1(model, rel_model)
        self.run_query2(model, rel_model)

        # If the target_id is a text column then we need to use a special model
        if rel_model._meta.fields['target_id'].field_type == 'TEXT':
            dup_final_model = dup_final_text
            dup2m = dup2text
        else:
            dup_final_model = dup_final
            dup2m = dup2

        dup_final_fields = [dup_final_model._meta.fields[field]
                            for field in dup_final_model._meta.fields]
        ins_query = (dup_final_model.insert_from(dup2m
                                                 .select(dup2m.input_table, dup2m.catalogid,
                                                         dup2m.target_id, dup2m.version_id,
                                                         dup2m.distance, dup2m.best,
                                                         dup2m.dup_catalogid, dup2m.dup_targetid,
                                                         dup2m.dup_type_catid,
                                                         dup2m.dup_type_tarid),
                                                 dup_final_fields)
                     .returning())
        n_ins = ins_query.execute()
        tf = time.time()
        self.log.info(f'Completely processed {tablename} in {(tf-ti):.2f} seconds '
                      f'ingesting {n_ins:,} entries into {dup_final._meta.table_name}')
        return n_ins

    def run(self, starting_table=None, ending_table=None, overwrite=True):

        """ This method goes over the different tables to store the results of each
        in the final output table. If overwrite is True (default) it creates a new output
        table deleting the old one if it exists. By default it loops over all the tables
        crossmatched in the corresponding plan but starting_table and ending_table
        can be used to restrict the list of tables in case the code needs to be tested
        or restarted.
        """

        self.log.info(f'Calculating duplicates ids for Plan: {self.plan}')
        self.log.info(f'Corresponding to version id {self.version_id}')
        total_entries = 0
        ti = time.time()
        tables_to_process = self.xmatched_tables
        if ending_table is not None:
            tables_to_process = tables_to_process[:tables_to_process.index(ending_table) + 1]
        if starting_table is not None:
            tables_to_process = tables_to_process[tables_to_process.index(starting_table):]

        rel_models = [self.database.models[f'catalogdb.catalog_to_{name}']
                      for name in tables_to_process]
        pk_types = [rel_model._meta.fields['target_id'].field_type for rel_model in rel_models]

        if overwrite is True:
            if 'INT' in pk_types or 'BIGINT' in pk_types:  # Processing tables with numbers in PK
                if dup_final.table_exists():
                    database.drop_tables([dup_final])
                    self.log.info(f'Dropped table {dup_final._meta.table_name}')
                dup_final.create_table()
                self.log.info(f'Creating table {dup_final._meta.table_name}')

            if 'TEXT' in pk_types:  # This means we are processing tables with text in PK
                if dup_final_text.table_exists():
                    database.drop_tables([dup_final_text])
                    self.log.info(f'Dropped table {dup_final_text._meta.table_name}')
                dup_final_text.create_table()
                self.log.info(f'Creating table {dup_final_text._meta.table_name}')

        for tablename in tables_to_process:
            n_ins = self.process_table(tablename)
            total_entries += n_ins
        tf = time.time()
        self.log.info(' ')
        self.log.info(f'Finished checking duplicates ingesting {total_entries:,} total entries '
                      f'in {(tf-ti):.2f} seconds')
        self.log.info(' ')


"""
EXAMPLE OF USAGE
dupobj = Duplicates('xmatch_v1.yml', '1.0.0', database)
dupobj.run(starting_table='tycho2', ending_table='tycho2', overwrite=False)
"""
