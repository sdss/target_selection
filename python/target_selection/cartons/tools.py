#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2021-04-29
# @Filename: tools.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import os
import tempfile

import numpy
import peewee
from astropy.table import Table

from sdssdb.utils.ingest import create_model_from_table

from target_selection.cartons import BaseCarton
from target_selection.exceptions import TargetSelectionError
from target_selection.utils import vacuum_table


def get_file_carton(filename):
    """Returns a carton class that creates a carton based on a FITS file.

    The FITS file is located in the ``open_fiber_path`` which is specified in
    ``python/config/target_selection.yml``.

    The list of FITS files to be loaded is specified in the
    file ``open_fiber_file_list.txt`` which is in the directory ``open_fiber_path``.

    """

    # Import this here to prevent this module not being importable if the database
    # connection is not ready.
    from sdssdb.peewee.sdss5db.catalogdb import (
        Catalog,
        CatalogToGaia_DR2,
        CatalogToGaia_DR3,
        CatalogToLegacy_Survey_DR8,
        CatalogToLegacy_Survey_DR10,
        CatalogToPanstarrs1,
        CatalogToTwoMassPSC,
        Gaia_DR2,
        Gaia_DR3,
        Legacy_Survey_DR8,
        Legacy_Survey_DR10,
        Panstarrs1,
        TwoMassPSC,
    )

    class FileCarton(BaseCarton):
        can_offset = None  # Will be set in query.

        def __init__(self, targeting_plan, config_file=None, schema=None, table_name=None):
            self._file_path = filename

            self._table = Table.read(self._file_path)
            if self._table.masked:
                self._table = self._table.filled()

            self._run_sanity_checks()

            super().__init__(
                targeting_plan,
                config_file=config_file,
                schema=schema,
                table_name=table_name,
            )

            self._disable_query_log = True

        def _run_sanity_checks(self):
            """Runs a series of sanity checks on the FITS table."""

            unique_cartonname = numpy.unique(self._table["cartonname"])
            if len(unique_cartonname) == 1:
                self.name = unique_cartonname[0].lower()
            else:
                raise TargetSelectionError(
                    "Error in get_file_carton(): "
                    + filename
                    + " contains more than one cartonname"
                )

            unique_can_offset = numpy.unique(self._table["can_offset"])
            if len(unique_can_offset) > 1:
                raise TargetSelectionError(
                    "Error in get_file_carton(): "
                    + filename
                    + " contains more than one"
                    + " value of can_offset:"
                    + " can_offset values must be "
                    + " all 0 or all 1"
                )

            if (unique_can_offset[0] != 1) and (unique_can_offset[0] != 0):
                raise TargetSelectionError(
                    "Error in get_file_carton(): "
                    + filename
                    + " can_offset can only be 0 or 1."
                    + " can_offset is "
                    + str(unique_can_offset[0])
                )

            unique_inertial = numpy.unique(self._table["inertial"])
            if len(unique_inertial) > 2:
                raise TargetSelectionError(
                    "Error in get_file_carton(): "
                    + filename
                    + " contains more than two"
                    + " values of inertial:"
                    + " inertial values must be "
                    + " 0 or 1"
                )

            if (unique_inertial[0] != 1) and (unique_inertial[0] != 0):
                raise TargetSelectionError(
                    "Error in get_file_carton(): "
                    + filename
                    + " inertial can only be 0 or 1."
                    + " inertial is "
                    + str(unique_inertial[0])
                )

            # If there is only one inertial value then the above statement
            # is enough. Otherwise, we need to run the below check.
            if len(unique_inertial) == 2:
                if (unique_inertial[1] != 1) and (unique_inertial[1] != 0):
                    raise TargetSelectionError(
                        "Error in get_file_carton(): "
                        + filename
                        + " inertial can only be 0 or 1."
                        + " inertial is "
                        + str(unique_inertial[1])
                    )

            # The valid_program list is from the output of the below command.
            # select distinct(program) from targetdb.carton order by program;
            #
            # mwm_bin is for future mwm binary star cartons
            valid_program = [
                "bhm_aqmes",
                "bhm_csc",
                "bhm_filler",
                "bhm_rm",
                "bhm_spiders",
                "commissioning",
                "mwm_bin",
                "mwm_cb",
                "mwm_dust",
                "mwm_erosita",
                "mwm_filler",
                "mwm_galactic",
                "mwm_gg",
                "mwm_halo",
                "mwm_legacy",
                "mwm_magcloud",
                "mwm_ob",
                "mwm_planet",
                "mwm_rv",
                "mwm_snc",
                "mwm_tessob",
                "mwm_tessrgb",
                "mwm_validation",
                "mwm_wd",
                "mwm_yso",
                "open_fiber",
                "ops",
                "ops_sky",
                "ops_std",
                "SKY",
            ]

            # The valid_category list is from CartonImportTable.pdf
            valid_category = [
                "science",
                "standard_apogee",
                "standard_boss",
                "guide",
                "sky_boss",
                "sky_apogee",
                "standard",
                "sky",
                "veto location boss",
                "veto_location_apogee",
            ]

            # The valid_category list is from CartonImportTable.pdf
            valid_mapper = ["", "MWM", "BHM"]

            unique_category = numpy.unique(self._table["category"])
            if len(unique_category) == 1:
                self.category = unique_category[0].lower()
                if self.category not in valid_category:
                    raise TargetSelectionError(
                        "Error in get_file_carton(): "
                        + filename
                        + " contains invalid category = "
                        + self.category
                    )
            else:
                raise TargetSelectionError(
                    "Error in get_file_carton(): " + filename + " contains more than one category"
                )

            unique_program = numpy.unique(self._table["program"])
            if len(unique_program) == 1:
                self.program = unique_program[0].lower()
                if self.program not in valid_program:
                    raise TargetSelectionError(
                        "Error in get_file_carton(): "
                        + filename
                        + " contains invalid program = "
                        + self.program
                    )
            else:
                raise TargetSelectionError(
                    "Error in get_file_carton(): " + filename + " contains more than one program"
                )

            unique_mapper = numpy.unique(self._table["mapper"])
            if len(unique_mapper) == 1:
                # We do not use lower() for mapper because
                # allowed values for mapper are '' or 'MWM' or 'BHM'.
                self.mapper = unique_mapper[0]
                if self.mapper not in valid_mapper:
                    raise TargetSelectionError(
                        "Error in get_file_carton(): "
                        + filename
                        + " contains invalid mapper = "
                        + self.mapper
                    )
                if self.mapper == "":
                    self.mapper = None
            else:
                raise TargetSelectionError(
                    "Error in get_file_carton(): " + filename + " contains more than one mapper"
                )

            basename_fits = os.path.basename(filename)
            basename_parts = os.path.splitext(basename_fits)
            basename = basename_parts[0]
            carton_name_from_filename = basename.lower()

            if self.name != carton_name_from_filename:
                raise TargetSelectionError(
                    "filename parameter of get_file_carton() and "
                    + "cartonname in FITS file do not match."
                    + "\n"
                    + "carton_name_from_filename = "
                    + carton_name_from_filename
                    + " cartonname = "
                    + self.name
                )

        def copy_data(self, temp_table: str):
            """Copy the input file data to a temporary table.

            The schema of the file carton table is such that we can dump to
            a CSV file without issues.

            """

            temp_csv = tempfile.NamedTemporaryFile(suffix=".csv", delete=False)
            self._table.write(temp_csv.name, format="csv", overwrite=True)

            cursor = self.database.cursor()
            temp_csv.seek(0)
            cursor.copy_expert(f"COPY {temp_table} FROM STDOUT WITH CSV HEADER", temp_csv)
            self.database.commit()

        def build_query(self, version_id, query_region=None):
            self.log.debug(f"Processing file {self._file_path}.")

            # We need to copy the data to a temporary table so that we can
            # join on it. We could use a Peewee ValueList but for large tables
            # that will hit the limit of 1GB in PSQL.

            # Create model for temporary table from FITS table columns.
            # This works fine because we know there are no arrays.
            temp_table = f"inputs_{self.name.lower()}_temp"
            self.database.execute_sql(f"DROP TABLE IF EXISTS {temp_table};")
            temp = create_model_from_table(temp_table, self._table)
            temp._meta.database = self.database
            temp.create_table(temporary=True)

            # Copy the data
            self.copy_data(temp_table)

            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("Gaia_DR3_Source_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("Gaia_DR2_Source_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("LegacySurvey_DR10_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("LegacySurvey_DR8_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("PanSTARRS_DR2_ID")')
            self.database.execute_sql(f'CREATE INDEX ON "{temp_table}" ("TwoMASS_ID")')
            vacuum_table(self.database, temp_table, vacuum=False, analyze=True)

            inertial_case = peewee.Case(
                None,
                ((temp.inertial.cast("boolean").is_null(), False),),
                temp.inertial.cast("boolean"),
            )

            # List of columns and aliases for the final query table.
            query_columns = [
                Catalog.catalogid,
                temp.Gaia_DR3_Source_ID.alias("gaia_dr3_source_id"),
                temp.Gaia_DR2_Source_ID.alias("gaia_source_id"),
                temp.LegacySurvey_DR10_ID.alias("ls_id10"),
                temp.LegacySurvey_DR8_ID.alias("ls_id8"),
                temp.PanSTARRS_DR2_ID.alias("catid_objid"),
                temp.TwoMASS_ID.alias("designation"),
                Catalog.ra,
                Catalog.dec,
                temp.delta_ra.cast("double precision"),
                temp.delta_dec.cast("double precision"),
                inertial_case.alias("inertial"),
                temp.cadence,
                temp.priority,
                temp.instrument,
                temp.can_offset.cast("boolean").alias("can_offset"),
                peewee.Value(0).alias("value"),
            ]

            # Calculate number of rows in the table for each parent catalogue identifier and
            # run some sanity checks.
            len_table = len(self._table)
            len_gaia_dr3 = len(self._table[self._table["Gaia_DR3_Source_ID"] > 0])
            len_gaia_dr2 = len(self._table[self._table["Gaia_DR2_Source_ID"] > 0])
            len_legacysurvey_dr10 = len(self._table[self._table["LegacySurvey_DR10_ID"] > 0])
            len_legacysurvey_dr8 = len(self._table[self._table["LegacySurvey_DR8_ID"] > 0])
            len_panstarrs_dr2 = len(self._table[self._table["PanSTARRS_DR2_ID"] > 0])

            # TwoMass_ID corresponds to the designation column of
            # the table catalogdb.twomass_psc.
            # Since the designation column is a text column, below
            # we are comparing it to the string 'NA' and not the integer 0.
            len_twomass_psc = len(self._table[self._table["TwoMASS_ID"] != "NA"])

            # Make sure this is not an empty table.
            if len_table == 0:
                raise TargetSelectionError(
                    f"Error in get_file_carton(): {self._file_path} is an empty table"
                )

            # There must be exactly one non-zero id per row else raise an exception.
            if (
                len_gaia_dr3
                + len_gaia_dr2
                + len_legacysurvey_dr10
                + len_legacysurvey_dr8
                + len_panstarrs_dr2
                + len_twomass_psc
            ) != len_table:
                raise TargetSelectionError(
                    "Error in get_file_carton(): "
                    + "(len_gaia_dr3 + len_gaia_dr2 + "
                    + "len_legacysurvey_dr10 + len_legacysurvey_dr8 +"
                    + "len_panstarrs_dr2 + len_twomass_psc) != "
                    + "len_table"
                )

            # For each identifier that has non-zero targets, appends some information that
            # we need to create the subquery for that table. This includes the CatalogToX model
            # to which we need to join, the field on which to join in the parent catalogue,
            # and the column in the temporary table. In all cases except 2MASS, we join on the
            # parent catalogue primary key.
            model_data = []

            if len_gaia_dr3 > 0:
                model_data.append(
                    {
                        "catalog_to": CatalogToGaia_DR3,
                        "parent_field": Gaia_DR3.source_id,
                        "temp_column": "Gaia_DR3_Source_ID",
                    }
                )

            if len_gaia_dr2 > 0:
                model_data.append(
                    {
                        "catalog_to": CatalogToGaia_DR2,
                        "parent_field": Gaia_DR2.source_id,
                        "temp_column": "Gaia_DR2_Source_ID",
                    }
                )

            if len_legacysurvey_dr10 > 0:
                model_data.append(
                    {
                        "catalog_to": CatalogToLegacy_Survey_DR10,
                        "parent_field": Legacy_Survey_DR10.ls_id,
                        "temp_column": "LegacySurvey_DR10_ID",
                    }
                )

            if len_legacysurvey_dr8 > 0:
                model_data.append(
                    {
                        "catalog_to": CatalogToLegacy_Survey_DR8,
                        "parent_field": Legacy_Survey_DR8.ls_id,
                        "temp_column": "LegacySurvey_DR8_ID",
                    }
                )

            if len_panstarrs_dr2 > 0:
                model_data.append(
                    {
                        "catalog_to": CatalogToPanstarrs1,
                        "parent_field": Panstarrs1.catid_objid,
                        "temp_column": "PanSTARRS_DR2_ID",
                    }
                )

            if len_twomass_psc > 0:
                model_data.append(
                    {
                        "catalog_to": CatalogToTwoMassPSC,
                        "parent_field": TwoMassPSC.designation,
                        "temp_column": "TwoMASS_ID",
                    }
                )

            if len(model_data) == 0:
                raise TargetSelectionError(
                    "Error in get_file_carton(): no join model found for the file carton."
                )

            # Create a query for each join model. The final query will be the union of all.
            # For each join model, we need to account for cases in which an identifier is
            # associated with more than one catalogid via a window function (we select either
            # the phase 1 match, or the one with the smallest distance).

            queries = []

            for data in model_data:
                catalog_to_model = data["catalog_to"]
                parent_field = data["parent_field"]
                temp_column = data["temp_column"]

                # Get the model class field for the column in the temporary table with the
                # identifier for this case.
                temp_field = getattr(temp, temp_column)

                # Get the psrent model. The only reason why we need to join all the way to the
                # parent catalogue model is 2MASS for which the column TwoMASS_ID in the file
                # carton corresponds to the designation column in the TwoMassPSC model, which
                # is not the primary key.
                parent_model = parent_field.model

                # Create a subquery that ranks the rows by distance to the target. Since
                # temp is also used in the main query, we need to alias it.
                # Note that we are using ROW_NUMBER() and not RANK() because the latter would
                # assign the same rank to multiple rows with the same distance. This can happen
                # if two catalogids are associated with a target in phase 1 (an example of this
                # is a Gaia DR2 target that has been deblended into two Gaia DR3 targets).
                # We also order by catalogid to ensure that the query is deterministic but
                # ultimately we are randomly selecting one of the targets here.
                temp_alias = temp.alias("temp_alias")
                temp_alias_field = getattr(temp_alias, temp_column)

                distance_rank_partition = peewee.fn.row_number().over(
                    partition_by=[catalog_to_model.target_id],
                    order_by=[
                        peewee.fn.coalesce(catalog_to_model.distance, 0.0).asc(),
                        catalog_to_model.catalogid.asc(),
                    ],
                )

                sub_query = (
                    temp_alias.select(
                        catalog_to_model.catalogid,
                        parent_field.alias("target_id"),
                        distance_rank_partition.alias("distance_rank"),
                    )
                    .join(parent_model, on=(temp_alias_field == parent_field))
                    .join(catalog_to_model)
                    .where(
                        catalog_to_model.best >> True,
                        catalog_to_model.version_id == version_id,
                    )
                ).alias("distance_rank_subquery")

                # Now add the main query. We join the subquery to the temporary table to
                # get all the relevant columns, but keep only the entries with distance_rank=1.
                queries.append(
                    Catalog.select(*query_columns, sub_query.c.distance_rank)
                    .join(
                        sub_query,
                        on=(Catalog.catalogid == sub_query.c.catalogid),
                    )
                    .join(
                        temp,
                        on=(temp_field == sub_query.c.target_id),
                    )
                    .where(sub_query.c.distance_rank == 1)
                    .distinct(Catalog.catalogid)
                )

            # Union all queries.
            query_union = queries[0]
            for query in queries[1:]:
                query_union = query_union | query

            # It seems to work better if we disable sequential scans and force the use of the
            # indices.
            self.database.execute_sql("SET LOCAL enable_seqscan = off;")

            # Now just distinct on catalogid for all the unions. Although we already have a
            # distinct in each query, they can yield the same catalogid from different queries.
            query_union = query_union.cte("query_union")

            return (
                query_union.select(query_union.__star__)
                .distinct(query_union.c.catalogid)
                .with_cte(query_union)
            )

        def post_process(self, model, **kwargs):
            """Runs sanity checks on the output of the query."""

            n_file_table = len(self._table)
            n_query = model.select().count()

            if n_file_table != n_query:
                self.log.warning(
                    f"The number of rows in the file table ({n_file_table}) does not "
                    f"match the number of rows returned by the query ({n_query})."
                )

    return FileCarton


def create_table_as(
    query,
    table_name,
    schema=None,
    temporary=False,
    database=None,
    execute=True,
    overwrite=False,
    indices=[],
    analyze=True,
):
    """Creates a table from a query.

    Parameters
    ----------
    query
        A Peewee ``ModelSelect`` or a string with the query to create a table from.
    table_name
        The name of the table to create.
    schema
        The schema in which to create the table. If ``table_name`` is in the
        form ``schema.table``, the schema parameter is overridden by ``table_name``.
    temporary
        Whether to create a temporary table instead of a persistent one.
    database
        The database connection to use to execute the query. If not passed and the
        query is a ``ModelSelect``, the database will be inherited from the query
        model.
    execute
        Whether to actually execute the query. Requires ``database`` to be passed.
    overwrite
        If `True`, the table will be create even if a table with the same name
        already exists. Requires ``database`` or will be ignored.
    analyze
        Whether to ``VACUUM ANALIZE`` the new table.
        Only relevant if ``execute=True``.
    indices
        List of columns to create indices on. Only relevant if ``execute=True``.

    Returns
    -------
    create_query
        A tuple in which the first element is a Peewee ``Table`` for the created table
        (the table is bound to ``database`` if passed), and the ``CREATE TABLE AS``
        query as a string.

    """

    if "." in table_name:
        schema, table_name = table_name.split(".")

    if schema is None and temporary is False:
        schema = "public"
    elif temporary is True:
        schema = None

    path = f"{schema}.{table_name}" if schema else table_name
    create_sql = f'CREATE {"TEMPORARY " if temporary else ""}TABLE {path} AS '

    if database is None and isinstance(query, peewee.ModelSelect):
        database = query.model._meta.database

    if overwrite and database:
        database.execute_sql(f"DROP TABLE IF EXISTS {path};")

    query_sql, params = database.get_sql_context().sql(query).query()
    cursor = database.cursor()
    query_str = cursor.mogrify(query_sql, params).decode()

    if execute:
        if database is None:
            raise RuntimeError("Cannot execute query without a database.")

        with database.atomic():
            database.execute_sql(create_sql + query_sql, params)

        for index in indices:
            if isinstance(index, (list, tuple)):
                index = ",".join(index)
            database.execute_sql(f"CREATE INDEX ON {path} ({index})")

        if analyze:
            database.execute_sql(f"VACUUM ANALYZE {path}")

    table = peewee.Table(table_name, schema=schema).bind(database)

    create_str = create_sql + query_str
    return table, create_str
