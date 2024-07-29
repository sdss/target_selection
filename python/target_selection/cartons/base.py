#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2020-02-24
# @Filename: base.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import abc
import datetime
import inspect
import re
import warnings

import numpy
import peewee
from astropy import table

from sdssdb.peewee import BaseModel
from sdssdb.peewee.sdss5db import catalogdb as cdb
from sdssdb.peewee.sdss5db import targetdb as tdb
from sdsstools import read_yaml_file
from sdsstools._vendor.color_print import color_text

from target_selection import __version__, config, log
from target_selection.exceptions import TargetSelectionError
from target_selection.utils import Timer


EPOCH = 2016.0


class BaseCarton(metaclass=abc.ABCMeta):
    """A base class for target cartons.

    This class is not intended for direct instantiation. Instead, it must be
    subclassed and the relevant class attributes overridden with the values
    corresponding to the carton.

    Parameters
    ----------
    targeting_plan : str
        The target selection plan version.
    config_file : str
        The path to the configuration file to use. If undefined, uses the
        internal ``target_selection.yml`` file.
    schema : str
        Schema in which the temporary table with the results of the
        query will be created. If `None`, tries to use the ``schema`` parameter
        from the configuration file for this plan of target selection. If
        the parameter is not set, defaults to ``'sandbox'``.
    table_name : str
        The name of the temporary table. Defaults to ``temp_<name>`` where
        ``<name>`` is the target class name.

    Attributes
    ----------
    name : str
        The name of the carton (required).
    cadence : str
        The label of the cadence rule for this carton.
    category : str
        The category of targets for this carton.
    mapper : str
        The mapper with which this carton is associated.
    orm : str
        The ORM library to be used, ``peewee`` or ``sqlalchemy``.
    tag : str
        The version of the ``target_selection`` code used.
    query_region : tuple
        A tuple defining the region over which the query should be performed,
        with the format ``(ra, dec, radius)`` in degrees. This will append a
        ``q3c_radial_query`` condition to the query.
    load_magnitudes : bool
        Whether to load target magnitudes. In general this must be `True`
        except for cartons for which it's known the magnitudes will not be
        used, e.g., skies.

    """

    name = None
    cadence = None
    category = None
    program = None
    mapper = None
    priority = None
    value = None
    instrument = None
    can_offset = False

    load_magnitudes = True

    query_region = None

    def __init__(self, targeting_plan, config_file=None, schema=None, table_name=None):
        assert self.name, "carton subclass must override name"
        assert self.category, "carton subclass must override category"
        assert self.program, "carton subclass must override program"

        self.plan = targeting_plan
        self.tag = __version__ if self.plan != "0.0.0-test" else None

        if config_file:
            this_config = read_yaml_file(config_file)
        else:
            this_config = config.copy()

        if self.plan not in this_config:
            raise TargetSelectionError(f"({self.name}): cannot find plan {self.plan!r} in config.")

        self.config = this_config[self.plan]

        if "parameters" in self.config:
            self.parameters = self.config["parameters"].get(self.name, None)
        else:
            self.parameters = None

        try:
            self.xmatch_plan = self.config["xmatch_plan"]
        except KeyError:
            raise TargetSelectionError(f"({self.name}): xmatch_plan not found in config.")

        # Check the signature of build_query
        self._build_query_signature = inspect.signature(self.build_query)
        if "version_id" not in self._build_query_signature.parameters:
            raise TargetSelectionError("build_query does not accept version_id")

        self.database = tdb.database
        assert self.database.connected, "database is not connected."

        self.schema = schema or self.config.get("schema", None) or "sandbox"
        self.table_name = table_name or f"temp_{self.name}".replace("-", "_")
        self.RModel = None

        if self.cadence:
            ncad = tdb.Cadence.select().where(tdb.Cadence.label == self.cadence).count()
            assert ncad == 1, f"{self.cadence!r} does not exist in targetdb.cadence."

        self.log = log
        self.has_run = False
        self._disable_query_log = False

        # We cannot set temp_buffers multiple times if there are temporary tables
        # (as in add_optical_magnitudes) so we set it here.
        if "database_options" in self.config and "temp_buffers" in self.config["database_options"]:
            temp_buffers = self.config["database_options"].pop("temp_buffers")
            self.database.execute_sql(f"SET temp_buffers = '{temp_buffers}'")

    @property
    def path(self):
        """The schema-qualified path to the output table."""

        return self.schema + "." + self.table_name

    @abc.abstractmethod
    def build_query(self, version_id, query_region=None):
        """Builds and returns the query.

        The ORM query for the target class. Note that this must be the
        *un-executed* query, which will be executed in `.run`. The select
        statement must be run on ``catalogdb.catalog`` and needs to include,
        at least, ``catalogid``. Additional columns such as
        ``ra, dec, pmra, pmdec`` can be specified but are otherwise completed
        when calling `.run`. Magnitude columns can be included and propagated
        to ``targetdb`` but must be aliased as ``magnitude_<band>`` where
        ``<band>`` must be one of the columns in ``targetdb.magnitude``.

        The query returned must include a filter for ``version_id`` on
        ``catalogdb.catalog``. Normally this is applied to the ``Catalog``
        select as ``.where(Catalog.version_id == version_id)``.

        Parameters
        ----------
        version_id : int
            The id of the cross-match version to use.

        Returns
        -------
        query
            A :class:`peewee:Select` or :class:`peewee:ModelSelect` query.

        """

        pass

    def get_model(self):
        """Returns a Peewee model for the temporary table using reflection."""

        # peewee has a Model class, BaseModel class, and ModelBase class.
        # The below Model class  (after this comment)
        # is different from peewee Model class.
        # The below BaseModel class (after this comment) is from the
        # sdssdb.peewee package's __init__.py.
        # It is different from the peewee BaseModel class.
        #
        # The below line is from sdssdb/peewee/__init__.py
        # This Model class in the below line is the peewee Model class.
        #     class BaseModel(Model, metaclass=ReflectMeta):
        #
        # The above line, relates the below Model class to
        # the ReflectMeta metaclass via the BaseModel class.
        # Inside the ReflectMeta metaclass is the code in __new__()
        # which generates the Model class.
        #
        # So defining the Model class here leads to creation of
        # the target_selection Model class. Later below we have the statement
        # "return Model" which is the last statement of get_model().
        #
        class Model(BaseModel):
            catalogid = peewee.BigIntegerField(primary_key=True)
            selected = peewee.BooleanField()
            cadence = peewee.TextField()
            instrument = peewee.TextField()
            priority = peewee.IntegerField()
            selected = peewee.BooleanField()

            class Meta:
                database = self.database
                table_name = self.table_name
                schema = self.schema
                reflection_options = {
                    "skip_foreign_keys": True,
                    "use_peewee_reflection": False,
                    "force": True,
                }
                use_reflection = True

        if not Model.table_exists():
            raise TargetSelectionError(f"temporary table {self.path!r} does not exist.")

        return Model

    def run(
        self,
        query_region=None,
        overwrite=False,
        limit=None,
        add_optical_magnitudes=True,
        **post_process_kawrgs,
    ):
        """Executes the query and post-process steps, and stores the results.

        This method calls `.build_query` and runs the returned query. The
        output of the query is stored in a temporary table whose schema and
        table name are defined when the object is instantiated.

        After the query has run, the `.post_process` routine is called if the
        method has been overridden for the given carton.

        Parameters
        ----------
        query_region : tuple
            A tuple defining the region over which the query should be
            performed, with the format ``(ra, dec, radius)`` in degrees. This
            will append a ``q3c_radial_query`` condition to the query.
        overwrite : bool
            Whether to overwrite the intermediary table if already exists.
        limit : int or None
            Limit the query to this number of targets. Useful for testing.
            The ``LIMIT`` statement is added to the query returned by
            `.build_query` and the exact behaviour will depend on the query.
        post_process_args : dict
            Keyword arguments to be passed to `.post_process`.

        Returns
        -------
        model : :class:`peewee:Model`
            The model for the intermediate table.

        """

        query_region = query_region or self.query_region

        path = self.path
        execute_sql = self.database.execute_sql

        if self.database.table_exists(self.table_name, schema=self.schema):
            if overwrite:
                self.log.info(f"Dropping table {path!r}.")
                self.drop_table()
            else:
                raise RuntimeError(f"Temporary table {path!r} already exists.")

        self.log.info("Running query ...")
        version_id = self.get_version_id()

        with Timer() as timer:
            # If build_query accepts a query_region parameter, call with the query region.
            # Otherwise will add the radial query condition later.
            if "query_region" in self._build_query_signature.parameters:
                query = self.build_query(version_id, query_region=query_region)
            else:
                query = self.build_query(version_id)

            # Make sure the catalogid column is selected.
            if cdb.Catalog.catalogid not in query._returning:
                raise RuntimeError("catalogid is not being returned in query.")

            if query_region:
                if "query_region" in self._build_query_signature.parameters:
                    pass
                else:
                    # This may be quite inefficient depending on the query.
                    subq = query.alias("subq")
                    query = (
                        peewee.Select(columns=[peewee.SQL("subq.*")])
                        .from_(subq)
                        .join(cdb.Catalog, on=(cdb.Catalog.catalogid == subq.c.catalogid))
                        .where(
                            peewee.fn.q3c_radial_query(
                                cdb.Catalog.ra,
                                cdb.Catalog.dec,
                                query_region[0],
                                query_region[1],
                                query_region[2],
                            )
                        )
                    )

            if limit:
                query = query.limit(limit)

            query_sql, params = self.database.get_sql_context().sql(query).query()
            cursor = self.database.cursor()
            query_str = cursor.mogrify(query_sql, params).decode()

            if not self._disable_query_log:
                log_message = f"CREATE TABLE IF NOT EXISTS {path} AS " + query_str
                if self.log.rich_console:
                    self.log.debug(log_message, extra={"highlighter": None})
                else:
                    self.log.debug(color_text(log_message, "darkgrey"))
            else:
                self.log.debug("Not printing VERY long query.")

            with self.database.atomic():
                self.setup_transaction()
                execute_sql(f"CREATE TABLE IF NOT EXISTS {path} AS " + query_sql, params)

        self.log.info(f"Created table {path!r} in {timer.interval:.3f} s.")

        # Note that above execute_sql() is the same as below
        # self.database.execute_sql() since at the top of this method
        # we have the below statement.
        #        execute_sql = self.database.execute_sql

        self.database.execute_sql("COMMIT;")
        self.RModel = self.get_model()

        self.log.debug("Adding columns and indexes.")

        columns = [col.name for col in self.database.get_columns(self.table_name, self.schema)]

        if "selected" not in columns:
            execute_sql(f"ALTER TABLE {path} ADD COLUMN selected BOOL DEFAULT TRUE;")
            self.RModel._meta.add_field("selected", peewee.BooleanField())

        if "cadence" not in columns and self.cadence is None:
            execute_sql(f"ALTER TABLE {path} ADD COLUMN cadence VARCHAR;")
            self.RModel._meta.add_field("cadence", peewee.TextField())

        if "priority" not in columns and self.priority is None:
            execute_sql(f"ALTER TABLE {path} ADD COLUMN priority INTEGER;")
            self.RModel._meta.add_field("priority", peewee.IntegerField())

        if "instrument" not in columns and self.instrument is None:
            execute_sql(f"ALTER TABLE {path} ADD COLUMN instrument TEXT;")
            self.RModel._meta.add_field("instrument", peewee.TextField())

        execute_sql(f"ALTER TABLE {path} ADD PRIMARY KEY (catalogid);")
        execute_sql(f"CREATE INDEX ON {path} (selected);")
        execute_sql(f"ANALYZE {path};")

        n_rows = self.RModel.select().count()
        self.log.debug(f"Table {path!r} contains {n_rows:,} rows.")

        self.log.debug("Running post-process.")
        with self.database.atomic():
            self.setup_transaction()
            self.post_process(self.RModel, **post_process_kawrgs)

        n_selected = self.RModel.select().where(self.RModel.selected >> True).count()
        self.log.debug(f"Selected {n_selected:,} rows after post-processing.")

        if add_optical_magnitudes:
            self.log.debug("Adding optical magnitude columns.")
            self.add_optical_magnitudes()

        self.has_run = True

        return self.RModel

    def get_version_id(self):
        """Returns the version_id for the cross-match plan."""

        return cdb.Version.get(plan=self.xmatch_plan).id

    def add_optical_magnitudes(self):
        """Adds ``gri`` magnitude columns."""

        Model = self.RModel

        magnitudes = ["g", "r", "i", "z"]

        # Check if ALL the columns have already been created in the query.
        # If so, just return.
        if any([mag in Model._meta.columns for mag in magnitudes]):
            if not all([mag in Model._meta.columns for mag in magnitudes]):
                raise TargetSelectionError(
                    "Some optical magnitudes are defined in the query but not all of them."
                )
            if "optical_prov" not in Model._meta.columns:
                raise TargetSelectionError("optical_prov column does not exist.")
            self.log.warning("All optical magnitude columns are defined in the query.")
            return

        # First create the columns. Also create z to speed things up. We won't
        # use transformations for z but we can use the initial query to populate
        # it and avoid doing the same query later when loading the magnitudes.
        for mag in magnitudes:
            self.database.execute_sql(f"ALTER TABLE {self.path} ADD COLUMN {mag} REAL;")
            Model._meta.add_field(mag, peewee.FloatField())

        self.database.execute_sql(f"ALTER TABLE {self.path} ADD COLUMN optical_prov TEXT;")
        Model._meta.add_field("optical_prov", peewee.TextField())

        # Step 1: join with gaia_dr3_synthetic_photometry_gspc and use synthetic SDSS mags.

        SynPhot = cdb.Gaia_dr3_synthetic_photometry_gspc
        with self.database.atomic():
            self.database.execute_sql("DROP TABLE IF EXISTS " + self.table_name + "_g3xp")
            temp_table = peewee.Table(self.table_name + "_g3xp")

            (
                Model.select(
                    Model.catalogid,
                    SynPhot.g_sdss_mag,
                    SynPhot.r_sdss_mag,
                    SynPhot.i_sdss_mag,
                    SynPhot.z_sdss_mag,
                )
                .join(
                    cdb.CatalogToGaia_DR3,
                    on=(cdb.CatalogToGaia_DR3.catalogid == Model.catalogid),
                )
                .join(SynPhot, on=(cdb.CatalogToGaia_DR3.target_id == SynPhot.source_id))
                .where(Model.g.is_null() | Model.r.is_null() | Model.i.is_null())
                .where(Model.selected >> True)
                .where(
                    cdb.CatalogToGaia_DR3.best >> True,
                    cdb.CatalogToGaia_DR3.version_id == self.get_version_id(),
                )
                .where(
                    SynPhot.g_sdss_mag.is_null(False),
                    SynPhot.r_sdss_mag.is_null(False),
                    SynPhot.i_sdss_mag.is_null(False),
                    # Now removing the flag selection because it killed everything
                    # with r_sdss_mag < 14.2
                    # SynPhot.g_sdss_flag == 1,
                    # SynPhot.r_sdss_flag == 1,
                    # SynPhot.i_sdss_flag == 1,
                    # It's possible that in future, we might want to apply cuts on bp-rp
                    # or on the derived sdss colours
                    # Fairly arbitrary faint cut
                    SynPhot.r_sdss_mag < 15.0,
                    # Could potentially go fainter, or could cut on G,BP,RP etc instead
                    # but that would require a join to gaia_dr3_source table
                    #
                    # The following cut on SNR is completely unnecessary at G<15
                    # (SynPhot.r_sdss_flux_error <
                    #  SynPhot.r_sdss_flux / 30.0),
                )
                .create_table(temp_table._path[0], temporary=True)
            )

            nrows = (
                Model.update(
                    {
                        Model.g: temp_table.c.g_sdss_mag,
                        Model.r: temp_table.c.r_sdss_mag,
                        Model.i: temp_table.c.i_sdss_mag,
                        Model.optical_prov: peewee.Value("sdss_psfmag_from_g3xp"),
                    }
                )
                .from_(temp_table)
                .where(Model.catalogid == temp_table.c.catalogid)
                .where(
                    temp_table.c.g_sdss_mag.is_null(False)
                    & temp_table.c.r_sdss_mag.is_null(False)
                    & temp_table.c.i_sdss_mag.is_null(False)
                )
            ).execute()

        self.log.debug(f"{nrows:,} associated with Gaia DR3 XP magnitudes.")

        # Step 2: localise entries with empty magnitudes,
        # join with sdss_dr13_photoobj and use SDSS magnitudes.

        # https://www.sdss4.org/dr16/algorithms/photo_flags/
        # https://www.sdss4.org/dr16/algorithms/photo_flags_recommend/
        # we ignore SDSS photometry when any of these bits are set
        # TODO improve/expand the list of rejected bits
        sdss_quality_bitmask = (
            0  # just a placeholder
            + (1 << 7)  # NOPROFILE
            + (1 << 18)  # SATURATED
            + (1 << 22)  # BADSKY
            + (1 << 24)  # TOOLARGE
            + (1 << (16 + 32))  # TOO_FEW_GOOD_DETECTIONS
        )

        with self.database.atomic():
            self.database.execute_sql("DROP TABLE IF EXISTS " + self.table_name + "_sdss")
            temp_table = peewee.Table(self.table_name + "_sdss")

            SDSS_DR13 = cdb.SDSS_DR13_PhotoObj
            (
                Model.select(
                    Model.catalogid,
                    SDSS_DR13.psfmag_g,
                    SDSS_DR13.psfmag_r,
                    SDSS_DR13.psfmag_i,
                    SDSS_DR13.psfmag_z,
                )
                .join(
                    cdb.CatalogToSDSS_DR13_PhotoObj_Primary,
                    on=(cdb.CatalogToSDSS_DR13_PhotoObj_Primary.catalogid == Model.catalogid),
                )
                .join(
                    SDSS_DR13,
                    on=(cdb.CatalogToSDSS_DR13_PhotoObj_Primary.target_id == SDSS_DR13.objid),
                )
                .where(Model.g.is_null() | Model.r.is_null() | Model.i.is_null())
                .where(
                    cdb.CatalogToSDSS_DR13_PhotoObj_Primary.best >> True,
                    cdb.CatalogToSDSS_DR13_PhotoObj_Primary.version_id == self.get_version_id(),
                )
                .where(
                    SDSS_DR13.psfmag_r > 14.0,  # avoid stars close to saturation
                    SDSS_DR13.psfmag_g > -99.0,  # avoid bad mags
                    SDSS_DR13.psfmag_i > -99.0,  # avoid bad mags
                    SDSS_DR13.flags.bin_and(sdss_quality_bitmask) == 0,
                )
                .where(Model.selected >> True)
                .create_table(temp_table._path[0], temporary=True)
            )

            nrows = (
                Model.update(
                    {
                        Model.g: temp_table.c.psfmag_g,
                        Model.r: temp_table.c.psfmag_r,
                        Model.i: temp_table.c.psfmag_i,
                        Model.optical_prov: peewee.Value("sdss_psfmag"),
                    }
                )
                .from_(temp_table)
                .where(Model.catalogid == temp_table.c.catalogid)
                .where(
                    temp_table.c.psfmag_g.is_null(False)
                    & temp_table.c.psfmag_r.is_null(False)
                    & temp_table.c.psfmag_i.is_null(False)
                )
            ).execute()

        self.log.debug(f"{nrows:,} associated with SDSS magnitudes.")

        # Step 3: localise entries with empty magnitudes and use PanSTARRS1
        # transformations.

        # PS1 fluxes are in Janskys. We use stacked fluxes instead of mean
        # magnitudes since they are more complete on the faint end.
        ps1_g = 8.9 - 2.5 * peewee.fn.log(cdb.Panstarrs1.g_stk_psf_flux)
        ps1_r = 8.9 - 2.5 * peewee.fn.log(cdb.Panstarrs1.r_stk_psf_flux)
        ps1_i = 8.9 - 2.5 * peewee.fn.log(cdb.Panstarrs1.i_stk_psf_flux)

        # Use transformations to SDSS from Tonry et al. 2012 (section 3.2, table 6).
        x = ps1_g - ps1_r
        ps1_sdss_g = 0.013 + 0.145 * x + 0.019 * x * x + ps1_g
        ps1_sdss_r = -0.001 + 0.004 * x + 0.007 * x * x + ps1_r
        ps1_sdss_i = -0.005 + 0.011 * x + 0.010 * x * x + ps1_i

        # limit panstarrs photometry to r > 14
        max_flux_r = 10 ** ((14.0 - 8.9) / -2.5)

        with self.database.atomic():
            self.database.execute_sql("DROP TABLE IF EXISTS " + self.table_name + "_ps1")
            temp_table = peewee.Table(self.table_name + "_ps1")

            (
                Model.select(
                    Model.catalogid,
                    ps1_sdss_g.alias("ps1_sdss_g"),
                    ps1_sdss_r.alias("ps1_sdss_r"),
                    ps1_sdss_i.alias("ps1_sdss_i"),
                )
                .join(
                    cdb.CatalogToPanstarrs1,
                    on=(cdb.CatalogToPanstarrs1.catalogid == Model.catalogid),
                )
                .join(
                    cdb.Panstarrs1,
                    on=(cdb.CatalogToPanstarrs1.target_id == cdb.Panstarrs1.catid_objid),
                )
                .where(Model.g.is_null() | Model.r.is_null() | Model.i.is_null())
                .where(Model.selected >> True)
                .where(
                    cdb.CatalogToPanstarrs1.best >> True,
                    cdb.CatalogToPanstarrs1.version_id == self.get_version_id(),
                )
                .where(cdb.Panstarrs1.r_stk_psf_flux < max_flux_r)
                .where(
                    cdb.Panstarrs1.g_stk_psf_flux.is_null(False),
                    cdb.Panstarrs1.g_stk_psf_flux != "NaN",
                    (cdb.Panstarrs1.g_stk_psf_flux > 0),
                )
                .where(
                    cdb.Panstarrs1.r_stk_psf_flux.is_null(False),
                    cdb.Panstarrs1.r_stk_psf_flux != "NaN",
                    (cdb.Panstarrs1.r_stk_psf_flux > 0),
                )
                .where(
                    cdb.Panstarrs1.i_stk_psf_flux.is_null(False),
                    cdb.Panstarrs1.i_stk_psf_flux != "NaN",
                    (cdb.Panstarrs1.i_stk_psf_flux > 0),
                )
                .create_table(temp_table._path[0], temporary=True)
            )

            nrows = (
                Model.update(
                    {
                        Model.g: temp_table.c.ps1_sdss_g,
                        Model.r: temp_table.c.ps1_sdss_r,
                        Model.i: temp_table.c.ps1_sdss_i,
                        Model.optical_prov: peewee.Value("sdss_psfmag_ps1"),
                    }
                )
                .from_(temp_table)
                .where(Model.catalogid == temp_table.c.catalogid)
                .where(
                    temp_table.c.ps1_sdss_g.is_null(False)
                    & temp_table.c.ps1_sdss_r.is_null(False)
                    & temp_table.c.ps1_sdss_i.is_null(False)
                )
            ).execute()

        self.log.debug(f"{nrows:,} associated with PS1 magnitudes.")

        # Step 4: localise entries with empty magnitudes and use Gaia transformations
        # from Evans et al (2018).

        gaia_G = cdb.Gaia_DR3.phot_g_mean_mag
        gaia_BP = cdb.Gaia_DR3.phot_bp_mean_mag
        gaia_RP = cdb.Gaia_DR3.phot_rp_mean_mag

        x = gaia_BP - gaia_RP
        x2 = x * x
        x3 = x * x * x
        gaia_sdss_g = -1 * (0.13518 - 0.46245 * x - 0.25171 * x2 + 0.021349 * x3) + gaia_G
        gaia_sdss_r = -1 * (-0.12879 + 0.24662 * x - 0.027464 * x2 - 0.049465 * x3) + gaia_G
        gaia_sdss_i = -1 * (-0.29676 + 0.64728 * x - 0.10141 * x2) + gaia_G

        with self.database.atomic():
            self.database.execute_sql("DROP TABLE IF EXISTS " + self.table_name + "_gaia")
            temp_table = peewee.Table(self.table_name + "_gaia")

            (
                Model.select(
                    Model.catalogid,
                    gaia_sdss_g.alias("gaia_sdss_g"),
                    gaia_sdss_r.alias("gaia_sdss_r"),
                    gaia_sdss_i.alias("gaia_sdss_i"),
                )
                .join(
                    cdb.CatalogToGaia_DR3,
                    on=(cdb.CatalogToGaia_DR3.catalogid == Model.catalogid),
                )
                .join(cdb.Gaia_DR3)
                .where(Model.g.is_null() | Model.r.is_null() | Model.i.is_null())
                .where(Model.selected >> True)
                .where(
                    cdb.CatalogToGaia_DR3.best >> True,
                    cdb.CatalogToGaia_DR3.version_id == self.get_version_id(),
                )
                .where(cdb.Gaia_DR3.phot_g_mean_mag.is_null(False))
                .where(cdb.Gaia_DR3.phot_bp_mean_mag.is_null(False))
                .where(cdb.Gaia_DR3.phot_rp_mean_mag.is_null(False))
                .create_table(temp_table._path[0], temporary=True)
            )

            nrows = (
                Model.update(
                    {
                        Model.g: temp_table.c.gaia_sdss_g,
                        Model.r: temp_table.c.gaia_sdss_r,
                        Model.i: temp_table.c.gaia_sdss_i,
                        Model.optical_prov: peewee.Value("sdss_psfmag_gaia"),
                    }
                )
                .from_(temp_table)
                .where(Model.catalogid == temp_table.c.catalogid)
                .where(
                    temp_table.c.gaia_sdss_g.is_null(False)
                    & temp_table.c.gaia_sdss_r.is_null(False)
                    & temp_table.c.gaia_sdss_i.is_null(False)
                )
            ).execute()

        self.log.debug(f"{nrows:,} associated with Gaia DR3 mean magnitudes.")

        # Finally, check if there are any rows in which at least some of the
        # magnitudes are null.

        n_empty = (
            Model.select()
            .where(Model.g.is_null() | Model.r.is_null() | Model.i.is_null())
            .where(Model.selected >> True)
            .count()
        )

        if n_empty > 0:
            log.warning(f"Found {n_empty} entries with empty magnitudes.")

    def post_process(self, model, **kwargs):
        """Post-processes the temporary table.

        This method provides a framework for applying non-SQL operations on
        carton query. It receives the model for the temporary table and can
        perform any operation on it, including modifying the ``selected``
        column with a mask of targets to be used.

        This method can also be used to set the ``cadence`` column in the
        temporary table. This column will be used to set the target cadence if
        the carton `.cadence` attribute is not set.

        `.post_process` runs inside a database transaction so it's not
        necessary to create a new one, but savepoints can be added.

        Parameters
        ----------
        model : peewee:Model
            The model of the intermediate table.

        Returns
        -------
        mask : `tuple`
            The list of catalogids from the temporary table that should be
            selected as part of this carton. If `True` (the default), selects
            all the records.

        """

        return True

    def setup_transaction(self):
        """Setups the transaction locally modifying the datbase parameters.

        This method runs inside a transaction and can be overridden to set the
        parameters of the transaction manually. It applies to both `.run` and
        `.load`.

        """

        if "database_options" not in self.config:
            return

        for param, value in self.config["database_options"].items():
            self.database.execute_sql(f"SET LOCAL {param} = {value!r};")

    def drop_table(self):
        """Drops the intermediate table if it exists."""

        self.database.execute_sql(f"DROP TABLE IF EXISTS {self.path};")

    def write_table(self, filename=None, mode="results", write=True):
        """Writes the selection to a FITS file.

        Parameters
        ----------
        filename : str
            The file to which to write the table. Defaults to
            ``<name>_<plan>.fits``.
        mode : str
            Defines what data to write. If ``'results'``, writes the
            intermediate table (usually just the ``catalogid`` column). If
            ``'targetdb'``, writes all the relevant columns for the targets
            loaded to ``targetdb`` for this carton and plan (must be used after
            `.load` has been called).
        write : bool
            Whether to write the table to disk. If `False`, just returns the
            table object.

        Returns
        -------
        table : `~astropy.table.Table`
            A table object with the selected results.

        """

        if not filename:
            if mode == "results":
                filename = f"{self.name}_{self.plan}.fits.gz"
            else:
                filename = f"{self.name}_{self.plan}_targetdb.fits.gz"

        self.log.debug(f"Writing table to {filename}.")

        if not self.RModel:
            self.database.execute_sql("COMMIT;")
            self.RModel = self.get_model()

        if mode == "results":
            results_model = self.RModel
            assert results_model.table_exists(), "temporary table does not exist."

            write_query = results_model.select().order_by(results_model.catalogid)

            colnames = [field.name for field in write_query._returning]

        elif mode == "targetdb":
            mag_fields = [
                field
                for field in tdb.Magnitude._meta.fields.values()
                if field.name not in ["pk", "target_pk", "target"]
            ]

            write_query = (
                tdb.Target.select(
                    tdb.Target,
                    *mag_fields,
                    tdb.CartonToTarget.priority,
                    tdb.CartonToTarget.value,
                    tdb.CartonToTarget.instrument_pk,
                    tdb.CartonToTarget.inertial,
                    tdb.CartonToTarget.can_offset,
                    tdb.Cadence.label.alias("cadence"),
                )
                .join(tdb.CartonToTarget)
                .join(tdb.Magnitude)
                .join_from(tdb.CartonToTarget, tdb.Cadence, peewee.JOIN.LEFT_OUTER)
                .join_from(tdb.CartonToTarget, tdb.Carton)
                .join(tdb.Version)
                .where(
                    tdb.Carton.carton == self.name,
                    tdb.Version.plan == self.plan,
                    tdb.Version.target_selection >> True,
                )
                .order_by(tdb.Target.catalogid)
            )

            colnames = []
            for col in write_query._returning:
                if isinstance(col, peewee.ForeignKeyField):
                    colnames.append(col.column_name)
                elif isinstance(col, peewee.Alias):
                    colnames.append(col._alias)
                else:
                    colnames.append(col.name)

        else:
            raise ValueError('invalud mode. Available modes are "results" and "targetdb".')

        if not write_query.exists():
            raise TargetSelectionError("no records found.")

        results = write_query.tuples()
        results = (
            (col if col is not None else numpy.nan for col in row) for row in tuple(results)
        )

        warnings.filterwarnings("ignore", message=".*converting a masked element to nan.*")

        carton_table = table.Table(rows=results, names=colnames, masked=True)

        if write:
            carton_table.write(filename, overwrite=True)

        return carton_table

    def load(self, mode="fail", overwrite=False):
        """Loads the output of the intermediate table into targetdb.

        Parameters
        ----------
        mode : str
            The mode to use when loading the targets. If ``'fail'``, raises an
            error if the carton already exist. If ``'overwrite'``, overwrites
            the targets. If ``'append'``, appends the targets.
        overwrite : bool
            Equivalent to setting ``mode='overwrite'``. This option is deprecated and
            will raise a warning.

        """

        if overwrite:
            mode = "overwrite"
            log.warning(
                "The `overwrite` option is deprecated and will be removed in a future version. "
                'Use `mode="overwrite"` instead.'
            )

        if self.check_targets():
            if mode == "overwrite":
                log.warning(
                    f"Carton {self.name!r} with plan {self.plan!r} "
                    f"already has targets loaded. "
                    "Dropping carton-to-target entries."
                )
                self.drop_carton()
            elif mode == "append":
                pass
            elif mode == "fail":
                raise TargetSelectionError(
                    f"Found existing targets for "
                    f"carton {self.name!r} with plan "
                    f"{self.plan!r}."
                )
            else:
                raise ValueError(f'Invalid mode {mode!r}. Use "fail", "overwrite", or "append".')

        if self.RModel is None:
            self.database.execute_sql("COMMIT;")
            RModel = self.get_model()
        else:
            RModel = self.RModel

        if not RModel.table_exists():
            raise TargetSelectionError(
                f"No temporary table found {self.full}. Did you call run()?"
            )

        has_targets = RModel.select().where(RModel.selected >> True).exists()

        if not has_targets:
            log.warning("No targets found in intermediate table.")

        with self.database.atomic():
            self.setup_transaction()
            self._create_carton_metadata()
            self._load_targets(RModel)
            self._load_carton_to_target(RModel)
            if self.load_magnitudes:
                self._load_magnitudes(RModel)
            else:
                log.warning("Skipping magnitude load.")

            self.log.debug("Committing records and checking constraints.")

    def check_targets(self):
        """Check if data has been loaded for this carton and targeting plan."""

        has_targets = (
            tdb.CartonToTarget.select()
            .join(tdb.Carton)
            .join(tdb.Version)
            .where(
                tdb.Carton.carton == self.name,
                tdb.Version.plan == self.plan,
                tdb.Version.target_selection >> True,
            )
            .exists()
        )

        return has_targets

    def _create_carton_metadata(self):
        mapper_pk = None
        category_pk = None

        # Create targeting plan in tdb.
        version, created = tdb.Version.get_or_create(
            plan=self.plan,
            tag=self.tag,
            target_selection=True,
            robostrategy=False,
        )
        version_pk = version.pk

        if created:
            self.log.info(
                f"Created record in targetdb.version for {self.plan!r} with tag {self.tag!r}."
            )

        if (
            tdb.Carton.select()
            .where(tdb.Carton.carton == self.name, tdb.Carton.version_pk == version_pk)
            .exists()
        ):
            return

        # Create carton and associated values.
        if self.mapper:
            mapper, created_pk = tdb.Mapper.get_or_create(label=self.mapper)
            mapper_pk = mapper.pk
            if created:
                self.log.debug(f"Created mapper {self.mapper!r}")

        if self.category:
            category, created = tdb.Category.get_or_create(label=self.category)
            category_pk = category.pk
            if created:
                self.log.debug(f"Created category {self.category!r}")

        tdb.Carton.create(
            carton=self.name,
            category_pk=category_pk,
            program=self.program,
            mapper_pk=mapper_pk,
            version_pk=version_pk,
            run_on=datetime.datetime.now().isoformat().split("T")[0],
        ).save()

        self.log.debug(f"Created carton {self.name!r}")

    def _load_targets(self, RModel):
        """Load data from the intermediate table tp targetdb.target."""

        self.log.debug("loading data into targetdb.target.")

        (
            tdb.Target.insert_from(
                cdb.Catalog.select(
                    cdb.Catalog.catalogid,
                    cdb.Catalog.ra,
                    cdb.Catalog.dec,
                    cdb.Catalog.pmra,
                    cdb.Catalog.pmdec,
                    cdb.Catalog.parallax,
                    peewee.Value(EPOCH),
                )
                .join(RModel, on=(cdb.Catalog.catalogid == RModel.catalogid))
                .where(RModel.selected >> True)
                .where(
                    ~peewee.fn.EXISTS(
                        tdb.Target.select(peewee.SQL("1")).where(
                            tdb.Target.catalogid == RModel.catalogid
                        )
                    )
                ),
                [
                    tdb.Target.catalogid,
                    tdb.Target.ra,
                    tdb.Target.dec,
                    tdb.Target.pmra,
                    tdb.Target.pmdec,
                    tdb.Target.parallax,
                    tdb.Target.epoch,
                ],
            )
            .returning()
            .execute()
        )

        self.log.info("Inserted new rows into targetdb.target.")

        return

    def _load_magnitudes(self, RModel):
        """Load magnitudes into targetdb.magnitude."""

        self.log.debug("Loading data into targetdb.magnitude.")

        Magnitude = tdb.Magnitude

        magnitude_paths = self.config["magnitudes"]
        fields = [Magnitude.carton_to_target_pk]

        select_from = (
            RModel.select(tdb.CartonToTarget.pk)
            .join(tdb.Target, on=(RModel.catalogid == tdb.Target.catalogid))
            .join(tdb.CartonToTarget)
            .join(tdb.Carton)
            .join(tdb.Version)
            .where(RModel.selected >> True)
            .where(
                tdb.Carton.carton == self.name,
                tdb.Version.plan == self.plan,
                tdb.Version.tag == self.tag,
                tdb.Version.target_selection >> True,
            )
        )

        for mag, mpath in magnitude_paths.items():
            fields.append(getattr(Magnitude, mag))
            if hasattr(RModel, mag):
                select_from = select_from.select_extend(getattr(RModel, mag))
                continue

            select_from = select_from.switch(tdb.Target)

            # For each node in the join list we check if the node model has
            # already been join and if so, switch the pointer to that model.
            # Otherwise do a LEFT OUTER join because we want all the rows in
            # the temporary table even if they don't have associated
            # magnitudes. We make sure we only use the "best" match from the
            # cross-match. No need to apply filters on cross-match version_id
            # because catalogid is unique across versions.

            for node in mpath:
                column = None
                if node == mpath[-1]:
                    node, column = node.split(".")
                node_model = self.database.models.get("catalogdb." + node)
                joins = [model[0] for join in select_from._joins.values() for model in join]
                if node_model in joins:
                    select_from = select_from.switch(node_model)
                else:
                    if node.startswith("catalog_to_"):
                        select_from = select_from.join(
                            node_model,
                            peewee.JOIN.LEFT_OUTER,
                            on=(tdb.Target.catalogid == node_model.catalogid),
                        )
                        select_from = select_from.where(
                            (
                                (node_model.best >> True)
                                & (node_model.version_id == self.get_version_id())
                            )
                            | (node_model.catalogid >> None)
                        )
                    else:
                        select_from = select_from.join(node_model, peewee.JOIN.LEFT_OUTER)
                if column:
                    select_from = select_from.select_extend(getattr(node_model, column))

        # Add gri from the temporary table.
        for mag in ["g", "r", "i", "z"]:
            select_from = select_from.select_extend(RModel._meta.columns[mag])
            fields.append(Magnitude._meta.columns[mag])

        if "optical_prov" not in RModel._meta.columns:
            raise TargetSelectionError("optical_prov column not found in temporary table.")

        select_from = select_from.select_extend(RModel._meta.columns["optical_prov"])
        fields.append(Magnitude.optical_prov)

        Magnitude.insert_from(select_from, fields).returning().execute()

        self.log.info("Inserted new rows into targetdb.magnitude.")

    def _load_carton_to_target(self, RModel):
        """Populate targetdb.carton_to_target."""

        self.log.debug("Loading data into targetdb.carton_to_target.")

        version_pk = tdb.Version.get(
            plan=self.plan,
            tag=self.tag,
            target_selection=True,
        )
        carton_pk = tdb.Carton.get(carton=self.name, version_pk=version_pk).pk

        Target = tdb.Target
        CartonToTarget = tdb.CartonToTarget

        select_from = (
            RModel.select(Target.pk, carton_pk)
            .join(Target, on=(Target.catalogid == RModel.catalogid))
            .where(RModel.selected >> True)
            .where(
                ~peewee.fn.EXISTS(
                    CartonToTarget.select(peewee.SQL("1"))
                    .join(tdb.Carton)
                    .where(
                        CartonToTarget.target_pk == Target.pk,
                        CartonToTarget.carton_pk == carton_pk,
                        tdb.Carton.version_pk == version_pk,
                    )
                )
            )
        )

        if self.cadence is not None:
            # Check that not both the carton cadence and the cadence column
            # are not null.
            if "cadence" in RModel._meta.fields:
                if RModel.select().where(~(RModel.cadence >> None)).exists():
                    raise TargetSelectionError(
                        "both carton cadence and target "
                        "cadence defined. This is not "
                        "allowed."
                    )

            cadence_pk = tdb.Cadence.get(label=self.cadence)
            select_from = select_from.select_extend(cadence_pk)

            if not self.value:
                # improve robustness of cadence name pattern matching slightly:
                try:
                    chunks = [s for s in self.cadence.split("_") if re.match(r"[0-9]+x[0-9]+", s)]
                    cadence_payload = chunks[0].split("x")
                except BaseException:
                    raise ("Uninterpretable cadence name: ", self.cadence)

                self.value = float(
                    numpy.multiply(
                        # *map(int, self.cadence.split('_')[-1].split('x'))
                        *map(int, cadence_payload)
                    )
                )

        else:
            # If all cadences are null we'll set that as a value and save us
            # a costly join.
            if not RModel.select().where(~(RModel.cadence >> None)).exists():
                select_from = select_from.select_extend(peewee.SQL("null"))
            else:
                select_from = (
                    select_from.select_extend(tdb.Cadence.pk)
                    .switch(RModel)
                    .join(
                        tdb.Cadence,
                        "LEFT OUTER JOIN",
                        on=(tdb.Cadence.label == RModel.cadence),
                    )
                )

        if self.priority is None:
            select_from = select_from.select_extend(RModel.priority)
        else:
            select_from = select_from.select_extend(self.priority)

        if self.can_offset is None:
            select_from = select_from.select_extend(RModel.can_offset)
        else:
            select_from = select_from.select_extend(self.can_offset)

        if self.value is not None:
            select_from = select_from.select_extend(self.value)
        else:
            # We will use the cadence to determine the value. First, if there is
            # not a user-defined value column, create it.
            if "value" not in RModel._meta.columns:
                self.database.execute_sql(f"ALTER TABLE {self.path} ADD COLUMN value REAL;")

                # We need to add the field like this and not call get_model() because
                # at this point the temporary table is locked and reflection won't work.
                RModel._meta.add_field("value", peewee.FloatField())

            # Get the value as the n_epochs * n_exposures_per_epoch. Probably this can
            # be done directly in SQL but it's just easier in Python. Note that because
            # we set value above in the case when cadence is a single value, if we
            # are here that means there is a cadence column.

            data = numpy.array(
                RModel.select(RModel.catalogid, RModel.cadence)
                .where(RModel.cadence.is_null(False))
                .tuples()
            )

            def split_cadence(cadence):
                chunks = [s for s in cadence.split("_") if re.match(r"[0-9]+x[0-9]+", s)]
                return chunks[0].split("x")

            if data.size > 0:
                values = tuple(
                    int(
                        numpy.multiply(
                            # *map(int, cadence.split('_')[-1].split('x'))))
                            # improve robustness of cadence name pattern matching slightly:
                            *map(int, split_cadence(cadence))
                        )
                    )
                    for cadence in data[:, 1]
                )

                catalogid_values = zip(map(int, data[:, 0]), values)

                vl = peewee.ValuesList(
                    catalogid_values,
                    columns=("catalogid", "value"),
                    alias="vl",
                )

                (
                    RModel.update(value=vl.c.value)
                    .from_(vl)
                    .where(RModel.catalogid == vl.c.catalogid)
                    .where(RModel.value.is_null())
                ).execute()

            select_from = select_from.select_extend(RModel.value)

        if "instrument" in RModel._meta.columns:
            select_from = (
                select_from.select_extend(tdb.Instrument.pk)
                .switch(RModel)
                .join(
                    tdb.Instrument,
                    "LEFT OUTER JOIN",
                    on=(tdb.Instrument.label == RModel.instrument),
                )
            )
        elif self.instrument is not None:
            select_from = select_from.select_extend(tdb.Instrument.get(label=self.instrument).pk)
        else:
            raise RuntimeError(f"Instrument not defined for carton {self.name}")

        for colname in ["delta_ra", "delta_dec", "inertial"]:
            if colname in RModel._meta.columns:
                select_from = select_from.select_extend(RModel._meta.columns[colname])
            else:
                if colname == "inertial":
                    select_from = select_from.select_extend(peewee.Value(False))
                else:
                    select_from = select_from.select_extend(peewee.Value(0.0))

        if "lambda_eff" in RModel._meta.columns:
            select_from = select_from.select_extend(RModel._meta.columns["lambda_eff"])
        else:
            if self.instrument is not None:
                instrument = self.instrument
            else:
                instrument = RModel.instrument
            select_from = select_from.select_extend(
                tdb.Instrument.select(tdb.Instrument.default_lambda_eff).where(
                    tdb.Instrument.label == instrument
                )
            )

        # Now do the insert
        (
            CartonToTarget.insert_from(
                select_from,
                [
                    CartonToTarget.target_pk,
                    CartonToTarget.carton_pk,
                    CartonToTarget.cadence_pk,
                    CartonToTarget.priority,
                    CartonToTarget.can_offset,
                    CartonToTarget.value,
                    CartonToTarget.instrument_pk,
                    CartonToTarget.delta_ra,
                    CartonToTarget.delta_dec,
                    CartonToTarget.inertial,
                    CartonToTarget.lambda_eff,
                ],
            )
            .returning()
            .execute()
        )

        self.log.info("Inserted rows into targetdb.carton_to_target.")

    def drop_carton(self):
        """Drops the entry in ``targetdb.carton``."""

        version = tdb.Version.select().where(
            tdb.Version.plan == self.plan,
            tdb.Version.tag == self.tag,
            tdb.Version.target_selection >> True,
        )

        if version.count() == 0:
            return

        tdb.Carton.delete().where(
            tdb.Carton.carton == self.name, tdb.Carton.version == version
        ).execute()
