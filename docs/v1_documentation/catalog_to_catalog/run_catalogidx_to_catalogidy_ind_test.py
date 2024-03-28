from create_catalogidx_to_catalogidy import MetaXMatch, create_unique_from_region
from sdssdb.peewee.sdss5db.catalogdb import database

metax = MetaXMatch('catalogidx_to_catalogidy_ind_test.yml', database)

metax.run()
create_unique_from_region("catalogidx_to_catalogidy_catalog_to_allstar_dr17_synspec_rev1")
