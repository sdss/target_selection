from create_catalogidx_to_catalogidy import MetaXMatch, create_unique_from_region

from sdssdb.peewee.sdss5db.catalogdb import database


metax = MetaXMatch("catalogidx_to_catalogidy.yml", database)

metax.run()
create_unique_from_region("catalogidx_to_catalogidy")
