-- This script creates sequences for the primary keys (sdss_id and pk) in the two sdss_id tables

CREATE SEQUENCE catalogdb.sdss_id_stacked_sdss_id__seq;
ALTER TABLE catalogdb.sdss_id_stacked ALTER COLUMN sdss_id SET DEFAULT nextval('catalogdb.sdss_id_stacked_sdss_id__seq');
ALTER SEQUENCE catalogdb.sdss_id_stacked_sdss_id__seq OWNED BY catalogdb.sdss_id_stacked.sdss_id;
SELECT setval('catalogdb.sdss_id_stacked_sdss_id__seq', (SELECT MAX(sdss_id) FROM catalogdb.sdss_id_stacked));
 
CREATE SEQUENCE catalogdb.sdss_id_flat_pk__seq;
ALTER TABLE catalogdb.sdss_id_flat ALTER COLUMN pk SET DEFAULT nextval('catalogdb.sdss_id_flat_pk__seq');
ALTER SEQUENCE catalogdb.sdss_id_flat_pk__seq OWNED BY catalogdb.sdss_id_flat.pk;
SELECT setval('sandbox.sdss_id_flat_pk__seq', (SELECT MAX(pk) from catalogdb.sdss_id_flat));
