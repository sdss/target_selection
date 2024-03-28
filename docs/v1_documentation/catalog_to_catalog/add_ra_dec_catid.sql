create table sandbox.sdss_id_flat_new as
select sid.sdss_id, sid.catalogid, sid.version_id, sid.ra_sdss_id, sid.dec_sdss_id,
       sid.n_associated, sid.pk, cat.ra as ra_catalogid, cat.dec as dec_catalogid
from catalogdb.sdss_id_flat sid 
join catalogdb.catalog cat
on sid.catalogid=cat.catalogid;
