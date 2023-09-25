UPDATE sandbox.sdss_id_stacked sid
SET ra = cat.ra,
	dec = cat.dec
FROM catalog cat
WHERE sid.catalogid31=cat.catalogid
and sid.catalogid31 is not NULL;

UPDATE sandbox.sdss_id_stacked sid
SET ra = cat.ra,
	dec = cat.dec
FROM catalog cat
WHERE sid.catalogid25=cat.catalogid
and sid.catalogid31 is NULL
and sid.catalogid25 is not NULL;

UPDATE sandbox.sdss_id_stacked sid
SET ra = cat.ra,
	dec = cat.dec
FROM catalog cat
WHERE sid.catalogid21=cat.catalogid
and sid.catalogid31 is NULL
and sid.catalogid25 is NULL
and sid.catalogid21 is not NULL;
