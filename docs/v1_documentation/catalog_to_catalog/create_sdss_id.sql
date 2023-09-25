-- Process to create sdss_id_stacked --

-- Collect all catalog ids (from targetdb.target and catalogidx_to_catalogidy) into temp tables in each version
-- We just want a list of all possible catalogids from each version
CREATE TABLE sandbox.temp_catalogid_v21 as (
	select cat.catalogid as catalogid21 from catalog cat join target tar on tar.catalogid=cat.catalogid where cat.version_id=21);
INSERT INTO sandbox.temp_catalogid_v21 (catalogid21)
	select cc.catalogid1 as catalogid21 from sandbox.catalogidx_to_catalogidy_unique cc where cc.version_id1=21;
INSERT INTO sandbox.temp_catalogid_v21 (catalogid21)
	select cc.catalogid2 as catalogid21 from sandbox.catalogidx_to_catalogidy_unique cc where cc.version_id2=21;

CREATE TABLE sandbox.temp_catalogid_v25 as (
	select cat.catalogid as catalogid25 from catalog cat join target tar on tar.catalogid=cat.catalogid where cat.version_id=25);
INSERT INTO sandbox.temp_catalogid_v25 (catalogid25)
	select cc.catalogid1 as catalogid25 from sandbox.catalogidx_to_catalogidy_unique cc where cc.version_id1=25;
INSERT INTO sandbox.temp_catalogid_v25 (catalogid25)
	select cc.catalogid2 as catalogid25 from sandbox.catalogidx_to_catalogidy_unique cc where cc.version_id2=25;

CREATE TABLE sandbox.temp_catalogid_v31 as (
	select cat.catalogid as catalogid31 from catalog cat join target tar on tar.catalogid=cat.catalogid where cat.version_id=31);
INSERT INTO sandbox.temp_catalogid_v31 (catalogid31)
	select cc.catalogid1 as catalogid31 from sandbox.catalogidx_to_catalogidy_unique cc where cc.version_id1=31;
INSERT INTO sandbox.temp_catalogid_v31 (catalogid31)
	select cc.catalogid2 as catalogid31 from sandbox.catalogidx_to_catalogidy_unique cc where cc.version_id2=31;

-- Make the stacked sdss_id table using the list of catalogids above
-- Ids are paired using the catalogidx_to_catalogidy table only on subsequent versions (21 -> 25 -> 31)

DROP TABLE sandbox.sdss_id_stacked;
CREATE TABLE sandbox.sdss_id_stacked as 
(with 
cat_ids21 as (select distinct tc.catalogid21, cat.lead from sandbox.temp_catalogid_v21 tc join catalog cat on cat.catalogid=tc.catalogid21),
cat_ids25 as (select distinct tc.catalogid25, cat.lead from sandbox.temp_catalogid_v25 tc join catalog cat on cat.catalogid=tc.catalogid25),
cat_ids31 as (select distinct tc.catalogid31, cat.lead from sandbox.temp_catalogid_v31 tc join catalog cat on cat.catalogid=tc.catalogid31),
sq21_25 as (select cc.catalogid1, cc.catalogid2 from sandbox.catalogidx_to_catalogidy_all cc 
			join catalog cat on cat.catalogid=cc.catalogid1 and cat.lead=cc.lead
			where cc.version_id1=21 and cc.version_id2=25),
sq25_31 as (select cc.catalogid1, cc.catalogid2 from sandbox.catalogidx_to_catalogidy_all cc 
			join catalog cat on cat.catalogid=cc.catalogid1 and cat.lead=cc.lead
			where cc.version_id1=25 and cc.version_id2=31),
left21_to_25 as (select cat_ids21.catalogid21, sq21_25.catalogid2 as catalogid25 from cat_ids21 left join sq21_25 on cat_ids21.catalogid21=sq21_25.catalogid1),
add_outer25 as (select left21_to_25.catalogid21, cat_ids25.catalogid25 as catalogid25 from left21_to_25 full outer join cat_ids25 on left21_to_25.catalogid25=cat_ids25.catalogid25),
left_25_to_31 as (select add_outer25.catalogid21, add_outer25.catalogid25, sq25_31.catalogid2 as catalogid31 from add_outer25 left join sq25_31 on add_outer25.catalogid25=sq25_31.catalogid1),
add_outer31 as (select left_25_to_31.catalogid21, left_25_to_31.catalogid25, cat_ids31.catalogid31 from left_25_to_31 full outer join cat_ids31 on left_25_to_31.catalogid31=cat_ids31.catalogid31)
select * from add_outer31);
ALTER TABLE sandbox.sdss_id_stacked ADD COLUMN sdss_id SERIAL PRIMARY KEY;

-- Clean up

DROP TABLE sandbox.temp_catalogid_v21;
DROP TABLE sandbox.temp_catalogid_v25;
DROP TABLE sandbox.temp_catalogid_v31;


--  Process to create sdss_id_flat --


DROP TABLE sandbox.sdss_id_flat;

CREATE TABLE sandbox.sdss_id_flat (sdss_id bigint, catalogid bigint);

INSERT INTO sandbox.sdss_id_flat (sdss_id, catalogid)
select sdss_id, catalogid21 as catalogid from sandbox.sdss_id_stacked where catalogid21 IS NOT NULL;

INSERT INTO sandbox.sdss_id_flat (sdss_id, catalogid)
select sdss_id, catalogid25 as catalogid from sandbox.sdss_id_stacked where catalogid25 IS NOT NULL;

INSERT INTO sandbox.sdss_id_flat (sdss_id, catalogid)
select sdss_id, catalogid31 as catalogid from sandbox.sdss_id_stacked where catalogid31 IS NOT NULL;
