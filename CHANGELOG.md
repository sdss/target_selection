# Changelog

## Next version

### Fixed

##  0.3.7 - November 19, 2021

### New

* Add the below new cartons.
ops_sky_apogee_best
ops_sky_apogee_good
ops_sky_boss_best
ops_sky_boss_good

##  0.3.6 - October 8, 2021

### New

* Add the below new cartons.
ops_tycho2_brightneighbors
ops_2mass_psc_brightneighbors
ops_gaia_brightneighbors

##  0.3.5 - September 27, 2021

### New

* Add the below cartons to target_selection.yml so they can be rerun.
bhm_aqmes_bonus_bright
bhm_aqmes_bonus_core
bhm_aqmes_bonus_faint
bhm_aqmes_wide2
bhm_aqmes_wide2_faint
mwm_cb_300pc_apogee
mwm_cb_300pc_boss

##  0.3.4 - September 6, 2021

### 🔧 Fixed

* [#101](https://github.com/sdss/target_selection/pull/101) Remove the `sleep` statements in `base.py`. This should now work fine after sdssdb PR [#99](https://github.com/sdss/sdssdb/pull/99).


## 0.3.3 - August 27, 2021

### ✨ Improvements

* Modify get_file_carton() in cartons/tools.py for open fiber cartons.
* Modify mwm_erosita_stars, mwm_erosita_compact_gen, and mwm_erosita_compact_var cartons for extra cases in assigning instrument, cadence, priority.


## 0.3.2 - July 21, 2021

### 🔧 Fixed

* Changes to address new columns not showing in the model of the carton temporary table.


## 0.3.1 - July 16, 2021

### 🚀 New

* Remove lower magnitude limit for mwm_yso, mwm_ob, and mwm_halo cartons for `0.5.0`.


## 0.3.0 - June 22, 2021

### 🚀 New

* Full implementation of all cartons for `0.5.0`.
* Implement method to create a carton from a FITS file (see `get_file_carton`).


## 0.2.2 - March 29, 2021

### ✨ Improvements

* Various changes to `xmatch` plan `0.5.0`.
* Xmatch: Run `ANALYZE` on temporary table only when clustering.
* Xmatch: Exclude reject tables from `extra_nodes`.
* Xmatch: No need to analyze if phase 3 is skipped.


## 0.2.1 - March 5, 2021

### 🧹 Cleanup

* Pin `sdssdb==0.4.8`.


## 0.2.0 - March 5, 2021

### 🚀 New

* [39](https://github.com/sdss/target_selection/pull/39) Improvements for cross-match `v0.5`. Implements `run_id` and catalogid with `run_id` bit shifting, use of `ra_orig` and `dec_orig` from TIC_v8 when `posflag=gaia2`, reject extended sources when matching against TIC_v8, and several improvements to cross-matching performance.


## 0.1.4 - January 11, 2021

### 🧹 Cleanup

* Unpin `healpy` version.


## 0.1.3 - January 11, 2021

### 🧹 Cleanup

* Tag for release of `ops_std_boss_tic`. Also includes changes and new cartons in preparation for `v0.5`.


## 0.1.2 - December 9, 2020

### 🚀 New

* Add Yanny scrapper script.
* Add `ops_apogee_stds` carton and 0.1.2 target selection plan.


## 0.1.1 - November 10, 2020

### ✨ Improvements

* Restore HRD cut for CB UVEX. Add SRCNUM column to UVEX.
* Remove GUVCat column output in cartons not joined with GUVCat.
* [28](https://github.com/sdss/target_selection/pull/28) Tidying up and applying carton name changes throughout.
* [27](https://github.com/sdss/target_selection/pull/27) Carton-program alignment.
* [29](https://github.com/sdss/target_selection/pull/29) Add MWM Legacy IR2Opt carton.


## 0.1.0 - August 18, 2020

### 🚀 New

* Framework for implementing cartons against `catalogdb` and load them into `targetdb`.
* All the cartons for target selection `v0`.
* Xmatch: do not apply Q3C in phase 1.
* Xmatch: do not disable `seqscan` during phase 3.
* Xmatch: add several indexes to `Catalog` and the relational tables.
* XMatch: support weights to determine the join paths.
* Xmatch: Fix bug that prevented proper motions to be used in cross-matching.
* Xmatch: In phase 2, determine what table is larger and define `q3c_join` accordingly.
* Code to create `.create_sky_catalogue <sky catalogues>`.


## 0.1.0-alpha.1 - April 21, 2020

### 🚀 New

* Basic framework. Cross-matching tools work. Target selection tools still incomplete.
