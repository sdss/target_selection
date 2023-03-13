# Changelog

## 1.0.7 - March 13, 2023

### Fixed

* Update target_selection.yml for plan 1.0.7 for manual cartons

## 1.0.6 - March 10, 2023

### Fixed

* Update target_selection.yml for plan 1.0.6 cartons

## 1.0.5 - March 9, 2023

### Fixed

* Update target_selection.yml for plan 1.0.5 cartons

* Updated code for gaia dr3 XP synthetic magnitudes

## 1.0.4 - March 4, 2023

### Fixed

* Use gaia dr3 XP synthetic mags as supplier of griz for targets

## 1.0.3 - February 10, 2023

### Fixed

* Update for manual cartons can_offset


## 1.0.2 - February 8, 2023

### Fixed

* Update tools.py for manual cartons for Gaia DR3

## 1.0.1 - February 8, 2023

### Fixed

* Update mwm_yso cartons for v1.0

## 1.0.0 - February 6, 2023

### 🚀 New

* This tag includes the full code used for the cross-match for v1. Target selection code for v1 cartons **is not complete** as of this tag.


## 0.3.23 - November 15, 2022

### ✨ Improved

* Added handling of `can_offset` column for file cartons.


## 0.3.22 - August 6, 2022

### 🔧 Fixed

* tools.py: add Gaia_DR3_Source_ID for manual/openfiber FITS file.


## 0.3.21 - August 4, 2022

### 🔧 Fixed

* ops_std_apogee: add Gaia proper motion and parallax cut


## 0.3.20 - July 5, 2022

### 🔧 Fixed

* ops_std_apogee: remove the condition (TwoMassPSC.j_m - TwoMassPSC.k_m) < 0.5


## 0.3.19 - May 24, 2022

### 🔧 Fixed

* tools.py: valid_program: add mwm_validation


## 0.3.18 - February 11, 2022

### 🔧 Fixed

* mwm_cb.py: MWM_CB_300_Carton: add condition FUV > -999
* The above change affects the below cartons.
  * mwm_cb_300pc_apogee
  * mwm_cb_300pc_boss


## 0.3.17 - February 8, 2022

### 🔧 Fixed

* Update priority for the below cartons.
  * mwm_rv_long_fps
  * mwm_rv_short_fps


## 0.3.16 - February 2, 2022

### 🔧 Fixed

* Update priority to 2705 for the below cartons.
  * mwm_yso_cluster_apogee
  * mwm_yso_cluster_boss
  * mwm_yso_cmz_apogee
  * mwm_yso_disk_apogee
  * mwm_yso_disk_boss
  * mwm_yso_embedded_apogee
  * mwm_yso_nebula_apogee
  * mwm_yso_variable_apogee
  * mwm_yso_variable_boss


## 0.3.15 - January 31, 2022

### 🔧 Fixed

* Update below carton.
  * bhm_colr_galaxies_lsdr8


## 0.3.14 - January 27, 2022

### 🔧 Fixed

* Update below cartons.
  * bhm_csc_boss
  * bhm_csc_apogee


## 0.3.13 - January 25, 2022

### 🚀 New

* Update below cartons.
  * bhm_csc_boss
  * bhm_csc_apogee


## 0.3.12 - January 20, 2022

### 🚀 New

* New open fiber cartons FITS format.


## 0.3.11 - January 5, 2022

### 🚀 New

* Add the below new carton.
  * ops_sky_boss_fallback
* Update the priority variable in target_selection.yml for the below carton.
  * bhm_colr_galaxies_lsdr8


## 0.3.10 - December 6, 2021

### 🔧 Fixed

* Update target_selection.yml for the below carton.
  * mwm_tess_ob


## 0.3.9 - November 23, 2021

### 🔧 Fixed

* Change program variable of the below cartons.
  * ops_sky_apogee_best
  * ops_sky_apogee_good
  * ops_sky_boss_best
  * ops_sky_boss_good


## 0.3.8 - November 20, 2021

### 🔧 Fixed

* Change category variable of the below cartons.
  * ops_sky_apogee_best
  * ops_sky_apogee_good
  * ops_sky_boss_best
  * ops_sky_boss_good


## 0.3.7 - November 19, 2021

### 🚀 New

* Add the below new cartons.
  * ops_sky_apogee_best
  * ops_sky_apogee_good
  * ops_sky_boss_best
  * ops_sky_boss_good


## 0.3.6 - October 8, 2021

### 🚀 New

* Add the below new cartons.
  * ops_tycho2_brightneighbors
  * ops_2mass_psc_brightneighbors
  * ops_gaia_brightneighbors


## 0.3.5 - September 27, 2021

### 🚀 New

* Add the below cartons to target_selection.yml so they can be rerun.
  * bhm_aqmes_bonus_bright
  * bhm_aqmes_bonus_core
  * bhm_aqmes_bonus_faint
  * bhm_aqmes_wide2
  * bhm_aqmes_wide2_faint
  * mwm_cb_300pc_apogee
  * mwm_cb_300pc_boss


## 0.3.4 - September 6, 2021

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
