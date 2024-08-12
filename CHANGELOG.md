# Changelog

## Next version

## 1.3.13 - August 10, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.14

## 1.3.12 - August 1, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.13

## 1.3.11 - August 1, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.12

## 1.3.10 - July 30, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.11

### 🔧 Fixed

* [#481](https://github.com/sdss/target_selection/pull/481) Explicitely force the reflection cache to be update after the carton temporary table has been created.


### 📖 Documentation

* [#480](https://github.com/sdss/target_selection/pull/480) Update `README.md` with instructions for parallel development and tagging of `sdssdb` and `target_selection`.


## 1.3.9 - July 29, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.10


## 1.3.8 - July 29, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.9


## 1.3.7 - July 28, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.8


## 1.3.6 - July 27, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.7


## 1.3.5 - July 24, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.6


## 1.3.4 - July 23, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.5
* pyproject.toml: edit include parameter


## 1.3.3 - July 21, 2024

### 🚀 New

* add cartons to target_selection plan 1.2.4
* added bhm\_aqmes\_wide1 and bhm\_aqmes\_wide1\_faint cartons
  * identical to 'wide2' versions, but with a 'dark\_1x4' cadence and +1 added to priorities
* new '\_d3' versions of selected BHM cartons
  * these new cartons assign targets which prevously had cadence 'dark\_flexible\_4x1' a new cadence of 'dark\_flexible\_3x1'
  * these new versions only output the subset of targets having cadence (targets with easier cadences are omitted)
  * the '\_d3' versions have target priority incremented by 1
  * new '\_d3' cartons are:
    * bhm\_csc\_boss\_d3
    * bhm\_gua\_dark\_d3
    * bhm\_spiders\_clusters\_lsdr10\_d3
    * bhm\_spiders\_agn\_lsdr10\_d3
    * bhm\_spiders\_agn\_hard\_d3
    * bhm\_spiders\_agn\_gaiadr3\_d3
    * bhm\_spiders\_agn\_tda\_d3
    * bhm\_spiders\_agn\_sep\_d3
    * bhm\_colr\_galaxies\_lsdr10\_d3
    * mwm\_erosita\_compact\_boss\_d3
  * incremented priority of mwm\_erosita\_compact\_boss\_shallow to avoid conflicts with '\_d3'


## 1.3.2 - July 17, 2024

### ✨ Improved

* Add cartons to target_selection plan 1.2.3


## 1.3.1 - July 12, 2024

### ✨ Improved

* Add cartons to target_selection plan 1.2.2


## 1.3.0 - July 10, 2024

### ✨ Improved

* Add cartons to target_selection plan 1.2.1

* [#455](https://github.com/sdss/target_selection/pull/455) Significant rewrite of the `FileCarton` code to deal with instances of duplicate rows for a single input identifier.

### ⚙️ Engineering

* [#454](https://github.com/sdss/target_selection/pull/454) Adapt `target_selection` to using a PEP517 installer with Poetry backend. Lint and format using `ruff`. Added linting and release workflows and updated the RTDs build.


## 1.2.7 - June 20, 2024

### 🚀 New

* tools.py: add column LegacySurvey_DR10_ID for manual carton fits file

## 1.2.6 - June 3, 2024

### 🚀 New

* changed the cadence to bright_2x1_long for the cartons mwm_bin_gaia_astb_apogee and mwm_bin_gaia_astb_boss.


## 1.2.5 - May 17, 2024

### 🚀 New

* Added cross-match plan `1.1.4` for `mangatarget`
* Added cross-match plan `1.1.5` for `mastar_goodstars`


## 1.2.4 - May 7, 2024

### 🚀 New

* Added cross-match plan `1.1.3` for `sdss_dr17_specobj`


## 1.2.3 - May 2, 2024

### 🚀 New

* Added cross-match plan `1.1.2` for `marvels_dr12_star`


## 1.2.2 - April 30, 2024

### ✨ Improved

* Changed some warning to log messages.

### 🔧 Fixed

* Removed leftover print message.


## 1.2.1 - April 25, 2024

### 🚀 New

* Added cross-match plan `1.1.1` for `marvels_dr11_star`.
* A few additional improvements to the `XMatchPlanner` code.


## 1.2.0 - April 23, 2024

### ✨ Improved

This version makes several improvements to `XMatch`:

* Removes some hard-coded use of SQL files. These should not be necessary and if they were, we should get to the bottom of why they are ...
* Store intermediate results in a sandboxed `catalog_to_XXX` table. This allows to run a full dry-run cross-match for a single catalogue.
* Phase 1: distinct only on `model_pk`.
* Phase 3: Require using only best matches.

Additionally, this version adds the `too` carton.


## 1.1.1 - October 13, 2023

### 🔧 Fixed

* mwm_rv.py: mwm_bin_rv_short_rgb_apogee: remove downsampling.
* Update target_selection.yml for target_selection_plan 1.0.51.


## 1.1.0 - September 28, 2023

### 🚀 New

* Added support in the cross-match code to run "addendum" cross-match runs.
* Configuration for cross-match run 1.1.0.


## 1.0.50 - August 20, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.50 cartons


## 1.0.49 - August 17, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.49 cartons


## 1.0.48 - August 8, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.48 cartons


## 1.0.47 - August 3, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.47 cartons


## 1.0.46 - August 2, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.46 manual cartons


## 1.0.45 - July 31, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.45


## 1.0.44 - July 28, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.44


## 1.0.43 - July 26, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.43


## 1.0.42 - July 24, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.42


## 1.0.41 - July 23, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.41
* Update mwm_cb_*, mwm_bin_vis* carton names + priorities


## 1.0.40 - July 21, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.40
* Update mwm_erosita_* carton names + priorities


## 1.0.39 - July 20, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.39
* Update mwm_bin_rv*, mwm_ob*, mwm_halo* cartons


## 1.0.38 - July 12, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.38


## 1.0.37 - July 7, 2023

### 🔧 Fixed

* Modify several BHM (+MWM-eROSITA) cadences to be 'flexible' equivalents.
* Adjust method by which AQMES cartons aquire cadence choices
* Update target_selection.yml for plan 1.0.37


## 1.0.36 - July 4, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.36
* change priority for mwm_ob_core and mwm_ob_cepheids


## 1.0.35 - July 2, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.35 manual cartons


## 1.0.34 - June 30, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.34 cartons
* add cartons: mwm_bin_rv_short_mdwarf, mwm_bin_rv_short_subgiant, mwm_bin_rv_short_rgb
* remove carton mwm_bin_rv_short


## 1.0.33 - June 14, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.33 cartons
* Add mwm_yso_*_single cartons


## 1.0.32 - June 1, 2023

## Fixed

* Rerun ``mwm_magcloud_agb_apogee``, ``mwm_magcloud_rgb_boss``, ``mwm_cb_galex_vol``, ``mwm_cb_xmmom``, and ``mwm_cb_swiftuvot``.


## 1.0.31 - May 23, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.31 cartons


## 1.0.30 - May 21, 2023

### 🔧 Fixed

* Rerun ``mwm_cb_galex_vol``, ``mwm_cb_xmmom``, and ``mwm_cb_swiftuvot``.


## 1.0.29 - May 18, 2023

### 🔧 Fixed

* Rerun ``ops_std_boss_ps1dr2``.


## 1.0.28 - May 16, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.28 cartons


## 1.0.27 - May 9, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.27 cartons


## 1.0.26 - May 8, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.26 cartons


## 1.0.25 - May 8, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.25 cartons


## 1.0.24 - May 7, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.24 cartons


## 1.0.23 - May 6, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.23 cartons


## 1.0.22 - May 5, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.22 cartons


## 1.0.21 - May 2, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.21 cartons


## 1.0.20 - May 2, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.20 cartons


## 1.0.19 - April 26, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.19 cartons


## 1.0.18 - April 17, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.18 cartons


## 1.0.17 - April 15, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.17 cartons


## 1.0.16 - April 9, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.16 cartons


## 1.0.15 - April 8, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.15 cartons


## 1.0.14 - April 7, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.14 cartons


## 1.0.13 - April 6, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.13 cartons


## 1.0.12 - April 3, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.12 cartons


## 1.0.11 - March 31, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.11 cartons


## 1.0.10 - March 28, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.10 for manual cartons


## 1.0.9 - March 28, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.9 cartons


## 1.0.8 - March 17, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.8 cartons


## 1.0.7 - March 13, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.7 for manual cartons


## 1.0.6 - March 10, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.6 cartons


## 1.0.5 - March 9, 2023

### 🔧 Fixed

* Update target_selection.yml for plan 1.0.5 cartons
* Updated code for gaia dr3 XP synthetic magnitudes


## 1.0.4 - March 4, 2023

### 🔧 Fixed

* Use gaia dr3 XP synthetic mags as supplier of griz for targets


## 1.0.3 - February 10, 2023

### 🔧 Fixed

* Update for manual cartons can_offset


## 1.0.2 - February 8, 2023

### 🔧 Fixed

* Update tools.py for manual cartons for Gaia DR3


## 1.0.1 - February 8, 2023

### 🔧 Fixed

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
