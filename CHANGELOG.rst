.. This changelog uses releases: https://releases.readthedocs.io/en/latest/

==========
Change Log
==========

* Full implementation of all cartons for ``0.5.0``.
* Implement method to create a carton from a FITS file (see `.get_file_carton`).

* :release:`0.2.2 <2021-03-29>`
* Various changes to ``xmatch`` plan ``0.5.0``.
* Xmatch: Run ``ANALYZE`` on temporary table only when clustering.
* Xmatch: Exclude reject tables from ``extra_nodes``.
* Xmatch: No need to analyze if phase 3 is skipped.

* :release:`0.2.1 <2021-03-05>`
* Pin ``sdssdb==0.4.8``.

* :release:`0.2.0 <2021-03-05>`
* :feature:`39` Improvements for cross-match ``v0.5``. Implements ``run_id`` and catalogid with ``run_id`` bit shifting, use of ``ra_orig`` and ``dec_orig`` from TIC_v8 when ``posflag=gaia2``, reject extended sources when matching against TIC_v8, and several improvements to cross-matching performance.

* :release:`0.1.4 <2021-01-11>`
* Unpin ``healpy`` version.

* :release:`0.1.3 <2021-01-11>`
* Tag for release of ``ops_std_boss_tic``. Also includes changes and new cartons in preparation for ``v0.5``.

* :release:`0.1.2 <2020-12-09>`
* Add Yanny scrapper script.
* Add ``ops_apogee_stds`` carton and 0.1.2 target selection plan.

* :release:`0.1.1 <2020-11-10>`
* Restore HRD cut for CB UVEX. Add SRCNUM column to UVEX.
* Remove GUVCat column output in cartons not joined with GUVCat.
* :support:`28` Tidying up and applying carton name changes throughout.
* :support:`27` Carton-program alignment.
* :feature:`29` Add MWM Legacy IR2Opt carton.

* :release:`0.1.0 <2020-08-18>`
* Framework for implementing cartons against ``catalogdb`` and load them into ``targetdb``.
* All the cartons for target selection ``v0``.
* Xmatch: do not apply Q3C in phase 1.
* Xmatch: do not disable ``seqscan`` during phase 3.
* Xmatch: add several indexes to ``Catalog`` and the relational tables.
* XMatch: support weights to determine the join paths.
* Xmatch: Fix bug that prevented proper motions to be used in cross-matching.
* Xmatch: In phase 2, determine what table is larger and define ``q3c_join`` accordingly.
* Code to create `.create_sky_catalogue <sky catalogues>`.

* :release:`0.1.0-alpha.1 <2020-04-21>`
* Basic framework. Cross-matching tools work. Target selection tools still incomplete.
