.. This changelog uses releases: https://releases.readthedocs.io/en/latest/

==========
Change Log
==========

:release:`0.1.0 <2020-08-18>`
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
