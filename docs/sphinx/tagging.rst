Tagging ``target-selection``
============================

Target selection (and cross-matching) runs rely on code being tagged before the run is executed. In ``catalogdb`` and ``targetdb`` a cross-match or new carton are associated with a version that include the run string and the code tag. It is not possible to modify an already existing run with a carton associated with a different tag.

The steps to implement a change and tag a new version/run are as follows:

- Ensure that you are working on a pre-release version. In ``setup.cfg`` the version should be something like ``0.7.1-alpha.0``. If it's not, you'll need to bump the version in ``main``.
- In a branch or branches, implement the changes that you want and merge them back to ``main``.
- For testing the new cartons, see the section below.
- Create a new target selection plan in ``python/target_selection/config/target_selection.yml``. Refer to :ref:`target-selection` for detail on the format of the run. Make sure the run name is unique. The new run can include a full rerun of all the cartons or only the ones that you have updated.

For writing the target_selection.yml entry for the new plan, you can copy and older plan and then modify it appropriately.
For example, below is the target_selection.yml for the case where we have modified six cartons. Note that the target_selection plan 0.9.1 and the xmatch plan 0.9.0 do not need to be the same.::

  '0.9.1':
   xmatch_plan: 0.9.0
   cartons:
     - mwm_yso_pms_apogee
     - mwm_yso_pms_boss
     - mwm_ob_core
     - mwm_ob_cepheids
     - mwm_halo_bb_boss
     - mwm_halo_sm_boss
   open_fiber_path: $CATALOGDB_DIR/../open_fiber/0.5.0/
   schema: sandbox
   magnitudes:
     h: [catalog_to_tic_v8, tic_v8, twomass_psc.h_m]
     j: [catalog_to_tic_v8, tic_v8, twomass_psc.j_m]
     k: [catalog_to_tic_v8, tic_v8, twomass_psc.k_m]
     bp: [catalog_to_tic_v8, tic_v8, gaia_dr2_source.phot_bp_mean_mag]
     rp: [catalog_to_tic_v8, tic_v8, gaia_dr2_source.phot_rp_mean_mag]
     gaia_g: [catalog_to_tic_v8, tic_v8, gaia_dr2_source.phot_g_mean_mag]


- Update the CHANGELOG.md. You can do this while modifying ``main`` under the section ``## Next release``. When you are ready to tag, rename ``## Next release`` to ``## {tag} - {data}``.
- Bump the version in ``setup.cfg`` to the release version. For example, if you have been working as ``0.7.1-alpha.0`` you can bump to ``0.7.1`` but if you have made significant changes you'll want to do ``0.8.0`` or even ``1.0.0``. Follow the rules of semantic versioning.
- Although it should not be required, it's best to reinstall the package with ``pip install .`` or ``pip install -e .`` to ensure the metadata is updated.
- Tag the new version with the same name as the version in ``setup.cfg``. Push the changes to GitHub.
- *Immediately* after this, bump the version in ``setup.cfg`` to the next pre-release.
- A GitHub Action should create a new GitHub release and push the new version to PyPI. Ideally, update the GitHub release with the same contents as the changelog for the new version.
- At Utah, pull the changes, checkout the new tag, and do the new target selection run as shown below. This should create a new entry in ``targetdb.version`` that matches the run version and code tag.

Below is the command to run target_selection for the plan 0.9.1

cd bin

python target_selection -u sdss run 0.9.1 &

If you want to overwrite the results of the previous run of the above command then run the below command.

Note that some options go before the 'run' and some go after the 'run'.
 
python target_selection  -u sdss run --overwrite 0.9.1 &

The above command runs all the cartons listed in target_selection.yml for 0.9.1. 
If you only want to run one carton (e.g. mwm_xyz) then run the below command.

python target_selection  -u sdss run --overwrite --include mwm_xyz  0.9.1 &

Run the below command for more information about target_selection.

python target_selection run --help

Running test cartons
--------------------

It is convenient to be able to run and rerun cartons and ingest them to ``targetdb`` during testing. For this you can use the special ``0.0.0-test`` run. For that run the requirement that you cannot change the tag for an existing run is not enforced. Just add your new cartons to that run definition in ``target_selection.yml`` and run ``0.0.0-test``. The test run is meant as a sandbox and it can change at any time, so nothing should be run with it that needs to be relied upon in the future.

Alternatively, you can create the new run definition in the YAML file and work with it. The associated runs in ``targetdb`` will be associated with the pre-release code version. Once you have tagged, you can manually remove that version (along with its associated cartons) and execute the run again, which should then be inserted again with the correct code tag.
