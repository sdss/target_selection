#!/bin/bash

# this script is copied+adapted from the dr18 version, and may contain some redundant info.
# In DR19 we just add the v0 bhm cartons that were used during the plate era

# So most of this is DR18 related
# run some queries against the database, format the output, paste the result into 'Custom HTML' block in wordpress.

# This generates the HTML for the tables of cartons on the following page:
# https://test.sdss.org/dr18/bhm/programs-dr19/cartons-dr19/

#alias sdssdb='psql -h localhost -p 7502 -U sdss_user -d sdss5db'
alias sdssdb='psql -h operations.sdss.org -U sdss_user -d sdss5db'

WORKINGDIR=~/SDSSV/dr19/minidb
cd $WORKINGDIR

if [[ "true" == "Ntrue" ]]; then
    cd $WORKINGDIR
    # No need to redo this for DR19
    Q="SELECT c.carton,c.target_selection_plan,v.tag,
              CASE WHEN (c.carton ~ 'bhm_colr' OR c.carton ~ 'bhm_gua')  THEN 'ancillary'
                   WHEN (cc.program = 'bhm_filler' OR c.carton ~ 'bhm_csc')  THEN 'non-core'
                   WHEN (cc.program = 'bhm_spiders' AND (c.carton ~ '_supercosmos' OR c.carton ~ '_efeds_stragglers' OR c.carton ~ '_skymapperdr2' OR c.carton ~ '_gaiadr2' OR c.carton ~ '_sep'))  THEN 'non-core'
                   ELSE 'core' END as status,
              cc.program
    FROM minidb.dr19_carton as c
    join carton as cc
      on c.carton = cc.carton
    join targetdb.version as v
      on cc.version_pk = v.pk
    where c.carton ~ 'bhm'
      and v.plan = c.target_selection_plan
    order by c.carton;"
    
    sdssdb -c "$Q" > q_result_v0.5.3.txt


    gawk -v gen='"v0.5.3"' 'BEGIN {
    printf("<figure class=\"wp-block-table is-style-stripes\">\n\
    <table>\n\
    <thead>\n\
    <tr>\
    <th>carton name</th>\
    <th class=\"has-text-align-center\" data-align=\"center\">science <br>program</th>\
    <th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>plan</th>\
    <th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>tag</th>\
    <th class=\"has-text-align-center\" data-align=\"center\">status</th>\
    </tr>\
    \n</thead>\n<tbody>\n");}
    $1~/^bhm_/{
    split($1,c,"_"); prog=c[2]; PROG=toupper(prog);
    if (prog~/gua|colr/) {prog="ancillary";PROG="Ancillary programs";};
    printf("<tr><td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a></td><td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/dr18/bhm/programs/%s\">BHM %s</a></td><td class=\"has-text-align-center\" data-align=\"center\">%s</td><td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td><td class=\"has-text-align-center\" data-align=\"center\">%s</td></tr>\n",
    $1, $3, $1, prog, PROG, $3, $5, $5, $7)}
    END {printf("</tbody></table><figcaption>A table of all BHM cartons from targeting generation %s</figcaption></figure>\n", gen)}' q_result_v0.5.3.txt > carton_table_block_generation_0.5.3.html
    
fi



if [[ "true" == "Ntrue" ]]; then
    cd $WORKINGDIR
    ######################################
    ## manually select the BHM plates cartons
    # check this logic when have db access
    
    Q="SELECT cc.carton,v.plan,v.tag,
              CASE WHEN (c.carton ~ 'bhm_colr' OR c.carton ~ 'bhm_gua')  THEN 'ancillary'
                   WHEN (cc.program = 'bhm_filler' OR c.carton ~ 'bhm_csc' OR c.carton ~ '-efeds')  THEN 'non-core'
                   ELSE 'core' END as status,
              cc.program
    FROM carton as cc
    join targetdb.version as v
      on cc.version_pk = v.pk
    where
      cc.carton ~ 'bhm'
      and v.plan ~ '0.1.0'
    order by cc.carton;"
    
    sdssdb -c "$Q" > q_result_SDSSV_plates.txt
    
    
    gawk -v gen='"SDSS-V plates"' 'BEGIN {
    printf("<figure class=\"wp-block-table is-style-stripes\">\n\
    <table>\n\
    <thead>\n\
    <tr>\
    <th>carton name</th>\
    <th class=\"has-text-align-center\" data-align=\"center\">science <br>program</th>\
    <th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>plan</th>\
    <th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>tag</th>\
    <th class=\"has-text-align-center\" data-align=\"center\">status</th>\
    </tr>\
    \n</thead>\n<tbody>\n");}
    $1~/^bhm_/{
    split($1,c,"_"); prog=c[2]; PROG=toupper(prog);
    if (prog~/gua|colr/) {prog="ancillary";PROG="Ancillary programs";};
    printf("<tr><td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a></td><td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/dr18/bhm/programs/%s\">BHM %s</a></td><td class=\"has-text-align-center\" data-align=\"center\">%s</td><td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td><td class=\"has-text-align-center\" data-align=\"center\">%s</td></tr>\n",
    $1, $3, $1, prog, PROG, $3, $5, $5, $7)}
    END {printf("</tbody></table><figcaption>A table of all BHM cartons from targeting generation %s</figcaption></figure>\n", gen)}' q_result_SDSSV_plates.txt > carton_table_block_generation_SDSSV_plates.html

fi


### now process the json
# see: proc_carton_desc.py
# e.g.
WEBBASE=~/SDSSV/gitwork/target_selection/docs/website_dr19/bhm_cartons/
cd $WEBBASE

HTMLOUT=${WEBBASE}/bhm_target_cartons.html
python proc_carton_desc.py > $HTMLOUT

echo "HTML output is here: $HTMLOUT"


# Then do the auto conversion to latex:
#
TEXOUT=${WEBBASE}/bhm_target_cartons.tex
pandoc --from=html --to=latex --output=$TEXOUT $HTMLOUT

# replace the references with bibtex equivalents - need to revise this for DR19:

perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2020ApJS..250....8L\/abstract}{Lyke\net al., 2020}/\\citealt{Lyke2020}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2014ApJ...785..104R\/abstract}{Rykoff\net al., 2014}/\\citealt{Rykoff2014}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2020MNRAS.499.4768I\/abstract}{Ider\nChitham et al., 2020}/\\citealt{IderChitham2020}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022A\\%26A...661A...3S\/abstract}{Salvato\net al., 2022}/\\citealt{Salvato2022}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022A\\%26A...661A...3S\/abstract}{Salvato\net al. 2022}/\\citealt{Salvato2022}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2019MNRAS.489.4741S\/abstract}{Shu\net al., 2019}/\\citealt{Shu2019}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2019MNRAS.489.4741S\/abstract}{Shu\net al., \(2019\)}/\\citet{Shu2019}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2021ApJS..253....8M\/abstract}{CatWISE2020}/CatWISE2020 \\citep{Marocco2021}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022arXiv220608989Y\/abstract}{Yang\nand Shen, \(2022\)}/\\citet{Yang2022}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022arXiv220608989Y\/abstract}{Yang\nand Shen \(2022\)}/\\citet{Yang2022}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022arXiv220608989Y\/abstract}{Yang\nand Shen, 2022}/\\citealt{Yang2022}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2019PASA...36...33O\/abstract}{SkyMapper-dr2}/SkyMapper-dr2 \\citep{Onken2019}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022A\\%26A...661A...2L\/abstract}{Liu\net al., 2022}/\\citealt{Liu2022}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022ApJS..259...35A\/abstract}{Abdurro'\''uf\net al., 2022}/\\citealt{Abdurrouf_2021_sdssDR17}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2011ApJ...729..141B\/abstract}{Bovy\net al., 2011}/\\citealt{Bovy2011}/g'  $TEXOUT
# perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/xxx\/abstract}{xxx\net al., xxxx}/\\citealt{xxxxx}/g'  $TEXOUT

# get the section numbering right
perl -0777 -pi -e 's/\\subsubsection{/\\subsection{/g'  $TEXOUT

# formatting
perl -0777 -pi -e 's/\\textbf{/\\noindent\\textbf{/g'  $TEXOUT

#replace some spurious stuff
perl -0777 -pi -e 's/σ/\$\\sigma\$/g'  $TEXOUT
perl -0777 -pi -e 's/→/\$\\rightarrow\$/g'  $TEXOUT

# perl -0777 -pi -e 's///g'  $TEXOUT
awk 'BEGIN {flag=1} NF==0 {flag=1} $0~/Catalogdb tables required/ {flag=0} flag==1 {print $0}' $TEXOUT >  ${TEXOUT}.temp
mv ${TEXOUT}.temp  $TEXOUT

# checks (should result in zero rows)
grep -A1 adsab $TEXOUT

echo "Latex output is here: $TEXOUT"
