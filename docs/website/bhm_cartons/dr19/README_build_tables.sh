#!/usr/bin/env sh

# This generates the HTML for the tables of cartons on the following page:
# https://testng.sdss.org/${DR}/bhm/programs/cartons/

# run some queries against the database, format the output, paste the result into 'Custom HTML' clock in wordpress.

DR=dr19
BASEDIR=~/SDSSV/gitwork/target_selection_again/docs/website/bhm_cartons/$DR

cd $BASEDIR

# MINIDB=minidb
MINIDB=minidb_$DR

# sdss5db=> select * from minidb_dr19.dr19_targeting_generation;
#  pk |         label         | first_release
# ----+-----------------------+---------------
#   1 | v0.plates             | dr19
#   2 | v0.5.epsilon-7-core-0 | dr19
#   3 | v0.5.2                | dr19
#   4 | v0.5.3                | dr18
#   5 | v0.5.5                | dr19
# (5 rows)

GEN_LIST="v0.plates v0.5.epsilon-7-core-0 v0.5.2 v0.5.3 v0.5.5"
GEN_LIST="v0.plates"

# change the below to your favoured output location
cd ~/SDSSV/gitwork/target_selection/docs/website/bhm_cartons/${DR}

# this one when running from home
alias sdssdb='psql -h localhost -p 7502 -U sdss_user -d sdss5db'

# this one when running on a Utah machine
#alias sdssdb='psql -h operations.sdss.org -U sdss_user -d sdss5db'

Q="\copy ( SELECT 
          c.carton,v.plan,v.tag,
          CASE WHEN (c.carton ~ 'bhm_colr' OR c.carton ~ 'bhm_gua')  THEN 'ancillary'
               WHEN (c.program = 'bhm_filler' OR c.carton ~ 'bhm_csc')  THEN 'non-core'
               WHEN (c.program = 'bhm_spiders' AND (c.carton ~ '_supercosmos' OR c.carton ~ '_efeds_stragglers' OR c.carton ~ '_skymapperdr2' OR c.carton ~ 'agn_gaia' OR c.carton ~ '_sep'))  THEN 'non-core'
               ELSE 'core' END as status,
          CASE WHEN (c.program ~ 'gua' OR c.program ~ 'colr') THEN 'ancillary' ELSE c.program END as program, 
          tg.label
FROM ${MINIDB}.${DR}_targeting_generation as tg
JOIN ${MINIDB}.${DR}_targeting_generation_to_carton as tg2c
ON tg.pk = tg2c.generation_pk
join ${MINIDB}.${DR}_carton as c
on tg2c.carton_pk = c.carton_pk 
join  ${MINIDB}.${DR}_targetdb_version as v
  on c.version_pk = v.pk
where c.carton ~ 'bhm'
  and v.plan = c.target_selection_plan
order by c.carton_pk)
to 'q_result.csv' with csv header"

sdssdb -c "$Q"

# this builds a table per targeting generation

for GEN in $GEN_LIST; do 
    gawk --field-separator=',' \
         -v dr=$DR \
         -v gen="$GEN" \
'BEGIN {
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
$1~/^bhm_/ && $6 == gen {
    prog=$5; PROG=toupper(prog);  gsub("bhm_", "", prog); tag=$3;plan=$2;
    if (prog=="filler") {prog="ancillary"; PROG="BHM Ancillary programs";};
    printf("<tr>\
<td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/%s/bhm/programs/%s\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/blob/%s/python/target_selection/config/target_selection.yml\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
</tr>\n",
    $1, plan, $1, dr, prog, PROG, plan, plan, tag, tag, $4)}
END {
printf("</tbody></table><figcaption>All BHM cartons from targeting generation \"%s\"</figcaption></figure>\n", gen);
}' q_result.csv > carton_table_block_generation_${GEN}.html

done



# Goup similar targeting generations in same table
# membership of targeting generations is done as set of ticks in bool columns

Q="\copy ( SELECT 
          c.carton,v.plan,v.tag,
          CASE WHEN (c.carton ~ 'bhm_colr' OR c.carton ~ 'bhm_gua')  THEN 'ancillary'
               WHEN (c.program = 'bhm_filler' OR c.carton ~ 'bhm_csc')  THEN 'non-core'
               WHEN (c.program = 'bhm_spiders' AND (c.carton ~ '_supercosmos' OR c.carton ~ '_efeds_stragglers' OR c.carton ~ '_skymapperdr2' OR c.carton ~ 'agn_gaia' OR c.carton ~ '_sep'))  THEN 'non-core'
               ELSE 'core' END as status,
          CASE WHEN (c.program ~ 'gua' OR c.program ~ 'colr') THEN 'ancillary' ELSE c.program END as program, 
          SUM(CASE WHEN tg.label = 'v0.5.epsilon-7-core-0' THEN 1 ELSE 0 END) in_tg1,
          SUM(CASE WHEN tg.label = 'v0.5.2' THEN 1 ELSE 0 END) in_tg2,
          SUM(CASE WHEN tg.label = 'v0.5.3' THEN 1 ELSE 0 END) in_tg3,
          SUM(CASE WHEN tg.label = 'v0.5.5' THEN 1 ELSE 0 END) in_tg4,
          ARRAY_AGG(distinct tg.label) as targeting_generations
FROM ${MINIDB}.${DR}_targeting_generation as tg
JOIN ${MINIDB}.${DR}_targeting_generation_to_carton as tg2c
ON tg.pk = tg2c.generation_pk
join ${MINIDB}.${DR}_carton as c
on tg2c.carton_pk = c.carton_pk 
join  ${MINIDB}.${DR}_targetdb_version as v
  on c.version_pk = v.pk
where c.carton ~ 'bhm'
  and v.plan = c.target_selection_plan
group by c.carton_pk,c.carton,v.plan,v.tag,c.program 
order by c.carton_pk)
to 'q2_result.csv' with csv header"

sdssdb -c "$Q"


#toomuchdetail gawk --field-separator=',' \
#toomuchdetail          -v dr=$DR \
#toomuchdetail          -v gens="v0.5.epsilon-7-core-0 v0.5.2 v0.5.3 v0.5.5" \
#toomuchdetail 'BEGIN {
#toomuchdetail     split(gens,agens," ");
#toomuchdetail     printf("<figure class=\"wp-block-table is-style-stripes\">\n\
#toomuchdetail <table>\n\
#toomuchdetail <thead>\n\
#toomuchdetail <tr>\
#toomuchdetail <th>carton name</th>\
#toomuchdetail <th class=\"has-text-align-center\" data-align=\"center\">science <br>program</th>\
#toomuchdetail <th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>plan</th>\
#toomuchdetail <th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>tag</th>\
#toomuchdetail <th class=\"has-text-align-center\" data-align=\"center\">status</th>\
#toomuchdetail <th class=\"has-text-align-center\" data-align=\"center\">%s</th>\
#toomuchdetail <th class=\"has-text-align-center\" data-align=\"center\">%s</th>\
#toomuchdetail <th class=\"has-text-align-center\" data-align=\"center\">%s</th>\
#toomuchdetail <th class=\"has-text-align-center\" data-align=\"center\">%s</th>\
#toomuchdetail </tr>\
#toomuchdetail \n</thead>\n<tbody>\n",
#toomuchdetail agens[1], agens[2], agens[3], agens[4]);}
#toomuchdetail $1~/^bhm_/ && $NF~/v0.5/ {
#toomuchdetail     prog=$5; PROG=toupper(prog); gsub("bhm_", "", prog); tag=$3;plan=$2;
#toomuchdetail     if (prog=="filler") {prog="ancillary"; PROG="BHM Ancillary programs";};
#toomuchdetail     printf("<tr>\
#toomuchdetail <td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a></td>\
#toomuchdetail <td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/%s/bhm/programs/%s\">%s</a></td>\
#toomuchdetail <td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/blob/%s/python/target_selection/config/target_selection.yml\">%s</a></td>\
#toomuchdetail <td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td>\
#toomuchdetail <td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
#toomuchdetail <td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
#toomuchdetail <td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
#toomuchdetail <td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
#toomuchdetail <td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
#toomuchdetail </tr>\n",
#toomuchdetail     $1, plan, $1, dr, prog, PROG, plan, plan, tag, tag, $4, 
#toomuchdetail     ($6 > 0 ? "&#9989;" : ""), 
#toomuchdetail     ($7 > 0 ? "&#9989;" : ""), 
#toomuchdetail     ($8 > 0 ? "&#9989;" : ""), 
#toomuchdetail     ($9 > 0 ? "&#9989;" : "") )}
#toomuchdetail END {
#toomuchdetail printf("</tbody></table><figcaption>All BHM cartons from targeting generation \"%s\"</figcaption></figure>\n", gen);
#toomuchdetail }' q2_result.csv > carton_table_block_multi_generation_v0.5.html


gawk --field-separator=',' \
         -v dr=$DR \
         -v gens="v0.5.epsilon-7-core-0 and v0.5.2/v0.5.3/v0.5.5" \
'BEGIN {
    split(gens,agens," and ");
    printf("<figure class=\"wp-block-table is-style-stripes\">\n\
<table>\n\
<thead>\n\
<tr>\
<th>carton name</th>\
<th class=\"has-text-align-center\" data-align=\"center\">science <br>program</th>\
<th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>plan</th>\
<th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>tag</th>\
<th class=\"has-text-align-center\" data-align=\"center\">status</th>\
<th class=\"has-text-align-center\" data-align=\"center\">%s</th>\
<th class=\"has-text-align-center\" data-align=\"center\">%s</th>\
</tr>\
\n</thead>\n<tbody>\n",
agens[1], agens[2]);}
$1~/^bhm_/ && $NF~/v0.5/ {
    prog=$5; PROG=toupper(prog);  gsub("bhm_", "", prog); tag=$3;plan=$2;
    if (prog=="filler") {prog="ancillary"; PROG="BHM Ancillary programs";};
    printf("<tr>\
<td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/%s/bhm/programs/%s\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/blob/%s/python/target_selection/config/target_selection.yml\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
<td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
<td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
</tr>\n",
    $1, plan, $1, dr, prog, PROG, plan, plan, tag, tag, $4, 
    ($6 > 0 ? "&#9989;" : ""), 
    ($7 > 0 || $8 > 0 || $9 > 0 ? "&#9989;" : "") )}
END {
printf("</tbody></table><figcaption>All BHM cartons from targeting generations %s</figcaption></figure>\n", gens);
}' q2_result.csv > carton_table_block_multi_generation_v0.5.html



# ######################################
# ## manually select the eFEDS plates cartons
# 
# Q="SELECT cc.carton,v.plan,v.tag,
#           'non-core' as status,
#           cc.program
# FROM carton as cc
# join targetdb.version as v
#   on cc.version_pk = v.pk
# where
#   cc.carton ~ '-efeds'
#   and v.plan ~ '0.1.0'
# order by cc.carton;"
# 
# sdssdb -c "$Q" > q_result_eFEDS_plates.txt
# 
# 
# gawk -v gen='"eFEDS plates"' 'BEGIN {
# printf("<figure class=\"wp-block-table is-style-stripes\">\n\
# <table>\n\
# <thead>\n\
# <tr>\
# <th>carton name</th>\
# <th class=\"has-text-align-center\" data-align=\"center\">science <br>program</th>\
# <th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>plan</th>\
# <th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>tag</th>\
# <th class=\"has-text-align-center\" data-align=\"center\">status</th>\
# </tr>\
# \n</thead>\n<tbody>\n");}
# $1~/^bhm_/{
# split($1,c,"_"); prog=c[2]; PROG=toupper(prog);
# if (prog~/gua|colr/) {prog="ancillary";PROG="Ancillary programs";};
# printf("<tr><td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a></td><td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/dr18/bhm/programs/%s\">BHM %s</a></td><td class=\"has-text-align-center\" data-align=\"center\">%s</td><td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td><td class=\"has-text-align-center\" data-align=\"center\">%s</td></tr>\n",
# $1, $3, $1, prog, PROG, $3, $5, $5, $7)}
# END {printf("</tbody></table><figcaption>A table of all BHM cartons from targeting generation %s</figcaption></figure>\n", gen)}' q_result_eFEDS_plates.txt > carton_table_block_generation_eFEDS_plates.html
# 


### now process the json
# see: proc_carton_desc.py
# e.g.
python proc_carton_desc_dr19.py > bhm_target_cartons.html

# now check for missing/surplus cartons
gawk 'ARGIND==1 && /h3 id/ {\
    split($0,a,"\"");\
    f=a[2];\
    r[f]++;\
    n[f]=0;} \
ARGIND>=2 && /href="#/ {\
    split($0,a,"\"");\
    f=a[2];\
    gsub("#","",f); \
    n[f]++;\
    ok="missing"; \
    if(r[f]>0){ok="ok"}; \
    printf("%-50s %s\n", f, ok)}\
END {\
    for (f in n){\
       if(n[f]!=1){\
           printf("%-50s %d!=1\n", f, n[f]);\
    }}\
}' bhm_target_cartons.html carton_table_block_generation_v0.plates.html carton_table_block_multi_generation_v0.5.html


# paste into HTML box in wordpress - these comands help
xclip -sel c < carton_table_block_generation_v0.plates.html
xclip -sel c < carton_table_block_multi_generation_v0.5.html
xclip -sel c < bhm_target_cartons.html


# Then do this auto conversion to latex:
#
TEXOUT=bhm_target_cartons.tex
pandoc --from=html --to=latex --output=$TEXOUT bhm_target_cartons.html

# replace the references with bibtex equivalents:

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
awk 'BEGIN {flag=1} NF==0 {flag=1} $0~/Catalogdb tables required/ {flag=0} flag==1 {print $0}' $TEXOUT >  temp.$TEXOUT
mv temp.$TEXOUT  $TEXOUT

# checks (should result in zero rows)
grep -A1 adsab $TEXOUT

