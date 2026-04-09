#!/usr/bin/env sh

# This generates the HTML for the tables of cartons on the following page:
# https://testng.sdss.org/${DR}/bhm/programs/cartons/

# run some queries against the database, format the output, paste the result into 'Custom HTML' clock in wordpress.

DR=dr20
BASEDIR=~/SDSSV/gitwork/target_selection/docs/website/bhm_cartons/$DR
mkdir -p $BASEDIR

cd $BASEDIR

# MINIDB=minidb
MINIDB=minidb_$DR

# sdss5db=> select * from minidb_dr20.dr20_targeting_generation;
#  pk |         label         | first_release
# ----+-----------------------+---------------
#   1 | v0.plates             | dr19
#   2 | v0.5.epsilon-7-core-0 | dr19
#   3 | v0.5.2                | dr19
#   4 | v0.5.3                | dr18
#   5 | v0.5.5                | dr19
#   6 | v1.0.2                | dr20
#   7 | v1.0.3                | dr20
#   8 | v1.0.4                | dr20
#   9 | v1.0.5                | dr20
#  13 | v1.0.5-boss-only      | dr20
#   (10 rows)

GEN_LIST="v0.plates v0.5.epsilon-7-core-0 v0.5.2 v0.5.3 v0.5.5 v1.0.2 v1.0.3 v1.0.4 v1.0.5 v1.0.5-boss-only dr20.manual"
#GEN_LIST="v0.plates"
#GEN_LIST="v1.0.2 v1.0.3 v1.0.4 v1.0.5 v1.0.5-boss-only"

# change the below to your favoured output location
OUTDIR=~/SDSSV/gitwork/target_selection/docs/website/bhm_cartons/${DR}
mkdir -p $OUTDIR
cd $OUTDIR

# this one when running from home
alias sdssdb='psql -h localhost -p 7502 -U sdss_user -d sdss5db'

# this one when running on a Utah machine
#alias sdssdb='psql -h operations.sdss.org -U sdss_user -d sdss5db'

#use this in preference# Q="\copy ( SELECT 
#use this in preference#           c.carton,v.plan,v.tag,
#use this in preference#           CASE WHEN (c.carton ~ 'bhm_colr' OR c.carton ~ 'bhm_gua')  THEN 'ancillary'
#use this in preference#                WHEN (c.program = 'bhm_filler' OR c.carton ~ 'bhm_csc')  THEN 'non-core'
#use this in preference#                WHEN (c.program = 'bhm_spiders' AND (c.carton ~ '_supercosmos' OR c.carton ~ '_efeds_stragglers' OR c.carton ~ '_skymapperdr2' OR c.carton ~ 'agn_gaia' OR c.carton ~ '_sep'))  THEN 'non-core'
#use this in preference#                WHEN (c.program = 'open_fiber')  THEN 'openfiber'
#use this in preference#                ELSE 'core' END as status,
#use this in preference#           CASE WHEN (c.program ~ 'gua' OR c.program ~ 'colr' OR c.program ~ 'open_fiber' OR c.program ~ 'manual') THEN 'ancillary' ELSE c.program END as program, 
#use this in preference#           tg.label
#use this in preference# FROM ${MINIDB}.${DR}_targeting_generation as tg
#use this in preference# JOIN ${MINIDB}.${DR}_targeting_generation_to_carton as tg2c
#use this in preference# ON tg.pk = tg2c.generation_pk
#use this in preference# join ${MINIDB}.${DR}_carton as c
#use this in preference# on tg2c.carton_pk = c.carton_pk 
#use this in preference# join  ${MINIDB}.${DR}_targetdb_version as v
#use this in preference#   on c.version_pk = v.pk
#use this in preference# where ( (c.carton ~ 'bhm' 
#use this in preference#          OR c.carton IN ('openfibertargets_nov2020_11', 'openfibertargets_nov2020_18', 'openfibertargets_nov2020_26', 'openfibertargets_nov2020_27', 'openfibertargets_nov2020_30', 'openfibertargets_nov2020_33') ) 
#use this in preference#   AND v.plan = c.target_selection_plan )
#use this in preference# order by c.carton_pk)
#use this in preference# to 'q_result.csv' with csv header"

# version of the above that uses the full targetdb version of the targeting generation in order to get around the missing 'dr20.manual' entry in minidb_dr20.dr20_targeting_generation
Q="\copy ( SELECT 
          c.carton,v.plan,v.tag,
          CASE WHEN (c.carton ~ 'bhm_colr' OR c.carton ~ 'bhm_gua')  THEN 'ancillary'
               WHEN (c.program = 'bhm_filler' OR c.carton ~ 'bhm_csc')  THEN 'non-core'
               WHEN (c.program = 'bhm_spiders' AND (c.carton ~ '_supercosmos' OR c.carton ~ '_efeds_stragglers' OR c.carton ~ '_skymapperdr2' OR c.carton ~ 'agn_gaia' OR c.carton ~ '_sep'))  THEN 'non-core'
               WHEN (c.program = 'open_fiber')  THEN 'openfiber'
               ELSE 'core' END as status,
          CASE WHEN (c.program ~ 'gua' OR c.program ~ 'colr' OR c.program ~ 'open_fiber' OR c.program ~ 'manual') THEN 'ancillary' ELSE c.program END as program, 
          tg.label
FROM targeting_generation as tg
JOIN targeting_generation_to_carton as tg2c
ON tg.pk = tg2c.generation_pk
join ${MINIDB}.${DR}_carton as c
on tg2c.carton_pk = c.carton_pk 
join targetdb.version as v
  on c.version_pk = v.pk
WHERE ( tg.first_release = '"$DR"' 
        AND (c.carton ~ 'bhm' 
         OR c.carton IN ('openfibertargets_nov2020_11', 'openfibertargets_nov2020_18', 'openfibertargets_nov2020_26', 'openfibertargets_nov2020_27', 'openfibertargets_nov2020_30', 'openfibertargets_nov2020_33') ) 
  AND v.plan = c.target_selection_plan )
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
$1~/bhm_/ && $6 == gen {
    prog=$5;  gsub("bhm_", "", prog); PROG=sprintf("BHM %s",toupper(prog)); tag=$3;plan=$2;
    if (prog=="filler") {prog="ancillary";};
    if (prog=="ancillary") {PROG="BHM Ancillary programs";};
    printf("<tr>\
<td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/%s/bhm/programs/%s\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/blob/%s/python/target_selection/config/target_selection.yml\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
</tr>\n",
    $1, plan, $1, dr, prog, PROG, tag, plan, tag, tag, $4)}
END {
printf("</tbody></table><figcaption>All BHM cartons from targeting generation \"%s\"</figcaption></figure>\n", gen);
}' q_result.csv > carton_table_block_generation_${GEN}.html

done



# Goup similar targeting generations in same table
# membership of targeting generations is done as set of ticks in bool columns
# when the minidb version is fully populated replace:
#
# FROM targeting_generation as tg
# JOIN targeting_generation_to_carton as tg2c
#
# with
#
# FROM ${MINIDB}.${DR}_targeting_generation as tg
# JOIN ${MINIDB}.${DR}_targeting_generation_to_carton as tg2c
#

Q="\copy ( SELECT 
          c.carton,v.plan,v.tag,
          CASE WHEN (c.carton ~ 'bhm_colr' OR c.carton ~ 'bhm_gua')  THEN 'ancillary'
               WHEN (c.program = 'bhm_filler' OR c.carton ~ 'bhm_csc')  THEN 'non-core'
               WHEN (c.program = 'bhm_spiders' AND (c.carton ~ '_supercosmos' OR c.carton ~ '_efeds_stragglers' OR c.carton ~ '_skymapperdr2' OR c.carton ~ 'agn_gaia' OR c.carton ~ '_sep'))  THEN 'non-core'
               WHEN (c.program = 'open_fiber')  THEN 'openfiber'
               ELSE 'core' END as status,
          CASE WHEN (c.program ~ 'gua' OR c.program ~ 'colr' OR c.program ~ 'open_fiber' OR c.program ~ 'manual') THEN 'ancillary' ELSE c.program END as program, 
          SUM(CASE WHEN tg.label = 'v0.5.epsilon-7-core-0' THEN 1 ELSE 0 END) AS in_tg1,
          SUM(CASE WHEN tg.label = 'v0.5.2' THEN 1 ELSE 0 END) AS in_tg2,
          SUM(CASE WHEN tg.label = 'v0.5.3' THEN 1 ELSE 0 END) AS in_tg3,
          SUM(CASE WHEN tg.label = 'v0.5.5' THEN 1 ELSE 0 END) AS in_tg4,
          SUM(CASE WHEN (tg.label = 'v1.0.2' or tg.label = 'v1.0.3' or tg.label = 'v1.0.4') THEN 1 ELSE 0 END) AS in_tg5,
          SUM(CASE WHEN (tg.label = 'v1.0.5' or tg.label = 'v1.0.5-boss-only') THEN 1 ELSE 0 END) AS in_tg6,
          SUM(CASE WHEN (tg.label = 'dr20.manual') THEN 1 ELSE 0 END) AS in_tg7,
          ARRAY_AGG(distinct tg.label) as targeting_generations
FROM targeting_generation as tg
JOIN targeting_generation_to_carton as tg2c
ON tg.pk = tg2c.generation_pk
join ${MINIDB}.${DR}_carton as c
on tg2c.carton_pk = c.carton_pk 
join  ${MINIDB}.${DR}_targetdb_version as v
  on c.version_pk = v.pk
where ( tg.first_release = '"$DR"'
        AND (c.carton ~ 'bhm' 
         OR c.carton IN ('openfibertargets_nov2020_11', 'openfibertargets_nov2020_18', 'openfibertargets_nov2020_26', 'openfibertargets_nov2020_27', 'openfibertargets_nov2020_30', 'openfibertargets_nov2020_33') ) 
  AND v.plan = c.target_selection_plan )
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


#too wide# gawk --field-separator=',' \
#too wide#          -v dr=$DR \
#too wide#          -v gens="v0.5.epsilon-7-core-0 and v0.5.2/v0.5.3/v0.5.5" \
#too wide# 'BEGIN {
#too wide#     split(gens,agens," and ");
#too wide#     printf("<figure class=\"wp-block-table is-style-stripes\">\n\
#too wide# <table>\n\
#too wide# <thead>\n\
#too wide# <tr>\
#too wide# <th>carton name</th>\
#too wide# <th class=\"has-text-align-center\" data-align=\"center\">science <br>program</th>\
#too wide# <th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>plan</th>\
#too wide# <th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>tag</th>\
#too wide# <th class=\"has-text-align-center\" data-align=\"center\">status</th>\
#too wide# <th class=\"has-text-align-center\" data-align=\"center\">%s</th>\
#too wide# <th class=\"has-text-align-center\" data-align=\"center\">%s</th>\
#too wide# </tr>\
#too wide# \n</thead>\n<tbody>\n",
#too wide# agens[1], agens[2]);}
#too wide# $1~/^bhm_/ && $NF~/v0.5/ {
#too wide#     prog=$5;  gsub("bhm_", "", prog); PROG=sprintf("BHM %s",toupper(prog)); tag=$3;plan=$2;
#too wide#     if (prog=="filler") {prog="ancillary";};
#too wide#     if (prog=="ancillary") {PROG="BHM Ancillary programs";};
#too wide#     printf("<tr>\
#too wide# <td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a></td>\
#too wide# <td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/%s/bhm/programs/%s\">%s</a></td>\
#too wide# <td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/blob/%s/python/target_selection/config/target_selection.yml\">%s</a></td>\
#too wide# <td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td>\
#too wide# <td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
#too wide# <td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
#too wide# <td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
#too wide# </tr>\n",
#too wide#     $1, plan, $1, dr, prog, PROG, tag, plan, tag, tag, $4, 
#too wide#     ($6 > 0 ? "&#9989;" : ""), 
#too wide#     ($7 > 0 || $8 > 0 || $9 > 0 ? "&#9989;" : "") )}
#too wide# END {
#too wide# printf("</tbody></table><figcaption>All BHM cartons from targeting generations %s</figcaption></figure>\n", gens);
#too wide# }' q2_result.csv > carton_table_block_multi_generation_v0.5.html


# try to make the table narrower so it can fit without a scrollbar
# compactify the information in the targeting_generation membership columns
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
</tr>\
\n</thead>\n<tbody>\n");}
$1~/^bhm_|^openfiber/ && $NF~/v0.5/ {
    prog=$5;  gsub("bhm_", "", prog); PROG=sprintf("BHM %s",toupper(prog)); tag=$3;plan=$2;
    if (prog=="filler") {prog="ancillary";};
    if (prog=="ancillary") {PROG="BHM Ancillary programs";};
    if (($6 > 0) && ($7 > 0 || $8 > 0 || $9 > 0)) {tgsym=""} 
    else if ($6>0) {tgsym=sprintf("<span title=\"only in targeting generation %s\">&Dagger;</span>",agens[1])} 
    else {tgsym=sprintf("<span title=\"only in targeting generations %s\">&clubs;</span>",agens[2])}; 
    printf("<tr>\
<td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a>%s</td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/%s/bhm/programs/%s\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/blob/%s/python/target_selection/config/target_selection.yml\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
</tr>\n",
    $1, plan, $1, tgsym, dr, prog, PROG, tag, plan, tag, tag, $4)}
END {
printf("</tbody></table><figcaption>All BHM cartons (with versions) from targeting generations %s. Carton-versions marked with a &Dagger; symbol are only present in targeting generation %s. Carton-versions marked with a &clubs; symbol are only present in targeting generations %s. </figcaption></figure>\n", gens, agens[1], agens[2]);
}' q2_result.csv > carton_table_block_multi_generation_v0.5_slim.html


# v1 version of above
gawk --field-separator=',' \
         -v dr=$DR \
         -v gens="v1.0.2/v1.0.3/v1.0.4 and v1.0.5/v1.0.5-boss-only" \
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
</tr>\
\n</thead>\n<tbody>\n");}
$1~/^bhm_|^openfiber/ && $NF~/v1.0/ {
    prog=$5;  gsub("bhm_", "", prog); PROG=sprintf("BHM %s",toupper(prog)); tag=$3;plan=$2;
    if (prog=="filler") {prog="ancillary";};
    if (prog=="ancillary") {PROG="BHM Ancillary programs";};
    if (($10>0) && ($11>0)) {tgsym=""} 
    else if ($10>0) {tgsym=sprintf("<span title=\"only in targeting generation %s\">&Dagger;</span>",agens[1])} 
    else {tgsym=sprintf("<span title=\"only in targeting generations %s\">&clubs;</span>",agens[2])}; 
    printf("<tr>\
<td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a>%s</td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/%s/bhm/programs/%s\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/blob/%s/python/target_selection/config/target_selection.yml\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td>\
<td class=\"has-text-align-center\" data-align=\"center\">%s</td>\
</tr>\n",
    $1, plan, $1, tgsym, dr, prog, PROG, tag, plan, tag, tag, $4)}
END {
printf("</tbody></table><figcaption>All BHM cartons (with versions) from targeting generations %s. Carton-versions marked with a &clubs; symbol are only present in targeting generations %s. </figcaption></figure>\n", gens, agens[2]);
}' q2_result.csv > carton_table_block_multi_generation_v1.0_slim.html


# paste into HTML box in wordpress - these comands help
xclip -sel c < carton_table_block_generation_v0.plates.html
xclip -sel c < carton_table_block_multi_generation_v0.5_slim.html
xclip -sel c < carton_table_block_multi_generation_v1.0_slim.html
xclip -sel c < carton_table_block_generation_dr20.manual.html


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

# build empty json structures with all easy info filled in 

# first collate the automatically extractable info from targetdb

Q="\copy ( SELECT 
          c.carton_pk,c.carton,v.plan,v.tag,
          array_agg(DISTINCT tg.label) as targeting_generations,
          count(DISTINCT c2t.target_pk) as ntargets,
          array_agg(DISTINCT cad.label) as cadences,
          array_agg(DISTINCT c2t.priority) as priority
FROM ${MINIDB}.${DR}_targeting_generation as tg
JOIN ${MINIDB}.${DR}_targeting_generation_to_carton as tg2c
ON tg.pk = tg2c.generation_pk
JOIN ${MINIDB}.${DR}_carton as c
on tg2c.carton_pk = c.carton_pk 
JOIN  ${MINIDB}.${DR}_targetdb_version as v
  on c.version_pk = v.pk
JOIN  ${MINIDB}.${DR}_carton_to_target as c2t
  ON c.carton_pk = c2t.carton_pk
JOIN  ${MINIDB}.${DR}_cadence as cad
  ON c2t.cadence_pk = cad.pk
WHERE ( (c.carton ~ 'bhm' 
         OR c.carton IN ('openfibertargets_nov2020_11', 'openfibertargets_nov2020_18', 'openfibertargets_nov2020_26', 'openfibertargets_nov2020_27', 'openfibertargets_nov2020_30', 'openfibertargets_nov2020_33') ) 
  AND v.plan = c.target_selection_plan )
GROUP BY c.carton_pk,c.carton,v.plan,v.tag
ORDER BY c.carton_pk)
to 'q3_result.csv' with DELIMITER '|' csv header"

sdssdb -c "$Q"


# now convert (post v1 cartons) into a json framework
gawk --field-separator='|' '$1!="carton_pk" && ($3 >= "1.0.0" || $2!~/openfiber/ ) {
crossmatch=($3>="1.0.0" ? "1.0.0" : ($3>="0.5.0" ? "0.5.0" : "0.1.0"));
printf("{\n    \"name\": \"%s\",\n    \"plan\": \"%s\",\n    \"tag\": \"%s\",\n    \"summary\": \"\",\n    \"selection\": \"\",\n    \"tables\": \"\",\n    \"crossmatch\": \"%s\",\n    \"cadences\": \"%s\",\n    \"priority\": \"%s\",\n    \"code\": \"\",\n    \"ntargets\": \"%d\"\n},\n", 
$2, $3, $4, crossmatch, $7, $8, $6)}' q3_result.csv

###################################################
###################################################
###################################################
###################################################
#Now edit JSON (carton_descriptions.json) manually
###################################################
###################################################
###################################################
###################################################


### now process the json
# see: proc_carton_desc.py
# e.g.
python proc_carton_desc_${DR}.py > bhm_target_cartons.html

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
}' bhm_target_cartons.html carton_table_block_generation_v0.plates.html carton_table_block_multi_generation_v0.5_slim.html carton_table_block_multi_generation_v1.0_slim.html carton_table_block_generation_dr20.manual.html


# paste carton descriptions into HTML box in wordpress - this command helps
xclip -sel c < bhm_target_cartons.html


##Now convert to latex

./convert_html_to_latex.sh
