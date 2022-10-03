#!/usr/bin/env sh

# run some queries against the database, format the output, paste the result into 'Custom HTML' clock in wordpress.

# This gnerates the HTML for the tables of cartons on the following page:
# https://testng.sdss.org/dr18/bhm/programs/cartons/

#alias sdssdb='psql -h localhost -p 7502 -U sdss_user -d sdss5db'
alias sdssdb='psql -h operations.sdss.org -U sdss_user -d sdss5db'

Q="SELECT c.carton,c.target_selection_plan,v.tag,
          CASE WHEN cc.program = 'bhm_filler' THEN 'non-core' ELSE 'core' END as status
FROM minidb.dr18_carton as c
join carton as cc
  on c.carton = cc.carton
join targetdb.version as v
  on cc.version_pk = v.pk
where c.carton ~ 'bhm'
  and v.plan = c.target_selection_plan
order by c.carton;"

sdssdb -c "$Q" > q_result.txt


gawk -v gen='"v0.5.3"' 'BEGIN {printf("<figure class=\"wp-block-table is-style-stripes\">\n<table>\n<thead>\n<tr><th>carton name</th><th class=\"has-text-align-center\" data-align=\"center\">science <br>program</th><th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>plan</th><th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>tag</th><th class=\"has-text-align-center\" data-align=\"center\">status</th><th>notes</th></tr>\n</thead>\n<tbody>\n")}
$1~/^bhm_/{
split($1,c,"_"); prog=c[2]; PROG=toupper(prog); if (prog~/gua|colr/) {prog="ancillary";PROG="Ancillary programs";};
printf("<tr><td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a></td><td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/dr18/bhm/programs/%s\">BHM %s</a></td><td class=\"has-text-align-center\" data-align=\"center\">%s</td><td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td><td class=\"has-text-align-center\" data-align=\"center\">%s</td><td></td></tr>\n",
$1, $3, $1, prog, PROG, $3, $5, $5, $7)}
END {printf("</tbody></table><figcaption>A table of all BHM cartons from targeting generation %s</figcaption></figure>\n", gen)}' q_result.txt > carton_table_block_generation_0.5.3.html




######################################
## manually select the eFEDS plates cartons

Q="SELECT cc.carton,v.plan,v.tag,
          'non-core' as status
FROM carton as cc
join targetdb.version as v
  on cc.version_pk = v.pk
where
  cc.carton ~ '-efeds'
  and v.plan ~ '0.1.0'
order by cc.carton;"

sdssdb -c "$Q" > q_result.txt

gawk -v gen='"eFEDS plates"' 'BEGIN {printf("<figure class=\"wp-block-table is-style-stripes\">\n<table>\n<thead>\n<tr><th>carton name</th><th class=\"has-text-align-center\" data-align=\"center\">science <br>program</th><th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>plan</th><th class=\"has-text-align-center\" data-align=\"center\">target <br>selection <br>tag</th><th class=\"has-text-align-center\" data-align=\"center\">status</th><th>notes</th></tr>\n</thead>\n<tbody>\n")}
$1~/^bhm_/{
split($1,c,"_"); prog=c[2]; PROG=toupper(prog); if (prog~/gua|colr/) {prog="ancillary";PROG="Ancillary programs";};
printf("<tr><td><a href=\"#%s_plan%s\" data-type=\"internal\">%s</a></td><td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://testng.sdss.org/dr18/bhm/programs/%s\">BHM %s</a></td><td class=\"has-text-align-center\" data-align=\"center\">%s</td><td class=\"has-text-align-center\" data-align=\"center\"><a href=\"https://github.com/sdss/target_selection/tree/%s/\">%s</a></td><td class=\"has-text-align-center\" data-align=\"center\">%s</td><td></td></tr>\n",
$1, $3, $1, prog, PROG, $3, $5, $5, $7)}
END {printf("</tbody></table><figcaption>A table of all BHM cartons from targeting generation %s</figcaption></figure>\n", gen)}' q_result.txt > carton_table_block_generation_eFEDS_plates.html



       
### now process the json
# see: proc_carton_desc.py
