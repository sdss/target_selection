#!/bin/bash

# can re-run the json -> html conversion again to make sure that we are working from most up to date .json

DR=dr20

### now process the json
# see: proc_carton_desc.py
# e.g.
python proc_carton_desc_${DR}.py > bhm_target_cartons.html



# Then do this auto conversion to latex:
#
TEXOUT=bhm_target_cartons.tex
pandoc --from=html --to=latex --output=$TEXOUT bhm_target_cartons.html

# replace the references with bibtex equivalents:

perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2020ApJS..250....8L\/abstract}{Lyke\net al., 2020}/\\citealt{Lyke2020}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2014ApJ...785..104R\/abstract}{Rykoff\net al., 2014}/\\citealt{Rykoff2014}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2020MNRAS.499.4768I\/abstract}{Ider\nChitham et al., 2020}/\\citealt{IderChitham2020}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022A\\%26A...661A...3S\/abstract}{Salvato\net al., 2022}/\\citealt{Salvato2022}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2019MNRAS.489.4741S\/abstract}{Shu\net al., 2019}/\\citealt{Shu2019}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2019MNRAS.489.4741S\/abstract}{Shu\net al. \(2019\)}/\\citet{Shu2019}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2021ApJS..253....8M\/abstract}{CatWISE2020}/CatWISE2020 \\citep{Marocco2021}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022arXiv220608989Y\/abstract}{Yang\nand Shen \(2022\)}/\\citet{Yang2022}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022arXiv220608989Y\/abstract}{Yang\nand Shen, 2022}/\\citealt{Yang2022}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2019PASA...36...33O\/abstract}{SkyMapper-dr2}/SkyMapper-dr2 \\citep{Onken2019}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022A\\%26A...661A...2L\/abstract}{Liu\net al., 2022}/\\citealt{Liu2022}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022ApJS..259...35A\/abstract}{Abdurro'\''uf\net al., 2022}/\\citealt{Abdurrouf_2021_sdssDR17}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2022ApJS..259...35A\/abstract}{Abdurro\\textquotesingle uf\net al., 2022}/\\citealt{Abdurrouf_2021_sdssDR17}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2011ApJ...729..141B\/abstract}{Bovy\net al., 2011}/\\citealt{Bovy2011}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2023A\\%26A...678A.151H\/abstract}{Hardcastle\net al., 2023}/\\citealt{Hardcastle2023}/g'  $TEXOUT

perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2020zndo...4279451S\/abstract}{Sanchez-Saez\net al., 2020}/\\citealt{SanchezSaez2020}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2023A\\%26A...675A.195S\/abstract}{Sanchez-Saez\net al. \(2023\)}/\\citealt{SanchezSaez2023}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2009ApJ...692..758G\/abstract}{Gibson\net al., 2009}/\\citealt{Gibson2009}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2021PASA...38...58H\/abstract}{Hale\net al., 2021}/\\citealt{Hale2021}/g'  $TEXOUT

perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2023MNRAS.521.1620D\/abstract}{Delaney\net al., 2023}/\\citealt{Delaney2023}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2020A\\%26A...641A.136W\/abstract}{Webb\net al., 2020}/\\citealt{Webb2020}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2018ApJ...868...99F\/abstract}{French\n\\& Zabludoff \(2018\)}/\\citealt{French2018}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2021ApJS..256...21N\/abstract}{Ni\net al. \(2021\)}/\\citealt{Ni2021}/g'  $TEXOUT
perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/2026A\\%26A...706A.144W\/abstract}{Waddell\net al. \(2026\)}/\\citealt{Waddell2026}/g'  $TEXOUT

# perl -0777 -pi -e 's/\\href{https:\/\/ui.adsabs.harvard.edu\/abs\/xxx\/abstract}{xxx\net al., xxxx}/\\citealt{xxxxx}/g'  $TEXOUT
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


# now select only cartons that were not previously reported in DR18, DR19

gawk '/plan1.0/ {flag=1} flag==1 {print $0}' $TEXOUT > ${TEXOUT%.*}_${DR}.tex

echo "Now upload ${TEXOUT%.*}_${DR}.tex into the DR20 overleaf
"

echo "Don't forget to paste carton descriptions into HTML box in wordpress - this command helps:
xclip -sel c < bhm_target_cartons.html
"
