#!/usr/bin/env python
import json

descfile = "carton_descriptions.json"
outfile = "bhm_cartons_table.tex"
url_ts = "https://github.com/sdss/target_selection"


with open(descfile) as f:
    j = json.load(f)


# write table preamble

s = '''
% --------------- BHM Cartons ------------------------%
\\begin{table}[ht]
\\centering
\\caption{New BHM (and extragalactic openfiber) cartons in DR20.
The plan column specifies the version of \\texttt{target\\_selection} code used to generate each
carton\\footnote{\\url{https://github.com/sdss/target_selection}{github.com/sdss/target\\_selection}}
$N_\\mathrm{targets}$ gives the number of targets in each carton.}
\\label{tab:bhm_cartons}
%\\begin{tabular}{lccr}
\\begin{tabular}{lcr}
%Carton name & plan & tag & $N_\\mathrm{targets}$ \\\\
Carton name & plan & $N_\\mathrm{targets}$ \\\\
\\hline
'''


with open(outfile, "wt") as of:
    print(s)
    of.write(s)

    for c in j["cartons"]:
        if c["name"] == "":
            continue
        if c["plan"] < "1.0.0":
            continue
        name = c["name"]
        plan = c.get("plan", "")
        tag = c.get("tag", "")
        crossmatch = c["crossmatch"]

        s = [
            '\\hyperref[' + c["name"] + '_plan' + c['plan'] + ']{\\texttt{' + c["name"].replace('_', '\\_') + '}}',
            c['plan'],
            # c['tag'],
            c['ntargets'],
        ]
        str_out = " & ".join(s) + " \\\\"
        print(str_out)
        of.write(str_out + '\n')

    s = '''
\\hline
\\end{tabular}
\\end{table}
% --------------- BHM Cartons ------------------------%
'''

    print(s)
    of.write(s)

    of.close()
