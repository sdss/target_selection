#!/usr/bin/env python
import json

descfile = "carton_descriptions.json"
outfile = "bhm_cartons_table.tex"
url_ts = "https://github.com/sdss/target_selection"


with open(descfile) as f:
    j = json.load(f)


# write table preamble

s = '''
% --------------- BHM Cartons ------------------------\%
\\begin{table}[ht]
\\centering
\\caption{New BHM cartons in DR20.
The tag and plan columns are used in versioning the code used to generate the carton.
$N_\\mathrm{targets}$ gives the number of targets in the carton.}
\\label{tab:boss_cartons}
\\begin{tabular}{lccr}
\\hline
carton name & tag & plan & $N_\\mathrm{targets}$ \\\\
\\hline'''


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
            "\\texttt{" + c["name"].replace("_", "\\_") + "}",
            c["plan"],
            c["tag"],
            c["ntargets"],
        ]
        str_out = " & ".join(s)
        print(str_out, "\\\\")
        of.write(str_out)

    s = '''\\hline
\\end{tabular}
\\end{table}
% --------------- BHM Cartons ------------------------%
'''

    print(s)
    of.write(s)

    of.close()
