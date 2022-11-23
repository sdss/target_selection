#!/usr/bin/env python
import json

'''
This is a routine to generate BHM target carton descriptions in HTML
using content from a json input file. The output can then be
manually pasted into a wordpress 'Custom HTML' block.

Specifically this generates the "BHM Carton descriptions" section of:
https://testng.sdss.org/dr18/bhm/programs/cartons/

Objects in JSON file must have following items:
  "name", "plan", "tag", "summary",
  "selection", "tables", "cadences",
  "priority", "code", "ntargets"

'''

descfile = "carton_descriptions.json"
outfilehtml = "carton_descriptions.html"
url_ts = "https://github.com/sdss/target_selection"
div_class = "wp-container-2 wp-block-group has-cyan-bluish-gray-background-color has-background"
div_style = "padding-top:3%;padding-right:3%;padding-bottom:3%;padding-left:3%"

with open(descfile) as f:
    j = json.load(f)

with open(outfilehtml, "wt") as of:

    for c in j['cartons']:
        if c['name'] == '':
            continue
        plan = c['plan']
        tag = c['tag']
        name = c['name']
        this_url_ts = f'{url_ts}/blob/{tag}/python/target_selection/cartons/{c["code"]}'
        s = [
            '\n<hr class="wp-block-separator has-alpha-channel-opacity">\n',
            f'<div class="{div_class} style="{div_style}">'
            f'<div class="wp-block-group__inner-container">',
            f'<h3 id="{name}_plan{plan}">{name}</h3>',
            f'<p><strong>target_selection plan:</strong> {plan}</p>',
            f'<p><strong>target_selection tag:</strong> '
            f'<a href="{url_ts}/tree/{tag}/">{tag}</a></p>',
            f'<p><strong>Summary:</strong> {c["summary"]}</p>',
            f'<p><strong>Simplified description of selection criteria:'
            f'</strong> {c["selection"]}</p>',
            f'<p><strong>Catalogdb tables required:</strong> {c["tables"]}</p>',
            f'<p><strong>Target priority options:</strong> {c["priority"]}</p>',
            f'<p><strong>Cadence options:</strong> {c["cadences"]}</p>',
            f'<p><strong>Implementation:</strong> '
            f'<a href="{this_url_ts}">{c["code"]}</a></p>',
            f'<p><strong>Number of targets:</strong> {c["ntargets"]}</p>',
            '</div></div>',
        ]
        str_out = "\n".join(s)
        print(str_out)
        of.write(str_out)

    of.close()
