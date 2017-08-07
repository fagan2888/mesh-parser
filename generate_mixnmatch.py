# https://tools.wmflabs.org/mix-n-match/import.php
# https://meta.wikimedia.org/wiki/Mix%27n%27match/Manual

import json

mesh_tree = {'A': 'Anatomy',
             'B': 'Organisms',
             'C': 'Diseases',
             'D': 'Chemicals and Drugs',
             'E': 'Analytical, Diagnostic and Therapeutic Techniques and Equipment',
             'F': 'Psychiatry and Psychology',
             'G': 'Phenomena and Processes',
             'H': 'Disciplines and Occupations',
             'I': 'Anthropology, Education, Sociology and Social Phenomena',
             'J': 'Technology, Industry, Agriculture',
             'K': 'Humanities',
             'L': 'Information Science',
             'M': 'Named Groups',
             'N': 'Health Care',
             'V': 'Publication Characteristics',
             'Z': 'Geographicals'}

d = json.load(open("mesh.json"))

for tree_start, tree_name in mesh_tree.items():
    s = {k:v for k,v in d.items() if 'tree' in v and any(tree.startswith(tree_start) for tree in v['tree'])}
    with open(tree_name + ".tsv", 'w') as f:
        print('\n'.join(['\t'.join([k,x['term'],x['note'] if 'note' in x else '']) for k,x in s.items()]), file=f)
