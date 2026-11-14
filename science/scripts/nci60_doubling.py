# Parse NCI-60 doubling times from Table 1 of Marshall et al. 2017 (PMC5654365), JATS XML.
import re, sys, csv, html
x = open(sys.argv[1], encoding='utf-8').read()
tabs = re.findall(r'<table-wrap\b.*?</table-wrap>', x, flags=re.S)
t = next(t for t in tabs if 'COLO 205' in t and 'Epithelial' in t)
cap = re.sub(r'<[^>]+>', ' ', re.search(r'<caption>(.*?)</caption>', t, re.S).group(1))
print('CAPTION:', ' '.join(cap.split()), file=sys.stderr)
rows = []
for tr in re.findall(r'<tr\b.*?</tr>', t, flags=re.S):
    cells = [' '.join(html.unescape(re.sub(r'<[^>]+>', ' ', c)).split()) for c in re.findall(r'<t[dh]\b[^>]*>(.*?)</t[dh]>', tr, flags=re.S)]
    rows.append(cells)
for r in rows[:4]: print(r, file=sys.stderr)
w = csv.writer(sys.stdout)
w.writerow(['cell_line', 'panel', 'doubling_time_h', 'source'])
panel = ''
for r in rows:
    nums = [c for c in r if re.fullmatch(r'\d+(\.\d+)?', c)]
    if len(r) >= 4 and nums:
        if r[0]: panel = r[0]
        name = r[1]
        # first numeric cell is doubling time (column order: panel, line, histology, type, doubling, year...)
        w.writerow([name, panel, nums[0], 'NCI-60 (DTP), via Marshall et al. 2017 Table 1, PMC5654365'])
