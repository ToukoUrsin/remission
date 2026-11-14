# Strip XML tags from a downloaded full-text file and print +-W chars around regex matches (data only).
import re, sys, html
t = open(sys.argv[1], encoding='utf-8', errors='replace').read()
t = html.unescape(re.sub(r'<[^>]+>', ' ', t)); t = ' '.join(t.split())
n = int(sys.argv[3]) if len(sys.argv) > 3 else 30
w = int(sys.argv[4]) if len(sys.argv) > 4 else 200
last = -1
for i, m in enumerate(re.finditer(sys.argv[2], t, flags=re.I)):
    if i >= n: break
    if m.start() < last: continue
    print('-', t[max(0, m.start() - w): m.end() + w]); last = m.end() + w
