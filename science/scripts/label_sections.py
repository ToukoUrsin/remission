# Print selected sections of an openFDA label JSON (data only).
import json, sys, re
d = json.load(open(sys.argv[1]))['results'][0]
keys = sys.argv[2].split(',')
print('##', d.get('openfda', {}).get('brand_name'), d.get('effective_time'), d.get('set_id'), 'version', d.get('version'))
for k in keys:
    for kk in d:
        if kk == k:
            for txt in d[kk]:
                print(f'--- {kk}:'); print(re.sub(r'\s+', ' ', txt)[:int(sys.argv[3]) if len(sys.argv) > 3 else 20000])
