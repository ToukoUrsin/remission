import pandas as pd, sys
d = pd.read_excel(sys.argv[1], sheet_name=None)
for k, v in d.items():
    print(k, v.shape); print(v.columns.tolist()); print(v.head(3).to_string())
