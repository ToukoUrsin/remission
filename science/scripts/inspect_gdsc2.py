import pandas as pd, sys
d = pd.read_excel(sys.argv[1])
d.to_parquet if False else None
print(d.shape); print(d.columns.tolist()); print(d.head(3).to_string())
names = sorted(d.DRUG_NAME.unique())
kw = ['dabraf','plx','vemur','encoraf','trametin','selumet','binimet','gefit','erlot','afat','cetux','alpel','pictil','taselis','navito','venet','wehi','a-1155463','a1155463','a-1331852','a-1331','obatoclax','idasanutlin','abt','lapatinib','osimert','buparl','gdc','byl']
for n in names:
    if any(k in n.lower() for k in kw): 
        s = d[d.DRUG_NAME==n]
        print(repr(n), s.DRUG_ID.unique().tolist(), s.PUTATIVE_TARGET.iloc[0], len(s), s.MIN_CONC.min(), s.MAX_CONC.max())
d.to_csv(sys.argv[2], index=False)
