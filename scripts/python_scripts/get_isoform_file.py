import pandas as pd
import numpy as np

columns = [
    "sequence",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

gtf_df = pd.read_csv('braker.gtf', sep='\t', names=columns)
name = gtf_df[gtf_df['feature']=='transcript']['attribute']
gene = name.str.split('.').str.get(0)
transcript = name.str.split('.').str.get(1)
df = pd.DataFrame(data={'gene':gene, 'name':name})
df.to_csv('afusca.braker3.isoforms.tsv',sep='\t', header=False,index=False)
#pivot = df.groupby('gene', sort=False).name.apply(list).apply(pd.Series).add_prefix('t').reset_index()
#pivot.to_csv('afusca.braker3.isoforms.tsv',sep='\t', header=False,index=False)