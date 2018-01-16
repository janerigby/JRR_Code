import pandas
df = pandas.DataFrame({'A' : [5,6,3,4], 'B' : [1,2,3,5]})
mylist = [3,6]
ndf =  df[df['A'].isin(mylist)]
ndf['sort_cat'] = pandas.Categorical(ndf['A'], categories=mylist, ordered=True)
ndf.sort_values('sort_cat', inplace=True)
ndf.reset_index(inplace=True)
print ndf
