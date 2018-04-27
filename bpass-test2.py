import pandas

def bpass_model_dir() : # Set this to where your models are
    return('/Volumes/Apps_and_Docs/JRR_Utils/BPASS_v2.1/')

def bpass_ages_setup() :  # Create a pandas dataframe of ages and column names ('coln'), as in BPASS v2.1 manual
    df_age = pandas.DataFrame([ 10**(6+0.1*(n-2)) for n in range(2,52+1)], columns=('age',))
    colnames = pandas.Series(["col"+str(n) for n in range(2,52+1)])
    df_age['colname'] = colnames
    return(df_age)

def find_closest_bpass_age(df_age, age_to_find):
    closest_age = df_age.iloc[(df_age['age']- age_to_find).abs().argsort()[:1]]
    return(closest_age) # Returns the row of df_age thats closest in age.  

def load_bpass_spectra(filename) :
    bpassdir = bpass_model_dir()
    df_age = bpass_ages_setup()
    df_bpass = pandas.read_table(bpassdir + filename, delim_whitespace=True, names=['wave'] + df_age['colname'].tolist())
    return(df_bpass)  # Wave is in Angstroms

    
print "Loading a BPASS2 file"
filename = 'BPASSv2.1_imf135_100/spectra.z020.dat.gz'
df_bpass = load_bpass_spectra(filename)

# Find just one age of that model
age_to_find = 8E6

# Show that the age setup is working
(closest_age) = find_closest_bpass_age(df_age, age_to_find)
print "Closest age was", closest_age, ", so the age column is", closest_age['colname'].values[0]
