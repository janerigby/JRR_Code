from bs4 import BeautifulSoup
import requests
import pandas

specials = [1433, 1424]  # odd numbering for a few proposals.  
regular  = list(range(1549,  2722))  # trying to avoid the commissioning proposals
fullindex = specials.copy()
fullindex.extend(regular)
exclude = [1602, 1628, 1629, 1630, 1631, 1632, 2586]  # These are Commissioning proposals with weird numbers.
newlist = [x for x in fullindex if x not in exclude]

N_gsfc_users = 0
N_proposals_w_gsfc = 0
URLs_w_gsfc = []
outfilename  = "scrape_stsci_JWSTcy1.txt"
wgetfilename = "wget_pdfs_JWSTcy1.txt"
outfile = open(outfilename, "w")
wgetfile = open(wgetfilename, 'w')

for ii in newlist :
    print("DEBUGGING, trying", ii)
    url = 'https://www.stsci.edu/cgi-bin/get-address-info?id=' + str(ii) + '&markupFormat=xml&observatory=JWST' 
    xml_data = requests.get(url).content
    soup = BeautifulSoup(xml_data, 'lxml')
    
    thisone=0
    for entry in soup.findAll('institution'):
        if  ('Goddard') in entry.text  or   ('GSFC') in entry.text :
            N_gsfc_users +=1
            thisone += 1
    if thisone:
        N_proposals_w_gsfc += 1
        URLs_w_gsfc.append(url)
        outfile.write(url + "\n")
        print("   HAS GSFC")
        wgetfile.write("wget https://www.stsci.edu/jwst/phase2-public/" + str(ii) + ".pdf\n")
print("Finished!  There were", N_proposals_w_gsfc, "proposals with GSFC investigators")
print("  and there were", N_gsfc_users, "GSFC proposers.  See outfile", outfilename, "for each URL.")

outfile.close()
wgetfile.close()
