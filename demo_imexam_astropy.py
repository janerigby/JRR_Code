import imexam
viewer = imexam.connect('ds9')
viewer.imexam()

# change defaults
viewer.set_plot_pars('r',"skyrad",20)

