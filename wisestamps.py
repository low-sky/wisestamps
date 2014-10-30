from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits

import numpy as np
import requests
import urllib
t = Table.read('bgps_v2.1.tbl',format='ipac')
t['WISE1']=np.zeros(len(t))+np.nan
t['WISE2']=np.zeros(len(t))+np.nan
t['WISE3']=np.zeros(len(t))+np.nan
t['WISE4']=np.zeros(len(t))+np.nan

correction = np.array([1.9350E-06,2.7048E-06,1.8326e-06,5.2269E-05])
wiseband = np.array(['WISE1','WISE2','WISE3','WISE4'])

for objind in np.arange(5):
	test_object = t[objind]
	queryurl =('http://irsa.ipac.caltech.edu/ibe/search/wise/allwise/p3am_cdd?POS={0},{1}').format(test_object['ra'],test_object['dec'])
	urllib.urlretrieve(queryurl,'votbl.tbl')

	#r = requests.get(queryurl,stream=True)
	#filename = 'votbl.tbl'
	#with open(filename, 'wb') as fd:
	#    for chunk in r.iter_content(chunk_size):
	#        fd.write(chunk)
	#filename.close()
	querytbl = Table.read('votbl.tbl',format='ipac')

	for ctr, bandnum in enumerate([1,2,3,4]):
		params = { 'coadd_id': querytbl['coadd_id'][0],
					   'band': bandnum,
					 }
		params['coaddgrp'] = params['coadd_id'][:2]
		params['coadd_ra'] = params['coadd_id'][:4]
		path = str.format(        	'{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id:s}-w{band:1d}-int-3.fits',
				**params)
		trailing = 	'?center={0},{1}&size=1arcmin'.format(test_object['ra'],test_object['dec'])
		url = 'http://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/' + path+ trailing

		urllib.urlretrieve(url,'stamp.fits.gz')
		img = fits.getdata('stamp.fits.gz')
		w = WCS('stamp.fits.gz')

		x,y = w.wcs_world2pix(test_object['ra'],test_object['dec'],0)
		flux = img[round(y),round(x)]
		print(objind,flux*correction[ctr])
		t[wiseband[ctr]][objind]=flux*correction[ctr]
t.write('bgps_plus_wise.fits',format='fits')