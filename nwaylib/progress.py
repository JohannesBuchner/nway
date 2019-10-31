import tqdm

def bar(**kwargs):
	return tqdm.tqdm
	
	
import inspect
import astropy.io.fits as pyfits
args = inspect.signature(pyfits.writeto).parameters if hasattr(inspect, 'signature') else inspect.getargspec(pyfits.writeto).args
if 'overwrite' in args:
	arg_overwrite = 'overwrite'
else:
	arg_overwrite = 'clobber'
kwargs_overwrite_true = {arg_overwrite:True}
kwargs_overwrite_false = {arg_overwrite:False}

