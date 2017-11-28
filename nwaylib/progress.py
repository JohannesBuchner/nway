import progressbar

class IncrementingProgressBar(progressbar.ProgressBar):
	def increment(self):
		if hasattr(self, 'value'):
			value  = self.value
		else:
			value  = self.currval
		self.update(value + 1)

def bar(ndigits=3, **kwargs):
	if progressbar.__version__ > '3':
		counterfmt = '%(value)'+str(ndigits)+'d'
	else:
		counterfmt = '%'+str(ndigits)+'d'
	
	pbar = IncrementingProgressBar(widgets=[
		progressbar.Percentage(), '|', progressbar.Counter(counterfmt),
		progressbar.Bar(), progressbar.ETA()], **kwargs)
	return pbar

import inspect
import astropy.io.fits as pyfits
if 'clobber' in inspect.getargspec(pyfits.writeto).args:
	arg_overwrite = 'clobber'
else:
	arg_overwrite = 'overwrite'
kwargs_overwrite_true = {arg_overwrite:True}
kwargs_overwrite_false = {arg_overwrite:False}

