#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

import sys
import warnings

from . import progress

class FakeProgressBar(object):
	def __init__(self, *args):
		pass
	def __call__(self, it):
		return it
	def start(self):
		return self
	def increment(self):
		pass
	def finish(self):
		pass

class NullOutputLogger(object):
	def __init__(self):
		pass
	def log(self, *msg):
		pass
	def warn(self, msg):
		warnings.warn(msg)
	def progress(self, *args, **kwargs):
		return FakeProgressBar()

class NormalLogger(object):
	def __init__(self):
		pass
	def log(self, msg):
		sys.stderr.write('%s\n' % msg)
	def warn(self, msg):
		warnings.warn(msg)
	def progress(self, ndigits=6, *args, **kwargs):
		return progress.bar(ndigits=ndigits, *args, **kwargs)



