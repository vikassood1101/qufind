#!/var/www/html/QUFIND/env/bin python3

import sys
import site

site.addsitedir('/var/www/html/QUFIND/env/lib/python3.9/site-packages')

sys.path.insert(0, '/var/www/html/QUFIND')

from app import app as application