#!/var/www/html/qufind/env/bin python3

import sys
import site

site.addsitedir('/var/www/html/qufind/env/lib/python3.6/site-packages')

sys.path.insert(0, '/var/www/html/qufind')

from app import app as application