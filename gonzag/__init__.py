"""
"""

#from geopy.geocoders import *  # noqa
#from geopy.location import Location  # noqa
#from geopy.point import Point  # noqa
#from geopy.timezone import Timezone  # noqa
#from geopy.util import __version__  # noqa

# geopy.geocoders.options must not be importable as `geopy.options`,
# because that is ambiguous (which options are that).
#del options  # noqa

# `__all__` is intentionally not defined in order to not duplicate
# the same list of geocoders as in `geopy.geocoders` package.


"""
Constants that are defined for gonzag
-------------------------------------
"""
# EARTH CONSTANTS

#: Volumetric mean radius (km)
VOLUMETRIC_MEAN_RADIUS: float = 6371.0

#: Convert degree to km
DEG_2_KM: float = 111.11

#: CELERITY
CELERITY: float = 299800000.0

#: Baseline (m)
BASELINE = 10

