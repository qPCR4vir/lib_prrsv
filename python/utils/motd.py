"""Screwing arround with MOTD
"""

import time
from gatherer import utahAQ as UAQ

class Atmosphere(object):
    """
    """

    def __str__(self):
        f=UAQ.AQFeed()
        data=f.next()
        lines=[
            time.strftime('%A, %B %d, %Y - %I:%M %p',data['date']),
            "\ttemp:%s F\twind: %s %s mph" %(
            data['temperature'],data['wind_dir'],data['wind_speed']),
            "\tozone: %s ppm\t(8hr avg: %s)\n\tPM2.5: %s ug/m^3\t(24 hr avg: %s)" %(
            data['ozone'],data['ozone_8hr_avg'],data['pm25'],data['pm25_24hr_avg'] ),
            "\tco: %s\tno2: %s\tnox: %s" % (data['co'],data['no2'],data['nox'])]
        
        return '\n'.join(lines)
        
        
