import math
import random
import utils

import logging
log = logging.getLogger(__name__)

class Migration:
    def __init__(self,in_days=[],destination=[]):
        self.in_days=in_days
        self.destination=destination

    def __str__(self):
        s  = 'destination=%s: ' % self.destination
        s += 'in_days=%f '      % self.in_days
        return s

    def update(self,dt):
        if self.in_days:
            self.in_days -= dt

    def migrating(self):
        return True if self.in_days and self.in_days<=0 else False

class MigrationInfo:
    def __init__(self,rates, day):
        self.day = day
        self.destinations=rates.keys()
        self.threshold_day = 10 * 365. # 8.5
        if self.day > self.threshold_day:
            self.relative_rates=utils.accumulate_cdf(rates.values())
            self.total_rate=sum(rates.values())
        else:
            self.relative_rates = 0.
            self.total_rate = 0.

    def __str__(self):
        s  = 'destinations=%s: '  % self.destinations
        s += 'relative rates=%s ' % self.relative_rates
        s += 'total_rate=%f '     % self.total_rate
        return s

    def pick_destination(self):
        return self.destinations[utils.weighted_choice(self.relative_rates)]

    def next_migration(self):
        if not self.total_rate or self.day < self.threshold_day:
            return Migration() # nowhere to migrate
        in_days=random.expovariate(self.total_rate)
        return Migration(in_days, self.pick_destination())

    def destinations_in_timestep(self,n_humans,dt):
        prob=1-math.exp(-self.total_rate*dt)
        migrants=utils.binomialApproxRandom(n_humans,prob)
        log.debug('(humans,prob,migrants)=(%d,%0.2f,%d)',n_humans,prob,migrants)
        return [self.pick_destination() for _ in range(migrants)]
