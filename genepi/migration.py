import random
import utils

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

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
    def __init__(self,rates):
        self.destinations=rates.keys()
        self.relative_rates=utils.accumulate_cdf(rates.values())
        self.total_rate=sum(rates.values())

    def __str__(self):
        s  = 'destinations=%s: '  % self.destinations
        s += 'relative rates=%s ' % self.relative_rates
        s += 'total_rate=%f '     % self.total_rate
        return s

    def pick_destination(self):
        return self.destinations[utils.weighted_choice(self.relative_rates)]

    def next_migration(self):
        if not self.total_rate:
            return Migration() # nowhere to migrate
        in_days=random.expovariate(self.total_rate)
        return Migration(in_days, self.pick_destination())

    def destinations_in_timestep(self,n_humans,dt):
        n_migrants=utils.poissonRandom(n_humans*self.total_rate*dt)
        return [self.pick_destination() for _ in range(n_migrants)]
