"""Accepted-step metadata: steps.csv starts with step 1, not an initial state."""
import math

def accepted_rows(records):
 return [row for row in records if int(row['step']) > 0]

def statistics(records):
 rows=accepted_rows(records)
 dt=[float(row['dt_s']) for row in rows]
 if any(not math.isfinite(x) or x<=0 for x in dt):
  raise ValueError('Accepted steps require finite positive time increments')
 total=math.fsum(dt)
 return dict(min=min(dt) if dt else None,max=max(dt) if dt else None,
             mean=total/len(dt) if dt else None,accepted=len(dt),sum_duration_s=total)
